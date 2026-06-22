#' ---
#' title: "MLP Beef Flavor"
#' author: "Anurag Banerjee"
#' date: "2026-04-05"
#' output: pdf_document
#' ---
#' 
## ----setup, include=FALSE---------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' 
#' 
## ----echo=FALSE, include=FALSE----------------------------------------------------------------------------------
library(tidyverse)
library(nnet)
library(caret)
library(vip)
library(ggplot2)

#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
df <- read_excel("Feature List.xlsx")
head(df)

#' 
## ---------------------------------------------------------------------------------------------------------------
outcome_col <- "Liking"
chem_cols <- colnames(df)[4:505]
dat <- df[, c(outcome_col, chem_cols)]
dat <- na.omit(dat)
dat[, chem_cols] <- log(dat[, chem_cols])
names(dat)[2:503] <- paste0("chemical_", 1:502)
chem_cols <- colnames(dat)[2:503]

# Safety checks
stopifnot(nrow(dat) >= 8)
stopifnot(length(chem_cols) >= 5)

#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
#Center and Scale the chemical variables
X <- dat[, chem_cols, drop = FALSE]
y <- dat[[outcome_col]]

X_scaled <- scale(X, center = T, scale = T) %>% as.data.frame()

dat_scaled <- bind_cols(
  data.frame(y = y),
  X_scaled
)

#' 
## ---------------------------------------------------------------------------------------------------------------
rmse_fun <- function(obs, pred) {
  sqrt(mean((obs - pred)^2))
}

r2_fun <- function(obs, pred) {
  if (sd(pred) == 0 || sd(obs) == 0) return(NA_real_)
  cor(obs, pred)^2
}

#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
#Fitting a conservative small MLP

fit_mlp <- function(train_df, seed, size = 3, decay = 0.2, maxit = 500) {
  set.seed(seed)
  nnet(
    y ~ .,
    data = train_df,
    size = size,
    linout = TRUE,
    decay = decay,
    maxit = maxit,
    trace = FALSE
  )
}

screen_train_variance <- function(train_df, feature_names, top_k = 40) {
  vars <- sapply(train_df[, feature_names, drop = FALSE], var, na.rm = TRUE)
  names(sort(vars, decreasing = TRUE))[1:min(top_k, length(vars))]
}

screen_train_dcorr <- function(train_df, feature_names, top_k = 40) {
  dcors <- sapply(feature_names, function(v) {
    x <- train_df[[v]]
    y <- train_df[["y"]]
    ok <- complete.cases(x, y)
    
    if (sum(ok) < 3) return(NA_real_)
    if (sd(x[ok]) == 0) return(NA_real_)
    
    energy::dcor(x[ok], y[ok])
  })
  
  dcors <- sort(dcors, decreasing = TRUE, na.last = NA)
  names(dcors)[1:min(top_k, length(dcors))]
}


#' 
#' 
#' Permutation Importance:
## ---------------------------------------------------------------------------------------------------------------
# Importance = increase in RMSE when one feature is permuted
perm_importance_fold <- function(model, train_df, test_df, selected_features,
                                 all_features, n_repeats = 10, seed = 1) {
  set.seed(seed)

  base_pred <- predict(model, newdata = test_df)
  base_rmse <- rmse_fun(test_df$y, base_pred)

  out <- map_dfr(selected_features, function(v) {
    delta <- numeric(n_repeats)

    for (b in seq_len(n_repeats)) {
      tmp <- test_df
      tmp[[v]] <- sample(tmp[[v]], replace = FALSE)#Randomly shuffle each observation for a particular chemical 
      p <- predict(model, newdata = tmp)
      delta[b] <- rmse_fun(tmp$y, p) - base_rmse #Change in RMSE after breaking the relationship with a chemical variable, this quantifies how important the particular chemical variable
    }

    tibble(
      feature = v,
      importance = mean(delta, na.rm = TRUE)
    )
  })

  # Add unselected features as zero/NA importance explicitly if desired
  missing_feats <- setdiff(all_features, selected_features)
  if (length(missing_feats) > 0) {
    out <- bind_rows(
      out,
      tibble(feature = missing_feats, importance = 0)
    )
  }

  out
}



#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
run_cv_seed <- function(dat_scaled, all_features, seed_id,
                        v = 5, repeats = 5, top_k = 40,
                        size = 3, decay = 0.2, maxit = 500,
                        perm_repeats = 10) {
  #seed_id <- 1
  set.seed(seed_id)
  #v =5; repeats = 5;
  folds <- createMultiFolds(dat_scaled$y, k = v, times = repeats)

  pred_store <- vector("list", length(folds))
  imp_store  <- vector("list", length(folds))

  fold_names <- names(folds)

  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(dat_scaled)), unique(train_idx))

    train_df <- dat_scaled[train_idx, , drop = FALSE]
    test_df  <- dat_scaled[test_idx, , drop = FALSE]

    #selected <- screen_train_variance(train_df, all_features, top_k = top_k)
    selected <- screen_train_dcorr(train_df, all_features, top_k = top_k)
    
    train_sub <- train_df[, c("y", selected), drop = FALSE]
    test_sub  <- test_df[, c("y", selected), drop = FALSE]

    model <- fit_mlp(
      train_df = train_sub,
      seed = seed_id + i,
      size = size,
      decay = decay,
      maxit = maxit
    )

    pred <- predict(model, newdata = test_sub)

    pred_store[[i]] <- tibble(
      obs = test_sub$y,
      pred = as.numeric(pred),
      fold = fold_names[i],
      seed = seed_id
    )

    imp_store[[i]] <- perm_importance_fold(
      model = model,
      train_df = train_sub,
      test_df = test_sub,
      selected_features = selected,
      all_features = all_features,
      n_repeats = perm_repeats,
      seed = seed_id + 1000 + i
    ) %>%
      mutate(
        fold = fold_names[i],
        seed = seed_id
      )
  }

  pred_df <- bind_rows(pred_store)
  imp_df  <- bind_rows(imp_store)

  perf_df <- pred_df %>%
    summarise(
      seed = first(seed),
      rmse = rmse_fun(obs, pred),
      r2 = r2_fun(obs, pred)
    )

  imp_summary <- imp_df %>%
    group_by(seed, feature) %>%
    summarise(
      mean_importance = mean(importance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(seed, desc(mean_importance)) %>%
    group_by(seed) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  list(
    performance = perf_df,
    predictions = pred_df,
    importance = imp_summary
  )
}



#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
# Jaccard overlap between top-k feature sets
jaccard_topk <- function(df1, df2, k = 20) {
  s1 <- df1 %>% arrange(rank) %>% slice_head(n = k) %>% pull(feature)
  s2 <- df2 %>% arrange(rank) %>% slice_head(n = k) %>% pull(feature)
  length(intersect(s1, s2)) / length(union(s1, s2))
}


#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
#Real Data - seed instability
seed_grid <- 1:30

real_runs <- purrr::map(seed_grid, ~ run_cv_seed(
  dat_scaled = dat_scaled,
  all_features = chem_cols,
  seed_id = .x,
  v = 5,
  repeats = 5,
  top_k = 50,
  size = 3,
  decay = 0.2,
  maxit = 500,
  perm_repeats = 10
))

real_perf <- map_dfr(real_runs, "performance")
real_preds <- map_dfr(real_runs, "predictions")
real_imp <- map_dfr(real_runs, "importance")


#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
# Top-20 selection frequency across seeds
real_top20_freq <- real_imp %>%
  group_by(seed) %>%
  slice_min(rank, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  count(feature, name = "top20_count") %>%
  mutate(top20_freq = top20_count / length(seed_grid)) %>%
  arrange(desc(top20_freq))

#' 
#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
imp_by_seed <- split(real_imp, real_imp$seed)
seed_pairs <- combn(names(imp_by_seed), 2, simplify = FALSE)

jaccard_real <- map_dbl(seed_pairs, function(pair) {
  jaccard_topk(imp_by_seed[[pair[1]]], imp_by_seed[[pair[2]]], k = 20)
})

real_jaccard_summary <- tibble(
  median_jaccard_top20 = median(jaccard_real),
  mean_jaccard_top20 = mean(jaccard_real),
  min_jaccard_top20 = min(jaccard_real),
  max_jaccard_top20 = max(jaccard_real)
)

#' 
#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
real_imp_overall <- real_imp %>%
  group_by(feature) %>%
  summarise(
    mean_importance = mean(mean_importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))

#' 
#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
run_one_permutation <- function(b, dat_scaled, feature_names) {
  set.seed(10000 + b)

  dat_perm <- dat_scaled
  dat_perm$y <- sample(dat_perm$y, replace = FALSE)

  perm_runs <- purrr::map(1:10, ~ run_cv_seed(
    dat_scaled = dat_perm,
    all_features = feature_names,
    seed = 1000 * b + .x,
    v = 5,
    repeats = 3,
    top_k = 40,
    size = 3,
    decay = 0.2,
    maxit = 500,
    perm_repeats = 5
  ))

  perf <- map_dfr(perm_runs, "performance")
  imp  <- map_dfr(perm_runs, "importance")

  top_imp <- imp %>%
    group_by(feature) %>%
    summarise(m = mean(mean_importance, na.rm = TRUE), .groups = "drop") %>%
    slice_max(m, n = 1, with_ties = FALSE)

  tibble(
    perm_id = b,
    mean_rmse = mean(perf$rmse, na.rm = TRUE),
    mean_r2 = mean(perf$r2, na.rm = TRUE),
    max_r2 = max(perf$r2, na.rm = TRUE),
    top_feature = top_imp$feature,
    top_importance = top_imp$m
  )
}

perm_results <- map_dfr(1:100, ~ run_one_permutation(.x, dat_scaled, chem_cols))


#' 
#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
real_summary <- real_perf %>%
  summarise(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse   = sd(rmse, na.rm = TRUE),
    mean_r2   = mean(r2, na.rm = TRUE),
    sd_r2     = sd(r2, na.rm = TRUE),
    max_r2    = max(r2, na.rm = TRUE)
  )

real_top_feature <- real_imp_overall %>% slice_max(mean_importance, n = 1, with_ties = FALSE)

pval_max_r2 <- mean(perm_results$max_r2 >= real_summary$max_r2, na.rm = TRUE)
pval_top_imp <- mean(perm_results$top_importance >= real_top_feature$mean_importance, na.rm = TRUE)

summary_table <- tibble(
  metric = c(
    "Real mean RMSE",
    "Real SD RMSE",
    "Real mean R2",
    "Real SD R2",
    "Real max R2",
    "Median Jaccard top20",
    "Mean Jaccard top20",
    "P(null max R2 >= real max R2)",
    "P(null top importance >= real top importance)"
  ),
  value = c(
    real_summary$mean_rmse,
    real_summary$sd_rmse,
    real_summary$mean_r2,
    real_summary$sd_r2,
    real_summary$max_r2,
    real_jaccard_summary$median_jaccard_top20,
    real_jaccard_summary$mean_jaccard_top20,
    pval_max_r2,
    pval_top_imp
  )
)

#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
write.csv(real_perf, "real_seed_performance_39x502.csv", row.names = FALSE)
write.csv(real_preds, "real_seed_predictions_39x502.csv", row.names = FALSE)
write.csv(real_imp, "real_seed_importance_39x502.csv", row.names = FALSE)
write.csv(real_imp_overall, "real_overall_importance_39x502.csv", row.names = FALSE)
write.csv(real_top20_freq, "real_top20_selection_frequency_39x502.csv", row.names = FALSE)
write.csv(real_jaccard_summary, "real_jaccard_summary_39x502.csv", row.names = FALSE)
write.csv(perm_results, "permutation_null_results_39x502.csv", row.names = FALSE)
write.csv(summary_table, "mlp_instability_summary_39x502.csv", row.names = FALSE)



#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
png("r2_real_vs_null_39x502.png", width = 900, height = 600)
hist(
  perm_results$max_r2,
  breaks = 25,
  main = "Null distribution of max R2 under permuted outcomes",
  xlab = "Max R2 under permutation",
  col = "grey80",
  border = "white"
)
abline(v = real_summary$max_r2, col = "red", lwd = 3)
dev.off()

#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
png("top20_freq_39x502.png", width = 1000, height = 700)
real_top20_freq %>%
  slice_max(top20_freq, n = 25) %>%
  ggplot(aes(x = reorder(feature, top20_freq), y = top20_freq)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "How often each chemical appears in MLP top-20 across seeds",
    x = "Chemical",
    y = "Selection frequency"
  )
dev.off()

#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
png("importance_top25_39x502.png", width = 1000, height = 700)
real_imp_overall %>%
  slice_max(mean_importance, n = 25) %>%
  ggplot(aes(x = reorder(feature, mean_importance), y = mean_importance)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  labs(
    title = "Average permutation importance of top 25 chemicals",
    x = "Chemical",
    y = "Mean CV permutation importance"
  )
dev.off()

#' 
#' 
#' Methodology:
#' 1. In order to compute importance in each run for each feature, we permute a single feature and see how much rmse increased because of this permutation of the feature
#' 2. To check whether similar features are being selected we look at the overlap frequence and the jaccard scoring (I'm not sure what the metric is called it's jaccard overlap but it's basically count(A intersection B)/count(A union B)), it measures the overlap)
#' 3. Then in order to see whether the model is picking up signals we permute the response and fit models on what is essentially noise, comparison between this model and the real model will tell us whether we are just picking up noise or is there some signal.
#' 4. First we compare the error metrics, R^2, RMSE
#' 5. Then we compare the feature importance between the real data runs v/s the null distribution runs, different features with more importance for real data should happen
#' 6. Then if the real data passes prior checks we will check if usually same features are being identified as the important ones (this is most likely where the model will fail)
#' 
#' Compare the error metrics from the model(s) using real data v/s the permuted data (in this we jumbled the response, so that this is all noise)
## ---------------------------------------------------------------------------------------------------------------
library(dplyr)

rmse_compare <- bind_rows(
  real_perf %>%
    transmute(rmse = rmse, source = "Real data"),
  perm_results %>%
    transmute(rmse = mean_rmse, source = "Permuted outcome")
)

ggplot(rmse_compare, aes(x = rmse, fill = source, color = source)) +
  geom_density(alpha = 0.25, adjust = 1.2) +
  labs(
    title = "RMSE distribution: real data vs permuted outcome",
    x = "RMSE",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

#' 
## ---------------------------------------------------------------------------------------------------------------
real_mean_rmse <- mean(real_perf$rmse, na.rm = TRUE)
real_min_rmse  <- min(real_perf$rmse, na.rm = TRUE)

ggplot(perm_results, aes(x = mean_rmse)) +
  geom_histogram(bins = 30, fill = "grey75", color = "white") +
  geom_vline(xintercept = real_mean_rmse, color = "red", linewidth = 1.2) +
  geom_vline(xintercept = real_min_rmse, color = "blue", linewidth = 1.2, linetype = "dashed") +
  labs(
    title = "Null RMSE distribution under permuted outcomes",
    subtitle = "Red = mean RMSE from real runs, Blue dashed = best real RMSE",
    x = "Mean RMSE under permutation",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

#' 
#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
r2_compare <- bind_rows(
  real_perf %>%
    transmute(r2 = r2, source = "Real data"),
  perm_results %>%
    transmute(r2 = mean_r2, source = "Permuted outcome")
)

ggplot(r2_compare, aes(x = r2, fill = source, color = source)) +
  geom_density(alpha = 0.25, adjust = 1.2) +
  labs(
    title = "R-squared distribution: real data vs permuted outcome",
    x = "R^2",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)


#' 
## ---------------------------------------------------------------------------------------------------------------
real_mean_r2 <- mean(real_perf$r2, na.rm = TRUE)
real_max_r2  <- max(real_perf$r2, na.rm = TRUE)

ggplot(perm_results, aes(x = mean_r2)) +
  geom_histogram(bins = 30, fill = "grey75", color = "white") +
  geom_vline(xintercept = real_mean_r2, color = "red", linewidth = 1.2) +
  geom_vline(xintercept = real_max_r2, color = "blue", linewidth = 1.2, linetype = "dashed") +
  labs(
    title = "Null R-squared distribution under permuted outcomes",
    subtitle = "Red = mean R-squared from real runs, Blue dashed = best real R-squared",
    x = "Mean R-squared under permutation",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

#' 
#' 
#' Seems like overall the real data does have some signal that is being captured by the MLP, the model trained on real data is clearly better, however there is some data leakage, I standardized the data using entire dataset and didn't include it in the model fitting pipeline.
#' 
## ---------------------------------------------------------------------------------------------------------------
real_top_imp <- real_imp_overall %>%
  slice_max(mean_importance, n = 1, with_ties = FALSE) %>%
  pull(mean_importance)

ggplot(perm_results, aes(x = top_importance)) +
  geom_histogram(bins = 30, fill = "grey75", color = "white") +
  geom_vline(xintercept = real_top_imp, color = "red", linewidth = 1.2) +
  labs(
    title = "Null distribution of top feature importance",
    subtitle = "Red line = top observed importance from real data",
    x = "Top importance under permuted outcomes",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

#' 
#' Fraction of null runs whose top importance is at least as large as the observed real top importance:
#' 
## ---------------------------------------------------------------------------------------------------------------
p_top_imp <- mean(perm_results$top_importance >= real_top_imp, na.rm = TRUE)
p_top_imp

#' 
## ---------------------------------------------------------------------------------------------------------------
ggplot(perm_results, aes(x = top_importance)) +
  geom_density(fill = "grey80", alpha = 0.6) +
  geom_vline(xintercept = real_top_imp, color = "red", linewidth = 1.2) +
  labs(
    title = "Observed top importance vs null distribution",
    x = "Top importance",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)


#' 
## ---------------------------------------------------------------------------------------------------------------
real_top10_imp <- real_imp_overall %>%
  slice_max(mean_importance, n = 10, with_ties = FALSE)

ggplot(perm_results, aes(x = top_importance)) +
  geom_histogram(bins = 30, fill = "grey80", color = "white") +
  geom_vline(
    data = real_top10_imp,
    aes(xintercept = mean_importance),
    color = "red",
    alpha = 0.5
  ) +
  labs(
    title = "Top 10 real feature importances vs null top-importance distribution",
    x = "Importance",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

#' It seems like other than the top feature the others have importances comparable or lower than the one captured through the null distribution which is problematic
#' 
## ---------------------------------------------------------------------------------------------------------------
null_95 <- quantile(perm_results$top_importance, 0.95, na.rm = TRUE)

real_top20_imp <- real_imp_overall %>%
  slice_max(mean_importance, n = 20, with_ties = FALSE) %>%
  mutate(feature = reorder(feature, mean_importance))

ggplot(real_top20_imp, aes(x = feature, y = mean_importance)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = null_95, color = "red", linewidth = 1.2, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Top 20 real feature importances",
    subtitle = "Dashed red line = 95th percentile of null top importance",
    x = "Feature",
    y = "Mean importance"
  ) +
  theme_minimal(base_size = 14)


#' 
#' A permutation null benchmarks asks whether there is a relationship in the real data by comparing it to a dataset where there is no relationship, as we can observe most of the importance scores from the real data overlaps with that of the null distribution, so the model is picking up signal but it is not producing evidence at the feature level that the model on real data is significantly better at identifying predictive features than model on null distribution. Importance ranking is weak relative to the null benchmark.
#' 
#' At most, there is limited evidence for one particularly influential feature, whereas the remaining ranked chemicals are not well distinguished from null importance levels.
#' 
#' 
#' Next we will look at the feature selection frequency into the top 20 features:
#' 
## ---------------------------------------------------------------------------------------------------------------
top_k <- 20

real_topk_freq <- real_imp %>%
  group_by(seed) %>%
  slice_min(rank, n = top_k, with_ties = FALSE) %>%
  ungroup() %>%
  count(feature, name = "count") %>%
  mutate(freq = count / n_distinct(real_imp$seed)) %>%
  arrange(desc(freq))


#' 
## ---------------------------------------------------------------------------------------------------------------
real_topk_freq %>%
  slice_max(freq, n = 30) %>%
  mutate(feature = reorder(feature, freq)) %>%
  ggplot(aes(x = feature, y = freq)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = paste("Selection frequency in top", top_k, "across runs"),
    x = "Feature",
    y = "Selection frequency"
  ) +
  theme_minimal(base_size = 14)

#' 
#' 
#' Does the model pick similar top k features in each run:
#' J(A, B) = |A and B|/|A or B|
#' 
## ---------------------------------------------------------------------------------------------------------------
imp_by_seed <- split(real_imp, real_imp$seed)
seed_pairs <- combn(names(imp_by_seed), 2, simplify = FALSE)

jaccard_real <- purrr::map_dbl(seed_pairs, function(pair) {
  jaccard_topk(imp_by_seed[[pair[1]]], imp_by_seed[[pair[2]]], k = top_k)
})

real_jaccard_summary <- tibble(
  mean_jaccard = mean(jaccard_real),
  median_jaccard = median(jaccard_real),
  min_jaccard = min(jaccard_real),
  max_jaccard = max(jaccard_real)
)

#' 
## ---------------------------------------------------------------------------------------------------------------
tibble(jaccard = jaccard_real) %>%
  ggplot(aes(x = jaccard)) +
  geom_histogram(bins = 25, fill = "grey70", color = "white") +
  geom_vline(xintercept = median(jaccard_real), color = "red", linewidth = 1.2) +
  labs(
    title = paste("Pairwise Jaccard overlap of top", top_k, "features"),
    subtitle = "Red line = median Jaccard",
    x = "Jaccard overlap",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

#' 
#' The top 20 features being selected across different runs are changing quite a lot, there's only a modest overlap with a median jaccard overlap of ~0.17, i.e. only about 6 shared features across different runs. It's possible only 6-8 chemicals that we saw with high empirical frequencies are the predictive ones, but I'm not too sure about that, it's possible that this model is varying across runs when it comes to identifying predictive features.
#' 
#' 
## ---------------------------------------------------------------------------------------------------------------
top_features <- real_topk_freq %>%
  slice_max(freq, n = 30) %>%
  pull(feature)

membership_df <- real_imp %>%
  group_by(seed) %>%
  slice_min(rank, n = top_k, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(selected = 1) %>%
  filter(feature %in% top_features) %>%
  dplyr::select(seed, feature, selected)

all_grid <- expand.grid(
  seed = unique(real_imp$seed),
  feature = top_features
)

membership_plot_df <- all_grid %>%
  left_join(membership_df, by = c("seed", "feature")) %>%
  mutate(selected = ifelse(is.na(selected), 0, selected))

ggplot(membership_plot_df, aes(x = feature, y = factor(seed), fill = factor(selected))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "grey90", "1" = "steelblue")) +
  labs(
    title = paste("Top", top_k, "selection membership across runs"),
    x = "Feature",
    y = "Seed",
    fill = "Selected"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#' 
#' 
