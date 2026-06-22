# OMICS_Beef_Flavor
Objective: Find the chemicals that are indicative of liking beef flavor and reinject them into beef sample.

This was a project from when I worked at Statistical Consulting Service at the Ohio State University. All of the statistical and data analysis was done by me and I took into consideration the regular feedback from my client.

There were 10 unique beef samples which were given to test subjects and they rated each sample on a hedonic scale of 1-9.

The objective was to reliably identify the top 20 chemicals out of the 500 which were predictive of liking beef flavor, due to the small n large p constraint this proved especially difficult. Several approaches were considered and tested, such as multi-layer perceptrons and Kernel PLS. However, it was not possible to get reliable results, I found that there was a signal when it came to prediction but not for feature (chemical) selection based on sampling tests. In order to test the same we used permutation tests, where we broke the signal and additionally created a null dataset and compared results with the original dataset. 

The rmd files are the codes for the 2 models mentioned above, and rest are figures and tables to portray the results. We found that the null distribution, at times, had higher feature importance. The file titled Beef_Sample_Report is the final report. 