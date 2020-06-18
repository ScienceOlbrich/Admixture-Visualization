# Admixture-Visualization
Some functions to facilitate plotting the output of PLINK in R.

This includes three segements pertaining to "ROH" (Runs of homozygosity), admixture analyses (barplots and pie-charts on a global map) and "PCA" (principal component analyses in 2/3D)

### ROH

Boxplot of chosen populations showing either sub-groups: short, medium, long (that can be selected via threshold) or the sum of all runs per individual larger than a chosen threshold

### Admixture analyses

Barplots of admixture components per sample, facet wrapped by population or dataset. Sorted by significance of clusters based on a main population to compare with. 

Pie charts of admixture components depicted at their coordinates on a world map. This code is based on the tutorial/git of https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial and adapted to work with our requirements.


### PCA

Three functions to produce interactive an static 3D-plots as well as a set of 2D-plots for the combinations of axes 1-4.

