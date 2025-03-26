D1vsD2_morphometric_comparison.R: perform exploration, analysis and visualziation of D1- and D2-MSNs in the P56 datasets, limiting to CP neurons.

Use the dataset that has 'surface neurons' excluded. Surface neurons are defined as located within the 30% area from the tissue cut surface. 

Generated files:
Plots:
- 1: Basic frequency distribution tables (Type, Sex, Striatal.Subregion, Brain..)
- 2: Distribution plots of Type in Striatal.Subregion, Straital.Community, Sex, Sex in Striatal.Subregion
- Type_distribution_in_community.pdf: Bar plot showing D1- and D2-MSN distribution across Striatal.Community
- 4: Individual histogram distribution of each of the 31 morphometrics
- 5: Correlation heatmap of 31 morphometrics
- 6: Hierarchical Clustering of 31 morphometris
- 7: PCA plot showing Type, Striatal.Subregion, Striatal.Community
- 7.1: PCA plot showing Type, Striatal.Subregion, Striatal.Community, highlighting one level at a time
- 8: UMAP plot showing Type, Striatal.Subregion, Striatal.Community
- 8.1: UMAP plot showing Type, Striatal.Subregion, Striatal.Community, highlighting one level at a time
- 9: Heatmap of normalized feature distribution across Striatal.Community
- 9.1: Heatmap of normalized feature distribution overall, across Striatal.Subregion and Striatal.Community
- 10: Heatmap of normalized feature comparison for pairwise Striatal.Subregions (CPr vs. CPi, CPr vs. CPc, CPi vs. CPc), by individual brain
- 11: Heatmap of normalized feature comparison for pairwise Striatal.Subregions (CPr vs. CPi, CPr vs. CPc, CPi vs. CPc), by D1 and D2
- 11.1: Heatmap of normalized feature comparison for pairwise Striatal.Subregions (CPr vs. CPi, CPr vs. CPc, CPi vs. CPc), by D1 and D2, only profound features
- 12: Heatmap of normalized feature comparison for pairwise Striatal.Subregions (CPr vs. CPi, CPr vs. CPc, CPi vs. CPc), by D1 and D2 in each brain
- 13: Heatmap of mean differences of morphometrics between D2 vs. D1 overall and across Striatal.Subregion
- 14: Heatmap of mean differences of morphometrics between D2 vs. D1 overall, across Striatal.Subregion and across Striatal.Community
- 15: Heatmap of percent change of morphometrics D2 vs. D1 overall, across Striatal.Subregion and across Striatal.Community
- 16: Heatmap of mean differences of morphometrics between pairwise Striatal.Subregions, by D1 and D2
- 17: Heatmap of morphometric distribution overall, across Striatal.Subregion and across Striatal.Community
- 18: Box plot showing morphometric distribution across Striatal.Subregion, for all neurons, D1-MSNs, and D2-MSNs for each morphometric
- 19: Box plot showing morphometric comparison D1 vs. D2, for all neurons and across Striatal.Subregion for each morphometric
