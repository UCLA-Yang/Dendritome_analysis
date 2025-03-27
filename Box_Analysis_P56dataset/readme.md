boxes_analysis.R: Box analysis of all neurons from P56 dataset

This script will generate (the numbers correspond to the result files): 
- 1: Bar plot showing the distribution of cell count per box by each brain 
- 2: Bar plot showing distribution of cell count per box for all brains combined
- 3: Dendrogram and similarity heatmap for different parameter combinations in the box clustering
- 4: The plot showing the changes of Silouette score for different parameter combinations in the box clustering
- 5: Dendrogram and similarity heatmap in the final box clustering
- 6: Bar plot showing cell count per Dendritic Module (DM)
- 7: Distribution of sizes of connected components in each DM
- 7.1: Heatmap showing distribution of boxes across the Striatal.Community and Morphological Territories (MTs)
- 7.2: Bar plot showing cell count per Morphological Territories (MTs)
- 8: Line plot showing the covariation of morphometrics of top 25 representative cells in each MT
- 9: Line plot showing the covariation of morphometrics of top 10 representative cells in each MT
- 10: Heatmap showing correlation heatmap of top 25 representative cells by each MT
- 11: Heatmap showing the distribution of median-normalized features across DMs, by each brain and all brains combined
- 12: Bar plot showing the distribution of D1- and D2-MSNs in each Dendritic Module (DM)
- 13: Bar plot showing the distribution of D1- and D2-MSNs in each Morphological Territories (MTs)
- 14: Heatmap showing D1 vs. D2 morphometric comparison in each box for each DM (boxes are filtered based on different cell count thresholds)
- 15: Heatmap showing the distribution of average morphometrics across DMs
- 16: Heatmap showing the distribution of average morphometrics across DMs and individual brain
- 17: Scatter plots showing box membership and cell count per box, for all DMs combined
- 18: Scatter plots showing box membership and cell count per box, for each DM
- 19: Heatmap showing the box membership in every DM 

Other .csv data files: box_cluster_bicor_D1D2.csv, box_cluster_bicor_D1D2_with_xyz.csv, box_cluster_bicor_D1D2_with_xyz_filtered_thresh7.csv, cell_rank_in_box_cluster_top25.csv, morpho_brain_adjusted_box_cluster_xyz_D1D2.csv, morpho_brain_un_ajusted_box_cluster_xyz_D1D2.csv, chisq_test_DM_component_distribution.csv
