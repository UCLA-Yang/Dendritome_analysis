##### lineplot_compare_regional_morphologies.R: <br>To generate line plots that describe the regional morphological distribution/differences for P56 WT MSNs, 12m WT MSNs and 12m HD MSNs #####

The code will generate:
1) line-plot_averaged_normalized_feature_subregion_pairwise_P56_12m.pdf: the line plot showing pairwise regional differences (CPr vs. CPi, CPr vs. CPc and CPi vs. CPc) for P56 WT MSNs, 12m WT MSNs and 12m HD MSNs by D1-MSNs and D2-MSNs
2) (Fig8d) lineplot_regional_differences_3sets.pdf: the line plot showing the distribution of median-normalized features across CPr, CPi and CPc for P56 WT MSNs, 12m WT MSNs and 12m HD MSNs by D1-MSNs and D2-MSNs

and other raw statistics for annotations

##### explorartory_analysis_surface_neurons.ipynb: <br>To perform exploratory analysis of the surface neurons in the P56 dataset #####

The code will generate plots (the numbers correspond to the figure numbers): 
1) (Fig3d). Line plot showing cumulative cell count by standardized Z location of the neurons within the tissue section
2) (Suppl Fig3). Average normalized values by standardized Z locations for each of the 31 features
3) (Fig3e). The bar plot showing the averaged normalized value of features at 10%, 20%, 30%, 40% area from the cut surface
4) The bar plot showing the proportion of surface neurons in D1- and D2-MSNs
5) The bar plot showing the proportion of surface neurons in striatal subregions (CPr, CPi and CPc)
6) (Suppl Fig8b). The bar plot showing the proportion of surface neurons in 7 dendritic modules (DMs)

##### additional_analysis_for_box_clustering.R: <br>To perform additional analysis for box-based clustering #####

The code will generate plots (the numbers correspond to the figure numbers): 
1) Dendrogram and similarity heatmap of box clustering based on different parameter combinations (excluded boxes with < 5 neurons)
2) The silouette score line plot based on different parameter combinations (excluded boxes with < 5 neurons)
3) (Suppl Fig8d). The final dendrogram and heatmap of box clustering (excluded boxes with < 5 neurons)
4) (Suppl Fig8c). The scatter plot showing feature CV by cell count per box (all boxes)
5) The box plot showing feature CV by categorized cell count per box (all boxes)
6) (Suppl Fig8a). Bar plot showing the proportion of brains across 7 DMs
7) Bar plot showing Cohen's D of differences between DM2 and DM6
8) Box plot showing feature differences between DM2 and DM6
9) LDA projection of DM2 and DM6
