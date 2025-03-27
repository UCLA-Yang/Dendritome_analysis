# work directory
setwd("C:/Users/yanyanming77/Desktop/Ming_Transition/Morphometric_Analysis/Additional_analysis_for_paper-Mar.25.2025/additional_for_response_box_anlaysis")

# source functions
func_dir = 'C:/Users/yanyanming77/Desktop/Ming_Transition/Morphometric_Analysis/Functions/'

source(paste0(func_dir, "individualAnalysis-General-010.R"));
source(paste0(func_dir, "GNVFunctions-018-02.R"));
source(paste0(func_dir, "networkFunctions-extras-20.R"));
source(paste0(func_dir, "labelPoints2-01.R"));
source(paste0(func_dir, "heatmap.wg.R"));

# import libraries
library(anRichmentMethods)
library(Cairo)
library(openxlsx)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(ggplotify)
library(gplots)
library(RColorBrewer)
library(ggsignif)
library(ggsci)
library(car)
library(ggtext)
library(tidyr)
library(circlize)
library(pheatmap)
library(ComplexHeatmap)
library(fossil)
library(igraph)
library(reshape2)
library(scales)
library(effsize)
library(MASS)
library(fossil)
library(mclust)

###############################################################################################
############################ read data ########################################################
###############################################################################################
data_dir = 'TME_morpho_box_May-20-2024_2466neurons.csv'
stats1 = read.csv(data_dir,check.names = FALSE)
colnames(stats1)
dim(stats1) # 2466
table(stats1$Type) # 1341 D1, 1125 D2

###############################################################################################
############################ Rename, Add columns, Process data ################################
###############################################################################################
# rename all variables
names(stats1)[names(stats1) == "box"] = "Box"

# remove columns
stats1 = subset(stats1, select = -c(Diameter, Fragmentation, Surface, Volume, Pk_classic, Contraction, Helix))
ncol(stats1) # 41

stats = list(stats1)

# create a dataset with all numeric statistics
firstStat = "Bif_ampl_local";
lastStat = "Convexity"
morf0 = lapply(stats, function(stats) setRownames(stats[match(firstStat, names(stats)):match(lastStat, names(stats))], 
                                                  make.unique(spaste(stats$`Brain`, ".", stats$`Reconstruction #`))));
numStats = lapply(morf0, dropConstantColumns)[[1]]

numericCols = lapply(list(numStats), colnames)[[1]]
nNumStats = sapply(list(numStats), ncol) # 31

# feature clustering
numStats.scaled = lapply(list(numStats), scale)[[1]]
tree = lapply(list(numStats), function(x) hclust(as.dist(1-bicor(x, use = 'p', maxPOutliers = 0.05)), method = "average"))
height = 0.2;
clusters = lapply(tree, cutree, h = height); 

# define levels of cell type and regions
order_subregion = c('CPr','CPi','CPc')
order_type = c('D1', 'D2')

stats1 = within(stats1, {
  Striatal.Subregion = factor(Striatal.Subregion, levels = order_subregion)
  Type = factor(Type, levels = order_type)
  Box = as.factor(Box)
})

stats[[1]] = within(stats[[1]], {
  Striatal.Subregion = factor(Striatal.Subregion, levels = order_subregion)
  Type = factor(Type, levels = order_type)
})

# find non-redundant features 
a = as.data.frame(clusters[[1]], stringsAsFactors = FALSE)
a$Feature = rownames(a)
colnames(a)[1] = 'group'
non_redundant_features = a %>%
  arrange(group) %>%
  group_by(group) %>%
  slice(1) %>%
  pull(Feature)
length(non_redundant_features) # 21 non-redundant features

####################################################################
################## Exclude box with < 5 neurons #################
####################################################################

cell_counts = stats1 %>%
  group_by(Box) %>%
  summarise(count = n(), .groups = 'drop') %>% data.frame()
boxes_lt5 = cell_counts[cell_counts$count < 5,]
dim(boxes_lt5) # 68 boxes have less than 5 cells

# exclude boxes with <5 cells
stats1_filtered = stats1[!stats1$Box %in% boxes_lt5$Box, ]
dim(stats1_filtered) # 2326 neurons, 140 neurons were excluded

length(unique(stats1_filtered$Box)) # 142 boxes (original: 210 boxes), 68 boxes were excluded

# create brain-adjusted data (filtered boxes)
stats1_filtered_adjusted = as.data.frame(empiricalBayesLM(stats1_filtered[, numericCols], removedCovariates = stats1_filtered[["Brain"]],
                                                 getEBadjustedData = FALSE)$adjustedData.OLS)
stats1_filtered_adjusted = cbind(stats1_filtered[, c('file_path', 'Brain', 'Type', 'Striatal.Subregion', 'Sex', 'Box')], stats1_filtered_adjusted)
dim(stats1_filtered_adjusted) # 2326, 37

###############################################################################################
############################ Box clustering using box-eigencell ###############################
######### using filtered boxes (boxes with low cell count are excluded) ######################
############################ Explore different param combination ###############################
###############################################################################################
# Select relevant columns
brain_data = stats1_filtered_adjusted[, c('Box', non_redundant_features)]
brain_data$Box = factor(brain_data$Box)

# order brain_data by 'Box'
brain_data = brain_data[order(brain_data$Box), ]

# Transpose and scale the non-redundant features
brain_data_t = t(scale(brain_data[, non_redundant_features])) %>% as.data.frame()

# Set column names based on ordered 'Box'
colnames(brain_data_t) = brain_data$Box

# Calculate eigencell by region
MEs = moduleEigengenes(brain_data_t, brain_data$Box, scale = FALSE)
box_eigencell = MEs$eigengenes

# Process box eigencell for clustering
box_eigencell = t(box_eigencell) %>% as.data.frame()
rownames(box_eigencell) = gsub("^ME", "", rownames(box_eigencell))
box_eigencell$Box = rownames(box_eigencell)

# Prepare data for clustering
box_eigencell_for_cluster = box_eigencell[, !names(box_eigencell) %in% 'Box']
prototypes = list(box_eigencell_for_cluster)

# 1) Calculate dissimilarity matrix using Euclidean distance
cellDist_eucD = lapply(prototypes, dist)

# 2) Calculate dissimilarity matrix using bicorrelation
cellDist_bicor = lapply(prototypes, function(prototype) {
  as.dist(1 - bicor(as.matrix(t(prototype)), use = 'p', maxPOutliers = 0.01))
})

# Set up lists of parameters for hierarchical clustering
hclust_method_list = c('average', 'ward.D2')
deepSplit_list = c(2)
minClusterSize_list = seq(10, 40, by = 5)

# Create combinations of parameters and initialize results
hclust_comb = expand.grid(Method = hclust_method_list, deepSplit = deepSplit_list, minClusterSize = minClusterSize_list)
hclust_comb_results = hclust_comb %>%
  mutate(Silouette_eucD = numeric(n()),
         numCluster_eucD = numeric(n()),
         Silouette_bicor = numeric(n()),
         numCluster_bicor = numeric(n()),
         Param = paste(Method, ", Split:", deepSplit, ", Size:", minClusterSize))

###############################################################################################
# Plot1: dendrogram and heatamap for different parameter combinations
pdf(file = "1.dendro_heatmap_different_params.pdf", width = 8, height = 6)
scpp(2.5)
par(mar = c(5, 5, 5, 5))

  # For each parameter combination
  for (i in seq_len(nrow(hclust_comb))) {
    hclust_method = hclust_comb$Method[i]
    deepSplit_val = hclust_comb$deepSplit[i]
    minClusterSize_val = hclust_comb$minClusterSize[i]
    
    for (method in c('eucD', 'bicor')) {
      cellDist = ifelse(method == 'eucD', cellDist_eucD, cellDist_bicor)
      
      cellTree = lapply(cellDist, hclust, method = hclust_method)
      cellClusters = mapply(cutreeDynamic, cellTree, distM = lapply(cellDist, as.matrix),
                            MoreArgs = list(deepSplit = deepSplit_val, minClusterSize = minClusterSize_val))
      
      sil_score = mean(silhouette(cellClusters, cellDist[[1]])[, 'sil_width'])
      hclust_comb_results[i, paste0('Silouette_', method)] = sil_score
      hclust_comb_results[i, paste0('numCluster_', method)] = length(unique(cellClusters))
      
      # Plot dendrogram and heatmap
      plotDendroAndColors(cellTree[[1]],
                          labels2colors(cellClusters), 
                          'Box clusters',
                          dendroLabels = FALSE,
                          rowText = cellClusters, 
                          rowTextAlignment = 'center',
                          main = paste("Box clustering -", method, hclust_method, deepSplit_val, minClusterSize_val))
      
      cellDist_sorted = as.matrix(cellDist[[1]])[cellTree[[1]]$order, cellTree[[1]]$order]
      print(pheatmap(cellDist_sorted, 
                     cluster_rows = FALSE, 
                     cluster_cols = FALSE,
                     color = colorRampPalette(c("red", "white", "blue"))(100),
                     show_rownames = FALSE, show_colnames = FALSE,
                     legend = FALSE))
    }
  }
dev.off()
###############################################################################################
###############################################################################################
# Plot2: plot for silouette score
pdf(file = "2.silouette_score.pdf", width = 8, height = 6)
scpp(2.5)
par(mar = c(5, 5, 5, 5))
  # Plot silhouette scores
  p = ggplot(hclust_comb_results, aes(x = Param)) +
    geom_line(aes(y = Silouette_eucD, colour = "Silouette_eucD", group = 1), size = 1) + 
    geom_text(aes(y = Silouette_eucD, label = numCluster_eucD), vjust = 0, hjust = 0.5) +
    geom_line(aes(y = Silouette_bicor, colour = "Silouette_bicor", group = 1), size = 1, linetype = "dashed") +
    geom_text(aes(y = Silouette_bicor, label = numCluster_bicor), vjust = 0, hjust = 0.5) +
    ggtitle('Silouette score by param combination\n') +
    ylab('Silouette score') + 
    theme_minimal() + 
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7, face = "bold")
    )
  print(p)
dev.off()
###############################################################################################

###############################################################################################
############################ Box clustering using box-eigencell ###############################
######### using filtered boxes (boxes with low cell count are excluded) ######################
############################ Final clustering ###############################
###############################################################################################
# use bicor, ward.D2, deepsplit = 2, mincluster = 10
# gives 7 clusters (DMs)
hclust_method = 'ward.D2'
deepSplit_val = 2
minClusterSize_val = 10

cluster_tree = list()
cluster_results = list()

# Loop through each method ('eucD' and 'bicor')
for (method in c('eucD', 'bicor')) {
  # Select the appropriate distance matrix
  cellDist = ifelse(method == 'eucD', cellDist_eucD, cellDist_bicor)
  
  # Perform hierarchical clustering
  cellTree = lapply(cellDist, hclust, method = hclust_method)
  
  # Cluster cells and compute silhouette score
  cellClusters = mapply(cutreeDynamic, cellTree, distM = lapply(cellDist, as.matrix),
                        MoreArgs = list(deepSplit = deepSplit_val, minClusterSize = minClusterSize_val))
  
  sil_score = mean(silhouette(cellClusters, cellDist[[1]])[, 'sil_width'])
  print(paste('Silhouette score for', method, 'is:', round(sil_score, 3)))
  
  # Store the results in lists
  cluster_tree[[method]] = cellTree[[1]]
  cluster_results[[method]] = cellClusters
}
# Silouette score for eucD: 0.107
# Silouette score for bicor: 0.227

###############################################################################################
# Plot3: dendrogram and dissimilarity heatmap for the clustering results
pdf(file = "3.dendro_heatmap_clustering.pdf", width = 8, height = 4)
scpp(2.5)
par(mar = c(6.3, 3, 2.5, 1))

  # Function to generate dendrogram and heatmap
  plot_dendro_heatmap = function(method, cellDist, cluster_tree, cluster_results) {
    # Generate color map for clusters
    color_map = setNames(brewer.pal(length(unique(cluster_results[[method]])), 'Set1'), unique(cluster_results[[method]]))
    cluster_colors = color_map[cluster_results[[method]]]
    
    # Plot dendrogram with colored clusters
    plotDendroAndColors(cluster_tree[[method]],
                        cluster_colors,
                        'Box clusters',
                        dendroLabels = FALSE,
                        rowText = cluster_results[[method]], 
                        rowTextAlignment = 'center',
                        main = paste0("Box clustering - All brains combined, ", method))
    
    # Plot sorted heatmap
    sorted = cluster_tree[[method]]$order
    cellDist_sorted = as.matrix(cellDist[[1]])[sorted, sorted]
    pheatmap(cellDist_sorted, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE,
             color = colorRampPalette(c("red", "white", "blue"))(100),
             show_rownames = FALSE, show_colnames = FALSE,
             legend = FALSE)
  }
  
  # Generate dendrograms and heatmaps for both methods
  plot_dendro_heatmap('eucD', cellDist_eucD, cluster_tree, cluster_results)
  plot_dendro_heatmap('bicor', cellDist_bicor, cluster_tree, cluster_results)

dev.off()
###############################################################################################

###############################################################################################
############################ Organize clustering results ###############################
###############################################################################################
# Add cluster results to stats1_adjusted
cluster_bicor_df = data.frame(Box = cluster_tree[['bicor']]$labels, Cluster = cluster_results[['bicor']])
merged_df = merge(stats1_filtered_adjusted, cluster_bicor_df, by = "Box")

# Add x, y, z location of the boxes for 3D visualization
# Brain-adjusted morpho, box, cluster, box xyz info
merged_df_xyz = merge(merged_df, stats1_filtered[, c('file_path', 'x', 'y', 'z', 'CP box')], by = 'file_path')

# Unadjusted morpho, box, cluster, box xyz info
unadjusted_merged_df_xyz = merge(merge(stats1_filtered, cluster_bicor_df, by = "Box"), 
                                 stats1_filtered[, c('file_path', 'x', 'y', 'z', 'CP box')], by = 'file_path')

# Prepare cluster bicor df with xyz
cluster_bicor_df_xyz = merge(cluster_bicor_df, stats1_filtered[, c('Box', 'x', 'y', 'z', 'CP box')], by = 'Box') %>% distinct()
cluster_bicor_df_xyz$Cluster = factor(cluster_bicor_df_xyz$Cluster)
table(cluster_bicor_df_xyz$Cluster)

###############################################################################################
#### Compare box clustering result (filtered boxes vs. all boxes) #############################
###############################################################################################
# read previous result, only keep the neurons in the filtered boxes
prev_merged_df_xyz = read.csv('morpho_brain_adjusted_box_cluster_xyz_D1D2.csv')
prev_merged_df_xyz = prev_merged_df_xyz[prev_merged_df_xyz$Box %in% merged_df$Box, ]
dim(prev_merged_df_xyz) # 2326

# use Rand index for comparison
rand.index(merged_df_xyz$Cluster, prev_merged_df_xyz$Cluster) # 0.91
adjustedRandIndex(merged_df_xyz$Cluster, prev_merged_df_xyz$Cluster) # 0.67


###############################################################################################
############## feature variability assessment (use all boxes) ############################
###############################################################################################
# create brain-adjusted data (all boxes)
stats1_adjusted = as.data.frame(empiricalBayesLM(stats1[, numericCols], removedCovariates = stats1[["Brain"]],
                                                 getEBadjustedData = FALSE)$adjustedData.OLS)
stats1_adjusted = cbind(stats1[, c('file_path', 'Brain', 'Type', 'Striatal.Subregion', 'Sex', 'Box')], stats1_adjusted)
dim(stats1_adjusted) # 2466, 37

# Function to calculate CV for each feature in a box
compute_cv = function(x) {
  if (length(x) > 1) {
    return(sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  } else {
    return(NA) # if only contains 1 neuron, do not calculate
  }
}

cv_df = stats1_adjusted %>%
  group_by(Box) %>%
  summarise(across(non_redundant_features, compute_cv, .names = "CV_{col}"),
            cell_count = n())  %>% data.frame()
head(cv_df)

# convert to long-format to look at CV of all features by each box
cv_long = cv_df %>%
  pivot_longer(cols = starts_with("CV_"), names_to = "Feature", values_to = "CV")

# define a custom theme for ggplot
custom_theme = function(){
  theme_classic() + 
    theme(
      axis.line = element_line(color = "black"),  
      axis.text = element_text(color = "black"),  
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
}

###############################################################################################
# Plot4: Feature CV
pdf(file = "4.Feature_CV.pdf", height = 5, width = 6)
par(mar = c(2, 2, 2, 2))
p = ggplot(cv_long, aes(x = cell_count, y = CV)) +
  geom_point(alpha = 0.5) +  # Scatter plot with transparency
  geom_smooth(method = "loess", color = "blue", se = TRUE) +  # Smoothed trend line
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, by = 0.2)) + 
  labs(
    title = "Feature Variability vs. Cell Count per Box",
    x = "Cell Count per Box",
    y = "Coefficient of Variation (CV)"
  ) +
  custom_theme()
print(p)
dev.off()
###############################################################################################

# categorize boxes by cell count = 5 and count = 10 as thresholds
cv_long = cv_long %>%
  mutate(Neuron_Group = case_when(
    cell_count < 5 ~ "<5 neurons",
    cell_count >= 5 & cell_count < 10 ~ "5-10 neurons",
    cell_count >= 10 ~ ">= 10 neurons"
  ))

cv_long$Neuron_Group = factor(cv_long$Neuron_Group, levels = c('<5 neurons', '5-10 neurons', '>= 10 neurons'))
neuron_group_sizes = table(cv_long$Neuron_Group) / length(non_redundant_features)

###############################################################################################
# Plot5: Feature CV by cell count group
pdf(file = "5.Feature_CV_by_cell_count_group_boxplot.pdf", height = 6, width = 5)
par(mar = c(2, 2, 2, 2))
p = ggplot(cv_long, aes(x = Neuron_Group, y = CV, fill = Neuron_Group)) +
  geom_boxplot() +
  labs(x = "Box Group", y = "Coefficient of Variation (CV)", 
       title = "Box Plot of CV by Neuron Group", fill = 'Box Group') +
  scale_x_discrete(labels = function(x) paste0('Box with ', x, " (n=", neuron_group_sizes[x], ")")) +
  theme_minimal() +
  custom_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
print(p)
dev.off()

# Perform test to compare the group means
print(summary(aov(cv_long$CV ~ Neuron_Group, data = cv_long))) # p = 0.662
###############################################################################################


###############################################################################################
###################### Brain distribution across 7 DMs ########################################
###############################################################################################
# look at whether DM distribution is random across brains (chi-square test)
prev_brain_adjusted = read.csv('morpho_brain_adjusted_box_cluster_xyz_D1D2.csv')
prev_brain_adjusted$Cluster = factor(prev_brain_adjusted$Cluster)

# chi-suqare test
contingency_table = table(prev_brain_adjusted$Brain, prev_brain_adjusted$Cluster)
chi_test = chisq.test(contingency_table)
print(chi_test) # p = 0.25

# stacked bar plot
brain_module_counts = prev_brain_adjusted %>%
  group_by(Cluster, Brain) %>%
  summarise(cnt = n(), .groups = 'drop') %>% ungroup() %>%
  group_by(Cluster) %>% mutate(proportion = cnt / sum(cnt)) %>% ungroup() %>% data.frame()

cluster_sample_size = brain_module_counts %>%
  group_by(Cluster) %>%
  summarise(total_count = sum(cnt), .groups = "drop")

brain_module_counts = brain_module_counts %>%
  left_join(cluster_sample_size, by = "Cluster") %>%
  mutate(Cluster_label = paste0(Cluster, " (n=", total_count, ")"))

brain_module_counts$Cluster = factor(brain_module_counts$Cluster)

# adjust brain name
brain_labels = c(
  "TME07-1" = "P56 WT#1",
  "TME08-1" = "P56 WT#2",
  "TME09-1" = "P56 WT#3",
  "TME10-1" = "P56 WT#4",
  "TME10-3" = "P56 WT#5"
)

###############################################################################################
# Plot6: Bar plot brain proportion in each DM
pdf(file = "6.Brain_proportion_in_DM.pdf", height = 6, width = 6)
par(mar = c(2, 2, 2, 2))
ggplot(brain_module_counts, aes(x = Cluster_label, y = proportion, fill = Brain)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 4, color = "white") +
  theme_minimal() +
  # scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(length(brain_labels), "Set1"),
                    labels = brain_labels) +
  labs(title = "Proportion of Brain Distribution in Each DM",
       x = "DM",
       y = "Proportion",
       fill = "Brain") +
  scale_y_continuous(labels = scales::percent) + 
  custom_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top")
dev.off()
###############################################################################################

###############################################################################################
###################### explore subtle difference between DM2 and DM6 ##########################
###############################################################################################
# 1) effect size - Cohen's D
effect_sizes = sapply(non_redundant_features, 
                      function(f) cohen.d(prev_brain_adjusted[[f]][prev_brain_adjusted$Cluster == '6'], prev_brain_adjusted[[f]][prev_brain_adjusted$Cluster == '2'])$estimate)
effect_size_df = data.frame(Feature = non_redundant_features, Cohen_d = effect_sizes)
effect_size_df = effect_size_df[order(abs(effect_size_df$Cohen_d), decreasing = TRUE),]
head(effect_size_df)

effect_size_df = effect_size_df %>%
  arrange(desc(abs(Cohen_d)))

# get the significance from the t-test comparisons between DM2 and DM6
# Create an empty dataframe to store results
t_test_results = data.frame(Feature = non_redundant_features, p_value = NA)

# Loop through each feature and perform a t-test
for (feature in non_redundant_features) {
  group_2 = prev_brain_adjusted[[feature]][prev_brain_adjusted$Cluster == '2']
  group_6 = prev_brain_adjusted[[feature]][prev_brain_adjusted$Cluster == '6']
  test_result = t.test(group_2, group_6)
  t_test_results$p_value[t_test_results$Feature == feature] = test_result$p.value
}

t_test_results$p_value_adjusted = t_test_results$p_value * length(non_redundant_features)
t_test_results = t_test_results %>%
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

effect_size_df = effect_size_df %>%
  left_join(t_test_results %>% select(Feature, Significance), by = "Feature")

###############################################################################################
# Plot7: Bar plot showing Cohen's D for each non-redundant feature between DM2 and DM6
pdf(file = "7.Barplot_CohenD_DM2vsDM6_Mar.5.2025.pdf", height = 4, width = 5)
par(mar = c(2, 2, 2, 2))
p = ggplot(effect_size_df, aes(x = reorder(Feature, abs(Cohen_d)), y = Cohen_d, fill = Cohen_d)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  ylim(-0.6, 0.6) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label = Significance, 
                hjust = ifelse(Cohen_d > 0, -0.3, 1.3)),  # Adjust position dynamically
            size = 5, color = "black") +  # Set text color to improve visibility
  theme_minimal() + 
  labs(title = "Cohen's d Effect Size by Feature\nDM6 vs. DM2",
       x = "Feature",
       y = "Cohen's d") +
  custom_theme()
print(p)
dev.off()
###############################################################################################


###############################################################################################
# Plot8: box plot of feature distribution between DM2 and DM6 for features with moderate cohen's D
df_subset = prev_brain_adjusted %>% filter(Cluster %in% c("2", "6"))
df_subset$Cluster = factor(df_subset$Cluster, levels = c('2', '6'))

pdf(file = "8.Box_plot_features_DM2vsDM6.pdf", height = 6, width = 4)
par(mar = c(2, 2, 2, 2))
for(feature in effect_size_df[abs(effect_size_df$Cohen_d) > 0.25, 'Feature']){
  p = ggplot(df_subset, aes(x = as.factor(Cluster), y = !!sym(feature))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "lightblue") +  
    geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  
    theme_minimal() +
    labs(title = paste("Box Plot of", feature, "of DM2 and DM6"),
         x = "Dendritic Modules (DMs)",
         y = feature) +
    custom_theme()
  print(p)
}
dev.off()
###############################################################################################

###############################################################################################
# 2) LDA seperation
# LDA reduces the feature space to a lower-dimensional space while maximizing class separation
lda_model = lda(Cluster ~ ., data = df_subset[, c("Cluster", non_redundant_features)])
lda_values = predict(lda_model)$x
df_subset$LDA1 = lda_values[,1]

# Plot9: LDA projection visualization
pdf(file = "9.LDA_projection_DM2vsDM6.pdf", height = 5, width = 6)
par(mar = c(2, 2, 2, 2))
ggplot(df_subset, aes(x = LDA1, fill = Cluster)) +
  geom_density(alpha = 0.7) +
  theme_minimal() +
  custom_theme() + 
  labs(title = "LDA Projection for DM2 vs. DM6", fill = 'DM')
dev.off()

lda_predictions = predict(lda_model)$class
# LDA prediction accuracy
print(mean(lda_predictions == df_subset$Cluster)) # 0.682
###############################################################################################
