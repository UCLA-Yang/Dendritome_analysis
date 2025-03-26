# work directory
setwd("C:/Users/yanyanming77/Desktop/Ming_Transition/Morphometric_Analysis/D1D2_box_analysis-Mar.25.2025")

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
library(ggpubr)

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
# Extract non-redundant features from the first cluster
cluster_df = clusters[[1]] %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "Feature") %>%
  rename(Group = 2) 

non_redundant_features = cluster_df %>%
  arrange(Group) %>%
  group_by(Group) %>%
  slice(1) %>%
  pull(Feature)

print(length(non_redundant_features)) # 21 non-redundant features

# create brain-adjusted data
stats1_adjusted = as.data.frame(empiricalBayesLM(stats1[, numericCols], removedCovariates = stats1[["Brain"]],
                                                 getEBadjustedData = FALSE)$adjustedData.OLS)
stats1_adjusted = cbind(stats1[, c('file_path', 'Brain', 'Type', 'Striatal.Subregion', 'Sex', 'Box')], stats1_adjusted)
dim(stats1_adjusted) # 2466, 37

###############################################################################################
################ Create brain median normalized feature ################################
###############################################################################################
# normalize features to its median by brain
normalized_df = stats1 %>% group_by(Brain) %>%
  mutate(across(all_of(numericCols), ~ . / median(., na.rm = TRUE))) %>% as.data.frame()

# calculate mean in region by brain using normalized df
complete_mean_df = normalized_df %>% group_by(Brain, Box) %>%
  summarise(across(all_of(numericCols), mean, na.rm = TRUE)) %>%
  ungroup()

# create a combined data for uniform color bar range
min_color = min(complete_mean_df[numericCols], na.rm = TRUE)
max_color = max(complete_mean_df[numericCols], na.rm = TRUE)

# calculate mean in box using normalized df, all brains combined
complete_mean_df_all = normalized_df %>% group_by(Box) %>%
  summarise(across(all_of(numericCols), mean, na.rm = TRUE)) %>%
  ungroup()

###############################################################################################
############################ non-empty Box distribution #######################################
###############################################################################################
n_box_all = length(unique(stats1$Box))
print(n_box_all) # 210 non-empty boxes
n_box_D1 = length(unique(stats1[stats1$Type == 'D1', 'Box'])) 
print(n_box_D1) # 196 non-empty boxes for D1 MSNs
n_box_D2 = length(unique(stats1[stats1$Type == 'D2', 'Box'])) 
print(n_box_D2)# 185 non-empty boxes for D2 MSNs

###############################################################################################
############################ Plot: Cell count distribution in boxes #################
###############################################################################################
# 1. Plot: cell count per box for each brain
# Define a function to generate normalized cell count in each box by brain
cnt_plot_by_brain = function(data, subset) {
  n_box = switch(subset,
                 'All' = n_box_all,
                 'D1' = n_box_D1,
                 'D2' = n_box_D2,
                 stop("Invalid subset"))
  
  if(subset != 'All') {
    data = subset(data, Type == subset)
  }
  
  # 1). Distribution of cell counts by each box
  cell_counts = data %>%
    group_by(Brain, Box) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(freq = count / sum(count))
  
  # Plot cell counts by box
  p = ggplot(cell_counts, aes(x = Box, y = count, fill = Brain)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~Brain, scales = "free_y", ncol = 1) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank(), 
      axis.text = element_text(color = "black", size = 10), 
      axis.title = element_text(color = "black", hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(x = "Box Number", y = "Frequency", title = paste0("Frequency of Box Counts by Brain, ", subset))
  
  print(p)
  
  # 2). Distribution of box count by cell count in the box
  count_frequencies = cell_counts %>%
    group_by(Brain, count) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(prop = n / n_box)
  
  # Plot box count frequencies
  p1 = ggplot(count_frequencies, aes(x = factor(count), y = prop, fill = factor(count))) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", 100 * prop)), vjust = -0.5, color = "black", size = 3) +
    facet_wrap(~ Brain, scales = "free_y", ncol = 1) +
    theme_minimal() +
    ylim(0, max(count_frequencies$prop) + 0.2) + 
    labs(x = "Cell count per box", y = "Frequency", 
         title = paste0("Frequency of boxes for each cell count, by brain - ", subset)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.position = "none",  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "black", hjust = 0.5),
          plot.title = element_text(hjust = 0.5)) 
  
  print(p1)
  
  return(cell_counts)
}

pdf(file = "1_Cell_count_box_by_brain.pdf", width = 7, height = 7)
scpp(2)
  cnt_data_list = list(
    'All' = cnt_plot_by_brain(stats1, 'All'),
    'D1' = cnt_plot_by_brain(stats1, 'D1'),
    'D2' = cnt_plot_by_brain(stats1, 'D2')
  )
dev.off()

# 2. Plot: box cell count all brains combined
cnt_plot_all_brains_combined = function(data) {
  # 1). Distribution of cell counts by each box
  cell_counts = data %>%
    group_by(Box) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(freq = count / sum(count))
  
  # Plot cell counts by box
  p = ggplot(cell_counts, aes(x = Box, y = count)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank(), 
      axis.text = element_text(color = "black"), 
      axis.title = element_text(color = "black", hjust = 0.5),
      axis.text.x = element_text(size = 2),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(x = "Box Number", y = "Frequency", title = "Frequency of Cell Counts All Brains Combined")
  
  print(p)
  
  # 2). Distribution of box count by cell count in the box
  count_frequencies = cell_counts %>%
    group_by(count) %>%
    summarise(n = n(), .groups = 'drop') %>% 
    mutate(prop = n / n_box_all)
  
  # Plot box count frequencies
  p1 = ggplot(count_frequencies, aes(x = factor(count), y = prop, fill = factor(count))) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", 100 * prop)), vjust = -0.5, color = "black", size = 3) +
    theme_minimal() +
    ylim(0, max(count_frequencies$prop) + 0.2) + 
    labs(x = "Cell count per box", y = "Frequency", 
         title = "Frequency of boxes for each cell count all brains combined") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          legend.position = "none",  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "black", hjust = 0.5),
          plot.title = element_text(hjust = 0.5)) 
  
  print(p1)
  
  # 3). CDF plot
  data_cdf = count_frequencies %>%
    arrange(count) %>%
    mutate(cumulative_n = cumsum(n),
           total_n = sum(n),
           cdf = cumulative_n / total_n)
  
  p2 = ggplot(data_cdf, aes(x = count, y = cdf)) +
    geom_line() +
    labs(title = "CDF of Box Cell Count, All Brains Combined",
         x = "Cell per Box",
         y = "Cumulative Distribution") +
    scale_x_continuous(breaks = seq(0, max(data_cdf$count), by = 2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.05)) +
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), 
          axis.title = element_text(color = "black", hjust = 0.5),
          plot.title = element_text(hjust = 0.5)) 
  
  print(p2)
}

pdf(file = "2_Cell_count_box_all_brain.pdf", width = 15, height = 7)
scpp(2)
  cnt_plot_all_brains_combined(stats1)
dev.off()

###############################################################################################
############################ Box clustering using box-eigencell ###############################
############################ Explore different param combination ###############################
###############################################################################################
# Select relevant columns
brain_data = stats1_adjusted[, c('Box', non_redundant_features)]
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

# 3. Plot: dendrogram and heatamap for different parameter combinations
pdf(file = "3_dendro_heatmap_different_params.pdf", width = 8, height = 6)
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

# 4. Plot: plot for silouette score
pdf(file = "4_silouette_score.pdf", width = 8, height = 6)
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
############################ Box clustering using box-eigencell ###############################
############################ Final clustering ###############################
###############################################################################################
# use bicor, ward.D2, deepsplit = 2, mincluster = 15 
hclust_method = 'ward.D2'
deepSplit_val = 2
minClusterSize_val = 15

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
# Silouette score for eucD: 0.095
# Silouette score for bicor: 0.167

# 5. Plot: dendrogram and dissimilarity heatmap for the clustering results
pdf(file = "5_dendro_heatmap_clustering.pdf", width = 8, height = 4)
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
############################ Organize clustering results ###############################
###############################################################################################
# Add cluster results to stats1_adjusted
cluster_bicor_df = data.frame(Box = cluster_tree[['bicor']]$labels, Cluster = cluster_results[['bicor']])
merged_df = merge(stats1_adjusted, cluster_bicor_df, by = "Box")

# Add x, y, z location of the boxes for 3D visualization
# Brain-adjusted morpho, box, cluster, box xyz info
merged_df_xyz = merge(merged_df, stats1[, c('file_path', 'x', 'y', 'z', 'CP box')], by = 'file_path')

# Unadjusted morpho, box, cluster, box xyz info
unadjusted_merged_df_xyz = merge(merge(stats1, cluster_bicor_df, by = "Box"), 
                                 stats1[, c('file_path', 'x', 'y', 'z', 'CP box')], by = 'file_path')

# Prepare cluster bicor df with xyz
cluster_bicor_df_xyz = merge(cluster_bicor_df, stats1[, c('Box', 'x', 'y', 'z', 'CP box')], by = 'Box') %>% distinct()
cluster_bicor_df_xyz$Cluster = factor(cluster_bicor_df_xyz$Cluster)
table(cluster_bicor_df_xyz$Cluster)

# Save CSV files
write.csv(cluster_bicor_df, 'box_cluster_bicor_D1D2.csv')
write.csv(merged_df_xyz, 'morpho_brain_adjusted_box_cluster_xyz_D1D2.csv', row.names = FALSE)
write.csv(unadjusted_merged_df_xyz, 'morpho_un_adjusted_box_cluster_xyz_D1D2.csv', row.names = FALSE)
write.csv(cluster_bicor_df_xyz, 'box_cluster_bicor_D1D2_with_xyz.csv', row.names = FALSE)

# Define cluster order for dendrogram and prepare barplot data
dendro_cluster_order = c('7','3', '5', '6','2', '4', '1')
cluster_bicor_df_xyz_barplot = cluster_bicor_df_xyz

# 6. Plot: bar plot distribution of boxes in each cluster
pdf(file = "6_box_count_per_DM_barplot.pdf", height = 5, width = 8)
par(mar = c(2, 2, 2, 2))
  p = ggplot(cluster_bicor_df_xyz_barplot, aes(x = Cluster, fill = Cluster)) +
    geom_bar() +
    geom_text(stat = 'count', aes(label = after_stat(count), y = after_stat(count)), vjust = -0.5, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, max(table(cluster_bicor_df_xyz_barplot$Cluster)) + 5), expand = c(0, 0)) +
    scale_fill_manual(values = setNames(brewer.pal(length(unique(cluster_bicor_df_xyz_barplot$Cluster)), "Set1"),
                                        paste0(1:length(unique(cluster_bicor_df_xyz$Cluster))))) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.line = element_line(color = "black"),
          panel.grid.major = element_blank(),
          axis.title = element_text(size = 14, face = 'bold'),
          axis.text.x = element_text(size = 14, face = 'bold'),
          axis.text.y = element_text(size = 14, face = 'bold')) +
    labs(title = "Box Count Distribution by Dendritic Modules (DMs)", 
         x = "Dendritic Modules (DMs)", y = "Box Count") + 
    guides(fill=guide_legend(title="DM"))
  print(p)
dev.off()

###############################################################################################
############################ Connected components in each cluster ###############################
###############################################################################################
# calculate distance between boxes and generate a graph-object
distance_matrix = as.matrix(dist(cluster_bicor_df_xyz[, c("x", "y", "z")]))
graph = graph_from_adjacency_matrix(distance_matrix <= 500, mode = "undirected", diag = FALSE)

V(graph)$name = cluster_bicor_df_xyz$Box
V(graph)$cluster = cluster_bicor_df_xyz$Cluster

# get the connected box group in each cluster, and store connected box count in each group
results = list()
clusters1 = seq(1, length(unique(cluster_bicor_df_xyz$Cluster)))
for (cluster in clusters1) {
  subgraph = induced_subgraph(graph, V(graph)$cluster == cluster)
  comps = components(subgraph)
  component_data = lapply(1:max(comps$membership), function(i) {
    list(
      size = sum(comps$membership == i),
      boxes = V(subgraph)$name[comps$membership == i]
    )
  })
  results[[cluster]] = component_data
}

# store the results to a dataframe
data_for_plot = data.frame(Cluster = character(), ComponentSize = integer(), ComponentID = character())

# Loop through each cluster to populate the data frame
for (cluster in clusters1) {
  cluster_data = results[[cluster]]
  data_for_plot = rbind(data_for_plot, 
                        data.frame(Cluster = paste0('Cluster ', cluster), 
                                   ComponentSize = sapply(cluster_data, function(x) x$size)))
}

level_cluster = paste("Cluster", seq_len(length(unique(cluster_bicor_df_xyz$Cluster))))
data_for_plot$Cluster = factor(data_for_plot$Cluster, levels = level_cluster)

data_for_plot = data_for_plot %>%
  arrange(Cluster, desc(ComponentSize)) %>%
  group_by(Cluster) %>%
  mutate(ComponentID = as.character(row_number())) %>%
  ungroup() %>%
  as.data.frame()

levels_componentID = as.character(seq_len(length(unique(data_for_plot$ComponentID))))
data_for_plot$ComponentID = factor(data_for_plot$ComponentID, levels = levels_componentID)

# generate bar plot distribution
# statistical test of distribution
# compare to uniform distribution
convert_pvalue = function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Perform chi-square goodness-of-fit test
chisq_test_result = sapply(clusters1, function(cluster) {
  observed = subset(data_for_plot, Cluster == paste0('Cluster ', cluster))$ComponentSize
  expected_prob = rep(1 / length(observed), length(observed))
  chisq_result = chisq.test(observed, p = expected_prob, simulate.p.value = TRUE, B = 2000)
  p = chisq_result$p.value
  paste0(sprintf("%.3e", p), '', convert_pvalue(p))
})
print(chisq_test_result)

chisq_test_statistic = sapply(clusters1, function(cluster) {
  observed = subset(data_for_plot, Cluster == paste0('Cluster ', cluster))$ComponentSize
  expected_prob = rep(1 / length(observed), length(observed))
  chisq_result = chisq.test(observed, p = expected_prob, simulate.p.value = TRUE, B = 2000)
  chisq_result$statistic
})

# write the test result to an excel file
chisq_df = data.frame('statistic' = chisq_test_statistic, 'p value' = chisq_test_result)
chisq_df$p.value = gsub('\\*', '', chisq_df$p.value)
chisq_df$DM = paste0('DM-', as.character(1:7))
chisq_df = chisq_df[, c(3, 1, 2)]

write.xlsx(chisq_df, "chisq_test_DM_component_distribution.xlsx")

# 7. Plot: bubble plot of distribution of connected components within each cluster
asterisks = c('Cluster 1' = '***', 'Cluster 2' = '***', 'Cluster 3' = '***', 
              'Cluster 4' = '   ', 'Cluster 5' = '** ', 
              'Cluster 6' = '***', 'Cluster 7' = '*  ')

pdf(file = paste0('7_distribution_of_connected_boxes_bubble.pdf'), height = 5, width = 7)
par(mar = c(2,2,2,2));
  p = ggplot(data_for_plot, aes(x = ComponentID, y = Cluster)) +
    geom_point(aes(size = ComponentSize), shape = 21, fill = "#3399ff") +
    scale_size_area(max_size = 10) +  # Adjust max_size to control bubble size
    theme_minimal() +
    scale_y_discrete(labels = function(labels) {
      sapply(labels, function(l) {
        paste0(l, ' <span style="color:red;">', asterisks[l], '</span>')
      })
    }) +
    labs(x = "Component ID", y = "Cluster", size = "Component Size") +
    theme(
      legend.position = "top", 
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      plot.title = element_text(face = "bold", hjust = 0.5),   
      axis.text.x = element_text(angle = 45, size = 10, face = 'bold'),
      axis.text.y = element_text(size = 10, face = 'bold', hjust = 0),
      axis.title.x = element_text(face = "bold"),  
      axis.title.y = element_text(face = "bold"), 
      strip.text = element_text(face = "bold"),    
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),          
      panel.grid.minor = element_blank()
    ) +
    theme(axis.text.y = element_markdown())
  print(p)
dev.off()

###############################################################################################
############################ Filter to only keep >= 7 boxes as MTs ###############################
###############################################################################################
filtered_boxes = unlist(lapply(clusters1, function(cluster) {
  valid_indices = sapply(results[[cluster]], function(group) group$size >= 7)
  unlist(lapply(results[[cluster]][valid_indices], function(group) group$boxes))
}))

print(length(filtered_boxes))  # 101 unique boxes

# Filter the dataframe based on the filtered boxes
df_filtered = subset(cluster_bicor_df_xyz, Box %in% filtered_boxes)

# write to a csv file
write.csv(df_filtered, 'box_cluster_bicor_D1D2_with_xyz_filtered_thresh7.csv', row.names = FALSE)

# Plot 7.1: bar plot distribution of box count per MT (filtered DM)
pdf(file = "7.1_box_count_per_MT_barplot.pdf", height = 5, width = 7)
  par(mar = c(2, 2, 2, 2))
  
  pal = setNames(brewer.pal(length(unique(cluster_bicor_df_xyz_barplot$Cluster)), "Set1"),
                 paste0(1:length(unique(cluster_bicor_df_xyz$Cluster))))
  pal_pop = pal[-4]
  
  p = ggplot(df_filtered, aes(x = Cluster, fill = Cluster)) +
    geom_bar() +
    geom_text(stat = 'count', aes(label = after_stat(count), y = after_stat(count)), vjust = -0.5, color = "black") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, max(table(cluster_bicor_df_xyz_barplot$Cluster)) + 5), expand = c(0, 0)) +
    scale_fill_manual(values = pal_pop) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.line = element_line(color = "black"),
          panel.grid.major = element_blank(),
          axis.title = element_text(size = 14, face = 'bold'),
          axis.text.x = element_text(size = 14, face = 'bold'),
          axis.text.y = element_text(size = 14, face = 'bold')) +
    labs(title = "Box Count Distribution by Morphological Territories (MTs)", 
         x = "Morphological Territories (MTs)", y = "Box Count") + 
    guides(fill=guide_legend(title="MT"))
  print(p)
dev.off()

# Plot 7.2: heatmap showing box count by MT and 11 striatal communities
cell_community = unadjusted_merged_df_xyz[unadjusted_merged_df_xyz$Box %in% df_filtered$Box, ]
dim(cell_community) # 1478 cells

# rename ACB2 to ACB
cell_community = cell_community %>% mutate(Striatal.Community = ifelse(Striatal.Community == 'ACB2', 'ACB', Striatal.Community))
table(cell_community$Striatal.Community, useNA = "ifany")

# define order of community
order_community = c("CPr_m", 'CPr_imd', "CPr_imv", "CPr_l",
                    "CPi_dm", "CPi_dl", "CPi_vm", "CPi_vl",
                    "CPc_d", "CPc_i", "CPc_v",
                    'ACB', 'Other non-CP')

cell_community$Striatal.Community = factor(cell_community$Striatal.Community, levels = order_community)
cell_community$Cluster = factor(cell_community$Cluster, levels = c(1,2,3,5,6,7))

# create data for heatmap use
tab = cell_community %>% 
  group_by(Cluster, Striatal.Community) %>% 
  summarise(unique_box_count = n_distinct(Box)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Striatal.Community, values_from = unique_box_count, values_fill = 0) %>% 
  as.data.frame() %>% mutate(Cluster = paste0('MT-', Cluster))
rownames(tab) = tab$Cluster
tab = tab[, -1]
tab = tab[, order_community]

pdf(file = "7.2_heatmap_box_count_in_MTs_by_community.pdf", wi = 8, he = 5)
par(mar = c(3,3,3,3));
  ht = Heatmap(as.matrix(tab), 
               col = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),  
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               name = "Box Count",
               column_names_rot = 45,
               row_names_side = "left",              
               column_names_side = "bottom",            
               width = unit(6, "inch"),          
               height = unit(3, "inch"),         
               heatmap_legend_param = list(
                 legend_direction = "horizontal"
               )
  )
  draw(ht, heatmap_legend_side = 'top')
dev.off()

###############################################################################################
############################ Find representative cell for each MT ###############################
###############################################################################################
cell_rank_in_box_cluster = function(){
  brain_data = merged_df_xyz[, c('file_path', 'Cluster', 'Box', non_redundant_features)]
  # filter boxes to keep only connected boxes
  connected_box_df = df_filtered
  connected_boxes = as.vector(connected_box_df[, 'Box'])
  
  brain_data = brain_data[brain_data$Box %in% connected_boxes, ]
  
  brain_data = brain_data[order(brain_data$Cluster), ]
  brain_data_t = as.data.frame(t(scale(brain_data[, non_redundant_features])))
  colnames(brain_data_t) = brain_data$file_path
  
  # compute loadings of each cell to the cluster first PC
  cell_loadings_df = data.frame(file_path = character(), 
                                Cluster = character(), 
                                Loading = numeric(), Rank_in_box_cluster = numeric())
 
  for(c in names(table(brain_data$Cluster))){
    
    # the cell names in each box
    moduleData = brain_data_t[, brain_data$Cluster == c] %>% as.data.frame()
    colnames(moduleData) = unlist(brain_data[brain_data$Cluster == c, 'file_path'])
    
    pcaResult = prcomp(moduleData, scale. = FALSE, center = TRUE)
    loadings = abs(pcaResult$rotation[, 1])
    named_loadings = setNames(loadings, colnames(moduleData))
    sorted_loadings = sort(named_loadings, decreasing = TRUE)
    
    ranked_loadings = data.frame(
      file_path = names(sorted_loadings),
      Cluster = c,
      Loading = sorted_loadings,
      Rank_in_box_cluster = rank(-sorted_loadings) 
    )
    
    cell_loadings_df = rbind(cell_loadings_df, ranked_loadings)
  }
  
  temp = t(brain_data_t)
  
  return(list(cell_loadings_df = cell_loadings_df, scaled_data = temp))
}

cell_rank_in_box_cluster_df = cell_rank_in_box_cluster()$cell_loadings_df
scaled_data = cell_rank_in_box_cluster()$scaled_data

# output the top 25 cells in each connected  box cluster
cell_rank_in_box_cluster_df_sub = subset(cell_rank_in_box_cluster_df, Rank_in_box_cluster<=25)
cell_rank_in_box_cluster_df_sub$Brain = sapply(strsplit(as.character(cell_rank_in_box_cluster_df_sub$file_path), "_"), function(x) x[2])
table(cell_rank_in_box_cluster_df_sub$Cluster, cell_rank_in_box_cluster_df_sub$Brain)

# write to a csv file
write.csv(cell_rank_in_box_cluster_df_sub, 'cell_rank_in_box_cluster_top25.csv', row.names = FALSE)

###############################################################################################
############################ Plot feature covariation of representative cells ###############################
###############################################################################################
# Filter scaled_data for rows that are in cell_rank_in_box_cluster_df_sub and convert to data frame
scaled_data_sub = scaled_data[rownames(scaled_data) %in% cell_rank_in_box_cluster_df_sub$file_path, , drop = FALSE] %>%
  as.data.frame() %>%
  mutate(file_path = rownames(.))

# Merge data by 'file_path'
joined = merge(cell_rank_in_box_cluster_df_sub, scaled_data_sub, by = "file_path", all.x = TRUE)

# Reshape the data from wide to long format
joined_long = joined %>%
  mutate(Cluster = as.integer(Cluster)) %>%
  pivot_longer(cols = all_of(non_redundant_features), names_to = "Variable", values_to = "Value")

# Define the desired order and grouping of the 21 non-redundant features
desired_nonredundant_feature_order = c("Bif_ampl_local", "Bif_ampl_remote", "Bif_tilt_local", "Bif_tilt_remote", "Bif_torque_local", "Bif_torque_remote",
                                       "Length", "Sum_EucDistance", "Max_EucDistance", "ABEL_All", "ABEL_Internal", "Height", "Width", "Depth",
                                       "N_stems", "Terminal_degree", "Branch_Order", "Fractal_Dim", "Partition_asymmetry", "Balancing_Factor", "Convexity")

# Reorder the 'Variable' factor
joined_long$Variable = factor(joined_long$Variable, levels = desired_nonredundant_feature_order)


colors_nonredundant_feature = c(rep('darkorange', 6), rep('lightseagreen', 8), rep('mediumorchid', 7))
names(colors_nonredundant_feature) = paste0(1:length(non_redundant_features))

color_map_bicor = setNames(brewer.pal(length(unique(cluster_results[['bicor']])), 'Set1'), unique(cluster_results[['bicor']]))
color_map_bicor[6] = 'gold1' # change yellow to gold for clearer view

# 8. Plot: the covariation of scaled features of 25 representative neurons
pdf(file = "8_scaled_features_top25_cells_per_cluster_connected_boxes.pdf", width = 7, height = 3)
scpp(2.5)
  for (c in c(1, 2, 3, 5, 6, 7)) {
    plot_df = subset(joined_long, Cluster == c & !Rank_in_box_cluster %in% c(1)) %>% as.data.frame()
    
    p = ggplot(plot_df, aes(x = Variable, y = Value, group = file_path)) +
      ylim(range(plot_df$Value)) +
      xlab(paste0('Cluster ', c)) + 
      geom_line(data = subset(plot_df, file_path != 'eigencell'), size = 0.4, color = color_map_bicor[c]) + 
      geom_point(data = subset(plot_df, file_path != 'eigencell'), size = 0.6, color = color_map_bicor[c]) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(color = colors_nonredundant_feature, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      scale_x_discrete(labels = paste0('', seq_along(non_redundant_features)))
    print(p)
}
dev.off()

# 9. Plot: the covariation of scaled features of 10 representative neurons
pdf(file = "9_scaled_features_top10_cells_per_cluster.pdf", width = 7, height = 3)
scpp(2.5)
for (c in c(1, 2, 3, 5, 6, 7)) {
  plot_df = subset(joined_long, Cluster == c & !Rank_in_box_cluster %in% c(1) & Rank_in_box_cluster <= 10)
  
  p = ggplot(plot_df, aes(x = Variable, y = Value, group = file_path)) +
    xlab(paste0('Cluster ', c)) + 
    geom_line(size = 0.2, color = color_map_bicor[c]) +
    geom_point(size = 0.4, color = color_map_bicor[c]) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  print(p)
}
dev.off()

# 10. Plot: heatmap visualization of correlation of top25 neurons and cluster eigencell
# Reshape and reorder the cell rank data
cell_rank_in_box_cluster_df_sub_reshaped = dcast(cell_rank_in_box_cluster_df_sub, Rank_in_box_cluster ~ Cluster, value.var = "Loading")
rownames(cell_rank_in_box_cluster_df_sub_reshaped) = paste0('Rank ', cell_rank_in_box_cluster_df_sub_reshaped$Rank_in_box_cluster)
cell_rank_in_box_cluster_df_sub_reshaped = cell_rank_in_box_cluster_df_sub_reshaped[, -1]

# Reorder columns according to dendrogram cluster order
dendro_cluster_order_filtered = c('7', '3', '5', '6', '2', '1')
cell_rank_in_box_cluster_df_sub_reshaped = cell_rank_in_box_cluster_df_sub_reshaped[, dendro_cluster_order_filtered]
colnames(cell_rank_in_box_cluster_df_sub_reshaped) = paste0('Box Cluster ', dendro_cluster_order_filtered)

# Create matrix and reverse row order rank (low to high)
mat = as.matrix(cell_rank_in_box_cluster_df_sub_reshaped)
mat_reversed = mat[nrow(mat):1, ]

# Generate heatmap
pdf(file = "10_corr_heatmap_top25neurons_connected_cluster_eigencell.pdf", width = 9, height = 12)
par(cex.lab = 0.6, cex.axis = 0.6, mar = c(2, 2, 2, 2))
  heatmap(mat_reversed,
          Rowv = NA, Colv = NA, scale = "column",
          main = "Heatmap correlation\ntop 25 neurons and connected cluster eigencell", 
          xlab = "", ylab = "Top 25 neurons")
dev.off()

# Reset graphical parameters
par(cex.lab = 1, cex.axis = 1, las = 1, mar = c(5, 4, 4, 2) + 0.1)

###############################################################################################
######### Heatmap: Normalized feature distribution by box, ordered by cluster #################
###############################################################################################
summary_type = stats1 %>% mutate(Box = as.character(Box)) %>% count(Box, Type) %>%
  spread(key = Type, value = n, fill = 0) 
rownames(summary_type) = summary_type$Box

summary_data = list(Type = summary_type)

# define 3 groups of features and assign colors
features_angle = c('Bif_ampl_local', 'Bif_ampl_remote', 'Bif_tilt_local', 'Bif_tilt_remote', 'Bif_torque_local', 'Bif_torque_remote', 'Centripetal_Bias')
features_length = c('Depth', 'Height', 'Width', 'Length', 'Sum_EucDistance', 'Sum_PathDistance', 'Max_EucDistance', 'Max_PathDistance', 'ABEL_All', 'ABEL_Internal', 'ABEL_Terminal', 'BAPL_All', 'BAPL_Internal', 'BAPL_Terminal')
features_complexity = c('N_bifs', 'N_branch', 'N_tips', 'N_stems', 'Branch_Order', 'Fractal_Dim', 'Partition_asymmetry', 'Terminal_degree', 'Balancing_Factor', 'Convexity')

feature_reorder = data.frame(Feature = c(features_angle, features_length, features_complexity),
                             group = rep(1:3, times = c(length(features_angle), length(features_length), length(features_complexity))))


draw_heatmap_distribution_by_cluster = function(data_for_heatmap, data_for_barplot, group_name, bar_name, dendro_cluster_order) {
  
  if(group_name != 'All'){
    data_for_heatmap = subset(data_for_heatmap, Brain == group_name & Box %in% data_for_barplot$Box)
    box_names = data_for_heatmap$Box
    data_for_heatmap = as.data.frame(t(subset(data_for_heatmap, select = -c(Brain, Box))))
    colnames(data_for_heatmap) = box_names
  }else{
    data_for_heatmap = subset(data_for_heatmap, Box %in% data_for_barplot$Box)
    box_names = data_for_heatmap$Box
    data_for_heatmap = as.data.frame(t(subset(data_for_heatmap, select = -c(Box))))
    colnames(data_for_heatmap) = box_names
  }
  
  # re-order feature by 3 main  groups
  # sort data matrix by re-ordered features
  data_for_heatmap = data_for_heatmap[match(feature_reorder$Feature, rownames(data_for_heatmap)), ]
  data_for_heatmap = round(data_for_heatmap,2)
  
  # add annotation
  # 1) row annotation (feature grouping)
  # add annotation
  annotation_df = data.frame(Cluster = feature_reorder[,'group'])
  rownames(annotation_df) = feature_reorder$Feature
  
  palette_colors = c("darkorange", "lightseagreen", "mediumorchid")
  annotation_color_mapping = setNames(palette_colors, c('1', '2', '3'))
  
  right_annotation = HeatmapAnnotation(df = annotation_df, 
                                     which = "row", 
                                     annotation_name_side = NULL,
                                     show_legend = FALSE,
                                     col = list(Cluster = annotation_color_mapping))
  
  # create color map for cells
  blue_white_red = colorRamp2(c(min_color, 1, max_color),
                            c("blue", "white", "red"))

  # create data matrix for text display
  text_matrix = as.matrix(data_for_heatmap)
  text_matrix[is.na(text_matrix)] = ''
  
  # 2. column annotation (boxes grouped and ordered by clusters)
  cluster_bicor_df_sub = cluster_bicor_df %>% filter(Box %in% box_names)
  cluster_bicor_df_sub$Cluster = factor(cluster_bicor_df_sub$Cluster, levels = dendro_cluster_order)
  ordered_boxes = as.vector(cluster_bicor_df_sub[order(cluster_bicor_df_sub$Cluster), ]$Box)
  
  data_for_heatmap = data_for_heatmap[, ordered_boxes]
  data_for_heatmap = round(data_for_heatmap, digits = 8)
  
  group_annotation = data.frame(
  Cluster = factor(paste0('Cluster', as.vector(cluster_bicor_df_sub[order(cluster_bicor_df_sub$Cluster), ]$Cluster)))
  )
  rownames(group_annotation) = colnames(data_for_heatmap) 
  
  color_mapping = setNames(brewer.pal(length(levels(group_annotation$Cluster)), 'Set1'), levels(group_annotation$Cluster))
  
  ha = HeatmapAnnotation(df = group_annotation, show_legend = TRUE, 
                       col = list(Cluster = color_mapping), annotation_name_side = 'left',
                       annotation_legend_param = list(
                         Cluster = list(
                           title = "Dendritic Moduel (DM)",    
                           labels = c("DM1", "DM2", "DM3", "DM4", "DM5", "DM6", "DM7")  
                         )
                       ))
  
  # bar plot annotation to show D1/D2 distribution
  rownames(data_for_barplot) = data_for_barplot$Box
  
  stacked_bar_annotation = columnAnnotation(`Type` = anno_barplot(data_for_barplot[ordered_boxes, c("D1", "D2")],
                                                                  which = "column", 
                                                                  bar_width = 1, 
                                                                  border = TRUE, 
                                                                  gp = gpar(fill = c('darkgreen', 'darkred'))))
  lgd = Legend(title = "Type", at =  c("D1", "D2"), legend_gp = gpar(fill = c('darkgreen', 'darkred')))
  
  # use prettified brain name
  if(group_name == "TME07-1"){
    lab = "P56_WT#1"
  }else if(group_name == "TME08-1"){
    lab = "P56_WT#2"
  }else if(group_name == "TME09-1"){
    lab = "P56_WT#3"
  }else if(group_name == "TME10-1"){
    lab = "P56_WT#4"
  }else if(group_name == "TME10-3"){
    lab = "P56_WT#5"
  }else if(group_name == 'All'){
    lab = 'All'
  }
  
  # generate heatmap
  ht = Heatmap(as.matrix(data_for_heatmap), 
             name = " ", 
             col = blue_white_red,
             cluster_rows = FALSE, 
             cluster_columns = FALSE, 
             show_row_names = TRUE, row_title = 'Morphometrics',
             show_column_names = TRUE, column_title = 'Box',
             width = unit(12, 'inch'), 
             height = unit(3, 'inch'),
             row_names_gp = gpar(fontsize= 8),
             column_names_gp = gpar(fontsize = 3), column_names_rot = 90,
             show_heatmap_legend = FALSE,
             top_annotation = stacked_bar_annotation,
             right_annotation = right_annotation,
             bottom_annotation = ha
  ) 
  draw(ht, heatmap_legend_side = 'top', annotation_legend_side = "left", annotation_legend_list = lgd)
  grid.text(paste0('Distribution of Normalized Features across Boxes: ', lab),
          x = unit(0.5, 'npc'),
          y = unit(0.95, 'npc'),
          gp = gpar(fontsize = 12, fontface = 'bold'))
}

# 11. Plot: heatmap visualization of normalized feature by box ordered by box clusters
pdf(file = "11_morpho_distribution_by_cluster_brain_ordered_pretty_name.pdf", wi = 18, he = 6)
scpp(2.5);
par(mar = c(3,5,3,5));
  for(brain in unique(complete_mean_df$Brain)){
    draw_heatmap_distribution_by_cluster(complete_mean_df, summary_data$Type, brain, 'Type', dendro_cluster_order)
  }
  draw_heatmap_distribution_by_cluster(complete_mean_df_all, summary_data$Type, 'All', 'Type', dendro_cluster_order)
dev.off()

# save image in high res for shiny dashboard
png(file = paste0('D1D2_box_cluster_shiny/www/', "morpho_dist_by_cluster.png"), width = 400, height = 150, units='mm', res = 600)
scpp(2.5);
par(mar = c(1,1,1,1));
  draw_heatmap_distribution_by_cluster(complete_mean_df_all, summary_data$Type, 'All', 'Type', dendro_cluster_order)
dev.off()

###############################################################################################
############# Stacked bar for D1 D2 distribution within cluster ######################
###############################################################################################
# barplot distribution of D1/D2, and with appropriate annotations of sample sizes and statistical test results
stacked_bar_distribution = function(data, by_var, by_var_name, by_var_levels, stack_var, set_name) {
  # Calculate overall proportions
  overall_proportions = data %>%
    count(!!sym(stack_var)) %>%
    mutate(
      !!sym(by_var) := "All", 
      Proportion = n / sum(n)
    )
  
  # Calculate cluster proportions
  cluster_totals = data %>%
    group_by(!!sym(by_var)) %>%
    summarise(total_per_cluster = n(), .groups = 'drop')
  
  cluster_proportions = data %>%
    group_by(!!sym(by_var), !!sym(stack_var)) %>%
    summarise(n = n(), .groups = 'drop') %>%
    left_join(cluster_totals, by = by_var) %>%
    mutate(Proportion = n / total_per_cluster)
  
  cluster_proportions$Cluster = factor(cluster_proportions$Cluster)
  
  combined_data = bind_rows(overall_proportions, cluster_proportions) %>%
    arrange(!!sym(by_var))
  
  # Reorder 'All' group to be first
  combined_data$Cluster = factor(combined_data$Cluster, levels = c('All', by_var_levels))

  # Group counts and labels
  grp_cnts = combined_data %>% group_by(!!sym(by_var)) %>% summarise(n = sum(n))
  grp_cnt_strings = paste(grp_cnts[[by_var]], " (", grp_cnts$n, ")", sep = "")
  
  # Chi-square test
  pval_chi = c("") # Placeholder for 'All' group
  expected_prop = prop.table(table(data[[stack_var]]))
  
  # Chi-square goodness of fit test for each cluster
  for (grp in by_var_levels) {
    ct = chisq.test(table(data[data[[by_var]] == grp, stack_var]), p = expected_prop)
    print(paste0('group ', grp, ' - ', ct$p.value))
    pval_chi = append(pval_chi, ifelse(ct$p.value <= 0.05, '*', ''))
  }
  
  # Add asterisks for p-values
  asterisks = rep(pval_chi, each = length(unique(data[[stack_var]])))
  combined_data = combined_data[order(combined_data$Cluster), ]
  combined_data$Asterisks = asterisks
  combined_data$Proportion = round(combined_data$Proportion, 3)
  
  # Create the stacked bar plot
  p = ggplot(combined_data, aes(x = get(by_var), y = Proportion, fill = get(stack_var), label = Asterisks)) +
    geom_bar(stat = "identity", width = 0.5) +
    labs(
      title = paste0("Proportion of ", stack_var, " by ", by_var_name, "\n(", set_name, ")"),
      x = by_var_name,
      y = "Proportion",
      fill = stack_var
    ) +
    scale_x_discrete(labels = grp_cnt_strings) +
    geom_text(aes(label = scales::percent(Proportion)), position = position_stack(vjust = 0.5), 
              size = 3.5, color = "white") +
    geom_text(aes(y = 1, label = Asterisks), vjust = 0, color = "red") +    
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() + 
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 14),  
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),  
      axis.text.y = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
    scale_fill_manual(values = c(`D1` = "darkgreen", `D2` = "darkred"))
  
  print(p)
}

# 12. Plot: bar plot distribution of Type in DM
pdf(paste0("12_Proportion_Type_DM.pdf"), wi = 7, he = 5)
  stacked_bar_distribution(merged_df, 'Cluster', 'Dendritic Modules (DMs)', 
                           as.character(seq(1:length(unique(merged_df[, 'Cluster'])))), 'Type', 'All')
dev.off()

# 12.1 Plot: bar plot distribution of Type in MT
filtered_merged_df = merged_df %>% filter(Box %in% filtered_boxes)
pdf(paste0("12.1_Proportion_Type_MT.pdf"), wi = 7, he = 5)
  stacked_bar_distribution(filtered_merged_df, 'Cluster', 'Morpholocial Territories (MTs)', 
                           c('1', '2', '3', '5', '6', '7'), 'Type', 'All')
dev.off()

###############################################################################################
############# Box plot comparison of morphometrics between clusters ######################
###############################################################################################
# 13. Plot:  feature distribution box plot across clusters by brain and all brains combined
pdf(file = "13_morpho_distribution_boxplot_by_cluster_brain.pdf", width = 6, height = 8)
scpp(2.5)
par(mar = c(3, 5, 3, 5))
  for (feature in numericCols) {
    # Boxplot by brain
    p = ggplot(merged_df_xyz, aes(x = factor(Cluster), y = !!sym(feature), fill = factor(Cluster))) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Set1") +
      labs(
        x = "Box Cluster", 
        y = feature, 
        title = paste0("Distribution of ", feature, " across box clusters"), 
        fill = 'Cluster'
      ) +
      facet_wrap(~ Brain, scales = "free", ncol = 1) + 
      theme_minimal() + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    print(p)
    
    # Boxplot for all brains combined
    kw_test = kruskal.test(as.formula(paste(feature, ' ~ Cluster')), data = merged_df_xyz)
    p_value_adjusted = p.adjust(kw_test$p.value, method = 'bonferroni', n = length(numericCols))
    
    annot_string = ifelse(p_value_adjusted <= 0.05, 
                          paste0('Kruskal-Wallis test: P = ', signif(p_value_adjusted, 3)), 
                          'Kruskal-Wallis test: NS.')
    
    p1 = ggplot(merged_df_xyz, aes(x = factor(Cluster), y = !!sym(feature), fill = factor(Cluster))) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Set1") +
      labs(x = "Box Cluster", 
        y = feature, 
        title = paste0("Distribution of ", feature, " across box clusters"), 
        subtitle = annot_string,
        fill = 'Cluster'
      ) +
      theme_minimal() + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5) 
      ) 
    print(p1)
  }
dev.off()

###############################################################################################
############# heatmap: D2 vs. D1 in each box, boxes grouped by cluster ######################
###############################################################################################
# 14. Plot: D2 vs. D1 in each box, boxes grouped by cluster
# only for boxes with certain amount of cells
initialize_empty_placeholders = function(){
  # list of boxes that have insufficient cell type
  insufficient_box_num <<- c()
  # list of boxes that have insufficient cell type during modeling
  insufficient_box_model_num <<- c()
  # list of boxes that have NA for p-vals that indicate potential issues during modeling
  p_val_NA_num <<- c()
  # list of boxes that have sufficent cells and have results
  sufficient_box_num <<- c()
}

color_gradient = function(coef, neg_log10p) {
  if (coef > 0) {
    return(colorRampPalette(c("white", "red"))(100)[min(ceiling(neg_log10p*10),100)])
  } else {
    return(colorRampPalette(c("white", "blue"))(100)[min(ceiling(neg_log10p*10),100)])
  }
}

colored_heatmap_by_type_in_box_kruskal_by_cluster = function(data, data_for_barplot, set_name, bar_name, dendro_cluster_order){
  data_matrix = matrix(NA, nrow = nNumStats)
  color_matrix = matrix(NA, nrow =nNumStats, ncol = 1)
  text_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  coef_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  
  for(box in sort(unique(data$Box))){
    print(box)
    
    results = data.frame(Feature = numericCols,
                         Coefficient = NA, 
                         P_Value = NA)
    
    for (feature in results$Feature) {
      
      sub_data = data[data$Box == box,]
      
      # exclude boxes that have fewer than min_cell_count_per_box cell count
      if(nrow(sub_data) < min_cell_count_per_box | length(unique(sub_data$Type)) < 2){
        if(!box %in% insufficient_box_num){
          insufficient_box_num = c(insufficient_box_num, box)
        }
      }else{
        
        sub_data[, 'Type'] = factor(sub_data[, 'Type'], levels = order_type)
        column_means = colMeans(sub_data[sub_data$Type == 'D1', numericCols], na.rm = TRUE)
        results$Baseline_mean = column_means
        
        # if encounters any errors when fitting the model
        tryCatch({
          
          kw_test = kruskal.test(sub_data[, feature] ~ sub_data[, 'Type'])
          p_value = kw_test$p.value
          
          # use difference in means D2 - D1
          coef = mean(sub_data[sub_data$Type == 'D2', feature]) - mean(sub_data[sub_data$Type == 'D1', feature])
          
          # some may have insufficient data that leads to p-value to be NaN
          if(is.na(p_value)){
            # if p_val is NA
            if(!box %in% p_val_NA_num){
              p_val_NA_num = c(p_val_NA_num, box)
            }
          }else{
            results[results$Feature == feature, "Coefficient"] = coef
            results[results$Feature == feature, "P_Value"] = p_value
          }
        }, error = function(e){
          if(!box %in% insufficient_box_model_num){
            insufficient_box_model_num <= c(insufficient_box_model_num, box)
          }
          print(paste0(box, ' has an error occured in the test'))
        })
      }
    }
    
    if(!box %in% c(insufficient_box_num, insufficient_box_model_num, p_val_NA_num)){
      # done when all features for a box is finished
      results$Fold_change = round(abs(results$Coefficient) / results$Baseline_mean,1)
      results$P_Value = signif(results$P_Value, 2)
      results[results$P_Value == 1, 'P_Value'] = 0.95 # account for P-value = 1
      
      results$neg_log10p = -log10(results$P_Value)
      results$Coefficient = signif(results$Coefficient,2)
      results$Color = mapply(color_gradient, results$Coefficient, results$neg_log10p)
      
      # if non-significant, change the color to light 
      results[results$P_Value >= 0.95 & results$Coefficient > 0, 'Color'] = '#FFE7E7' 
        results[results$P_Value >= 0.95 & results$Coefficient < 0, 'Color'] = '#E7E7FF'
          # if fold change < 0.05, set to white
        results[results$Fold_change < 0.05, 'Color'] = 'white'
          
        data_matrix = cbind(data_matrix, results$neg_log10p)
        color_matrix = cbind(color_matrix, results$Color)
        text_matrix = cbind(text_matrix, results$Coefficient)
        coef_matrix = cbind(coef_matrix, results$Coefficient)
        
        # add to column names
        if(!box %in% sufficient_box_num){
          sufficient_box_num = c(sufficient_box_num, box)
        }
    }
  }
  
  # add column names (box number) to the matrices
  matrices_list = list(data_matrix, color_matrix, text_matrix, coef_matrix)
  
  drop_first_column_and_rename = function(mat, col_name){
    mat = mat[, -1]
    colnames(mat) = col_name
    return(mat)
  }
  new_matrices_list = lapply(matrices_list, function(x) drop_first_column_and_rename(x, sufficient_box_num))  
  data_matrix = new_matrices_list[[1]]
  color_matrix = new_matrices_list[[2]]
  text_matrix = new_matrices_list[[3]]
  coef_matrix = new_matrices_list[[4]]
  
  # add * based on -log10P, if p<0.05 then -log10P should > 1.301
  text_matrix = matrix(as.character(text_matrix), nrow=nrow(text_matrix))
  colnames(text_matrix) = sufficient_box_num
  text_matrix[data_matrix <= 1.301] = ''
  text_matrix[data_matrix > 1.301] = '*'
  text_matrix[data_matrix > 2.301] = '**'
  text_matrix[data_matrix > 3.301] = '***'
  
  # set row and column names of the matrix
  colLabs = colnames(data_matrix)
  rowLabs = results$Feature
  rownames(data_matrix) = rownames(color_matrix)= rownames(text_matrix) = rownames(coef_matrix) = rowLabs

  # to map cell colors using color_matrix
  custom_color_fun = function(value, matrix = coef_matrix, colors = color_matrix) {
    color_index = match(value, matrix)
    return(colors[color_index])
  }
  
  data_matrix = data_matrix[feature_reorder$Feature, ]
  color_matrix = color_matrix[feature_reorder$Feature, ]
  text_matrix = text_matrix[feature_reorder$Feature, ]
  coef_matrix = coef_matrix[feature_reorder$Feature, ]

  annotation_df = data.frame(Cluster = feature_reorder[,'group'])
  rownames(annotation_df) = feature_reorder$Feature
  
  palette_colors = c("darkorange", "lightseagreen", "mediumorchid")
  annotation_color_mapping = setNames(palette_colors, c('1', '2', '3'))
  
  right_annotation = HeatmapAnnotation(df = annotation_df, 
                                       which = "row", 
                                       annotation_name_side = NULL,
                                       show_legend = FALSE,
                                       col = list(Cluster = annotation_color_mapping))
  
  # 2. column annotation
  # filter Boxes based on all_boxes dataset
  cluster_bicor_df_sub = cluster_bicor_df %>% filter(Box %in% colnames(data_matrix))
  cluster_bicor_df_sub$Cluster = factor(cluster_bicor_df_sub$Cluster, levels = dendro_cluster_order)
  ordered_boxes = as.vector(cluster_bicor_df_sub[order(cluster_bicor_df_sub$Cluster), ]$Box)
  
  # reorder all matrices
  data_matrix = data_matrix[, ordered_boxes]
  color_matrix = color_matrix[, ordered_boxes]
  text_matrix = text_matrix[, ordered_boxes]
  coef_matrix = coef_matrix[, ordered_boxes]
  
  group_annotation = data.frame(
    Cluster = factor(paste0('Cluster', as.vector(cluster_bicor_df_sub[order(cluster_bicor_df_sub$Cluster), ]$Cluster)))
  )
  rownames(group_annotation) = colnames(data_matrix) 
  
  color_mapping = setNames(brewer.pal(length(levels(group_annotation$Cluster)), 'Set1'), paste0('Cluster', 1:length(levels(group_annotation$Cluster))))
  
  ha = HeatmapAnnotation(df = group_annotation, show_legend = TRUE, 
                         col = list(Cluster = color_mapping), annotation_name_side = 'left')
  
  # bar plot annotation
  rownames(data_for_barplot) = data_for_barplot$Box
  stacked_bar_annotation = columnAnnotation(`Type` = anno_barplot(data_for_barplot[ordered_boxes, c("D1", "D2")],
                                                                  which = "column", 
                                                                  bar_width = 1, 
                                                                  border = TRUE, 
                                                                  gp = gpar(fill = c('darkgreen', 'darkred'))))
  lgd = Legend(title = "Type", at =  c("D1", "D2"), legend_gp = gpar(fill = c('darkgreen', 'darkred')))
  
  ht = Heatmap(coef_matrix, 
               col = custom_color_fun,  
               name = " ",  
               show_row_names = TRUE, 
               show_column_names = TRUE,  
               cluster_rows = FALSE,  
               cluster_columns = FALSE,  
               width = unit(12, 'inch'), 
               height = unit(3, 'inch'),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 6))
               },
               top_annotation = stacked_bar_annotation,
               right_annotation = right_annotation,
               bottom_annotation = ha,
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 2), column_names_rot = 90,
               show_heatmap_legend = FALSE
  )
  
  draw(ht, heatmap_legend_side = 'top', annotation_legend_side = "left", annotation_legend_list = lgd)
  grid.text(paste0("Coefficient and Significance Heatmap (D2 vs. D1), Non-Parametric Test, ", set_name),
            x = unit(0.5, "npc"),
            y = unit(0.98, "npc"),  
            gp = gpar(fontsize = 12, fontface = 'bold'))
  
  # Define the viewport to place the legend (adjust x and y as needed)
  viewport = viewport(x = 0.5, y = 0.94, width = 1, height = 0.5)
  pushViewport(viewport)
  
  # Add the custom legend
  grid.rect(gp = gpar(fill = "red"), x = 0.2, width = 0.06, height = 0.03)
  grid.text("D2 higher than D1", x = 0.24, just = "left", gp = gpar(fontsize = 8))
  grid.rect(gp = gpar(fill = "blue"), x = 0.5, width = 0.06, height = 0.03)
  grid.text("D2 lower than D1", x = 0.54, just = "left", gp = gpar(fontsize = 8))
  
  # Reset viewport
  popViewport()
}

pdf(file = "14_Heatmap-Association-Type-Blue-Red-D2vsD1-in-box.pdf", wi = 18, he = 7);
scpp(2.5);
par(mar = c(3,5,3,5));
  min_cell_count_per_box = 10
  initialize_empty_placeholders()
  colored_heatmap_by_type_in_box_kruskal_by_cluster(stats1, summary_data$Type,'Min cell count: 10', 'Type', dendro_cluster_order)
  
  min_cell_count_per_box = 20
  initialize_empty_placeholders()
  colored_heatmap_by_type_in_box_kruskal_by_cluster(stats1, summary_data$Type, 'Min cell count: 20', 'Type', dendro_cluster_order)
dev.off()

###############################################################################################
############# heatmap of mean distribution of each box cluster ######################
###############################################################################################
# 15. Plot: heatmap of mean distribution of each box cluster
pdf(file = "15_morphometric_mean_distribution_by_cluster.pdf", wi = 5, he = 8)
scpp(2.5);
par(mar = c(3,5,3,5));
  a = merged_df_xyz[, c('Cluster', numericCols)] %>%
    group_by(Cluster) %>%
    summarise(across(all_of(numericCols), mean, na.rm = TRUE), .groups = 'drop') %>%
    as.data.frame()
  a = t(a)
  colnames(a) = a[1,]
  a = a[-1,]
  a = a[, dendro_cluster_order]
  
  # add row annotation
  annotation_df = data.frame(Cluster = feature_reorder[,'group'])
  rownames(annotation_df) = feature_reorder$Feature
  
  palette_colors = c("darkorange", "lightseagreen", "mediumorchid")
  annotation_color_mapping = setNames(palette_colors, c('1', '2', '3'))
  
  right_annotation = HeatmapAnnotation(df = annotation_df, 
                                       which = "row", 
                                       annotation_name_side = NULL,
                                       show_legend = FALSE,
                                       col = list(Cluster = annotation_color_mapping))
  
  a = a[rownames(annotation_df), ]
  z_scores = apply(a, 1, scale)
  z_scores_matrix = matrix(z_scores, nrow = nrow(a), byrow = TRUE)
  rownames(z_scores_matrix) = rownames(a)
  colnames(z_scores_matrix) = colnames(a)
  
  ht = Heatmap(z_scores_matrix, 
          name = " ",  
          show_row_names = TRUE, 
          show_column_names = TRUE,  
          cluster_rows = FALSE,  
          cluster_columns = FALSE,  
          width = unit(2, 'inch'), 
          height = unit(6, 'inch'),
          right_annotation = right_annotation,
          row_names_gp = gpar(fontsize= 8),
          column_names_gp = gpar(fontsize = 8), column_names_rot = 90,
          show_heatmap_legend = TRUE,
  )
  draw(ht)
  grid.text(paste0('Morphometric mean distribution by box cluster'),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
dev.off()

# 16. Plot: heatmap of mean distribution of each box cluster, breaking down by brain
pdf(file = "16_morphometric_mean_distribution_by_cluster_by_brain.pdf", wi = 8, he = 8)
scpp(2.5);
par(mar = c(3,5,3,5));
  a = merged_df_xyz[, c('Cluster', 'Brain', numericCols)] %>%
    group_by(Cluster, Brain) %>%
    summarise(across(all_of(numericCols), mean, na.rm = TRUE), .groups = 'drop') %>%
    as.data.frame()
  a = a %>% arrange(factor(Cluster, levels = dendro_cluster_order), Brain)
  a$`Cluster-Brain` = paste(a$Cluster, a$Brain, sep = "-")
    a = t(a)
  colnames(a) = a['Cluster-Brain',]
  a = a[numericCols,]
  a = a %>% as.data.frame() %>% mutate(across(everything(), as.numeric))
  a = a[rownames(annotation_df), ]
  z_scores = apply(a, 1, scale)
  z_scores_matrix = matrix(z_scores, nrow = nrow(a), byrow = TRUE)
  rownames(z_scores_matrix) = rownames(a)
  colnames(z_scores_matrix) = colnames(a)
  
  ht = Heatmap(z_scores_matrix, 
          name = " ",  
          show_row_names = TRUE, 
          show_column_names = TRUE,  
          cluster_rows = FALSE,  
          cluster_columns = FALSE,  
          width = unit(4, 'inch'), 
          height = unit(6, 'inch'),
          right_annotation = right_annotation,
          row_names_gp = gpar(fontsize= 8),
          column_names_gp = gpar(fontsize = 8), column_names_rot = 90,
          show_heatmap_legend = TRUE
  )
  draw(ht)
  grid.text(paste0('Morphometric mean distribution by box cluster by brain'),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
dev.off()

#########################################################################
######################### Box Membership ################################
#########################################################################

################### Box module membership calculation ###########
# box_eigencell: the eigencell of each box, need to add the DM of each box to the box_eigencell
box_eigencell_DM = merge(box_eigencell, cluster_bicor_df_xyz, by = 'Box')

### Use Correlations ####
# calculate DM eigencell 
# order brain-adjusted morphometrics by DM
data_for_membership = merged_df_xyz
data_for_membership = data_for_membership[order(data_for_membership$Cluster), ]
data_for_membership_t = t(scale(data_for_membership[, non_redundant_features])) %>% as.data.frame()
colnames(data_for_membership_t) = data_for_membership$Cluster

# Calculate eigencell by DM
DM_MEs = moduleEigengenes(data_for_membership_t, data_for_membership$Cluster, scale = FALSE)
DM_eigencell = DM_MEs$eigengenes
DM_eigencell = t(DM_eigencell) %>% as.data.frame()
rownames(DM_eigencell) = gsub("^ME", "", rownames(DM_eigencell))
DM_eigencell$Cluster = rownames(DM_eigencell)

################### 1) Box membership scatter plot - All DMs ######################
# calucate correlation of Box eigencell vs. DM eigencell #####
# Initialize an empty data frame to store results
box_member_result = data.frame(Box = character(), correlation = numeric(), stringsAsFactors = FALSE)
# list of unique Boxes
boxes = unique(box_eigencell_DM$Box)

# Iterate through each Box
for (box in boxes) {
  
  # Extract the row for the current Box
  box_row = box_eigencell_DM[box_eigencell_DM$Box == box, ]
  # Extract the corresponding DM eigencell (ensure one row is selected)
  DM_row = DM_eigencell[DM_eigencell$Cluster == box_row$Cluster, non_redundant_features, drop = FALSE]
  # convert to vector
  dm_vector = as.numeric(DM_row)
  
  # Extract the Box eigencell as a numeric vector
  box_vector = as.numeric(box_row[, non_redundant_features])
  
  # Compute correlation
  corr_value = cor(box_vector, dm_vector)
  
  # Store the result in the data frame
  box_member_result = rbind(box_member_result, data.frame(Box = box, correlation = corr_value))
}

# calculate cell count per box
cell_counts_per_box = merged_df_xyz %>%
  group_by(Box) %>%
  summarise(cell_count = n(), .groups = 'drop') %>% as.data.frame()
# add the cell count of each box to box_member_result
box_member_result = merge(box_member_result, cell_counts_per_box, by = 'Box')
# add Cluster to box_member_result
box_member_result = merge(box_member_result, merged_df_xyz[, c('Box', 'Cluster')], by = 'Box')
# convert Cluster to factor
box_member_result$Cluster = factor(box_member_result$Cluster, levels = c("1", "2", "3", "4", "5", "6", "7"))

# Visualization settings
# Define the custom theme function
custom_theme = function() {
  theme(
    axis.line = element_line(color = "black"), 
    panel.grid = element_blank(),  
    axis.title = element_text(size = 10, face = "bold"),  
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),  
    panel.background = element_blank()
  )
}

# Create a named vector mapping cluster levels to colors
color_mapping = setNames(brewer.pal(length(as.character(1:7)), "Set1")  , as.character(1:7))
color_mapping['6'] = 'gold'

### 17. Plot: Box membership Scatter plot - all DMs at once
pdf(file = "17_box_membership_scatterplot_AllDMs.pdf", wi = 6, he = 5)
scpp(2.5);
par(mar = c(3,5,3,5));
p = ggplot(data = box_member_result, aes(x = cell_count, y = correlation, color = Cluster)) +
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values = color_mapping,
                     labels = paste0('DM', as.character(1:7))) + 
  # scale_color_brewer(palette = 'Set1')+
  labs(title = "Scatter Plot of Box Membership in DM",
       subtitle = paste("Correlation: ", round(cor(box_member_result[, 'correlation'], box_member_result$cell_count), 2)),
       x = "Cell count per box",
       y = ("Correlation in DM"),
       color = 'DM') +
  theme_minimal() + custom_theme()
print(p)
dev.off()

################### 2) Box membership scatter plot - by each DM ######################
### 18. Plot: Box membership Scatter plot - each DMs at a time
# Function to create scatter plots for each DM
plot_list = lapply(c("1", "2", "3", "4", "5", "6", "7"), function(cluster_level) {
  # Subset data for the specific DM
  subset_df = box_member_result[box_member_result$Cluster == cluster_level, ]
  
  # Calculate correlation
  cor_value = cor(subset_df[, 'correlation'], subset_df$cell_count, use = "complete.obs")
  
  # Create scatter plot
  ggplot(subset_df, aes(x = cell_count, y = correlation)) +
    geom_point(size = 3, alpha = 0.7, color = color_mapping[[cluster_level]]) +
    labs(title = paste("DM:", cluster_level),
         subtitle = paste("Correlation: ", round(cor_value, 2)),
         x = "Cell count per box",
         y = ("Correlation in DM")) +
    theme_minimal() + custom_theme()
})

# Arrange the 7 plots into a single giant plot
final_plot = ggarrange(plotlist = plot_list, ncol = 4, nrow = 2)  
final_plot = annotate_figure(final_plot, 
                             top = text_grob("Box membership based on correlation vs. cell count per box", 
                                             face = "bold", 
                                             size = 16, 
                                             hjust = 0.5))
pdf(file = "18_box_membership_scatterplot_each_DM.pdf", wi = 10, he = 6)
scpp(2.5);
par(mar = c(3,5,3,5));
  print(final_plot)
dev.off()

################### 3) Box membership in each DM, heatmap ######################
# calucate correlation of Box eigencell vs. DM eigencell for each DM
box_member_result_all = data.frame(Box = character(), 
                                   correlation = numeric(), 
                                   Cluster = character(),
                                   stringsAsFactors = FALSE)

# Iterate through each Box
for (box in boxes) {
  # Extract the row for the current Box
  box_row = box_eigencell_DM[box_eigencell_DM$Box == box, ]
  # Iterate through each DM
  for (cluster in c("1", "2", "3", "4", "5", "6", "7")){
    # Extract the corresponding DM eigencell (ensure one row is selected)
    DM_row = DM_eigencell[DM_eigencell$Cluster == cluster, non_redundant_features, drop = FALSE]
    # convert to vector
    dm_vector = as.numeric(DM_row)
    
    # Extract the Box eigencell as a numeric vector
    box_vector = as.numeric(box_row[, non_redundant_features])
    
    # Compute correlation
    corr_value = cor(box_vector, dm_vector)
    
    # Store the result in the data frame
    box_member_result_all = rbind(box_member_result_all, 
                                  data.frame(Box = box, Cluster = cluster, correlation = corr_value)) 
  }
}

# add the cell count of each box to box_member_result_all
box_member_result_all = merge(box_member_result_all, cell_counts_per_box, by = 'Box')
# convert Cluster to factor
box_member_result_all$Cluster = factor(box_member_result_all$Cluster, levels = c("1", "2", "3", "4", "5", "6", "7"))

# sort the boxes by Cluster and cell_count in ascending order, and place in columns, place Cluster in rows
# Convert long format to wide format: Clusters as rows, Boxes as columns
heatmap_matrix = box_member_result_all %>%
  select(Box, Cluster, correlation) %>%
  pivot_wider(names_from = Box, values_from = correlation) %>% 
  arrange(Cluster) %>% as.matrix()

# Convert to matrix and set row names
heatmap_matrix = as.matrix(heatmap_matrix[, -1])  # Remove Cluster column from data
heatmap_matrix = apply(heatmap_matrix, 2, as.numeric)
rownames(heatmap_matrix) = paste0("DM", c("1", "2", "3", "4", "5", "6", "7"))

# Re-order Boxes (columns) by Cluster (DM) and cell_count
box_eigencell_DM = merge(box_eigencell_DM, cell_counts_per_box, by = 'Box')
box_eigencell_DM$Cluster = factor(box_eigencell_DM$Cluster, levels = dendro_cluster_order)
box_order = box_eigencell_DM %>%
  distinct(Box, Cluster, cell_count) %>%  
  arrange(factor(box_eigencell_DM$Cluster), cell_count) %>% 
  pull(Box) 

# Reorder columns of heatmap_matrix based on sorted Box order
heatmap_matrix = heatmap_matrix[, box_order]
heatmap_matrix = heatmap_matrix[paste0('DM',dendro_cluster_order),]
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Annotation (bottom) - 1
group_annotation = box_eigencell_DM[, c('Box','Cluster')] %>% 
  mutate(Box = factor(Box, levels = box_order)) %>%
  arrange(Box) %>% select(Cluster)

ha_bottom = HeatmapAnnotation(df = group_annotation, show_legend = TRUE, 
                              annotation_name_gp = gpar(fontsize = 0),
                              col = list(Cluster = color_mapping), annotation_name_side = 'left',
                              annotation_legend_param = list(
                                Cluster = list(
                                  title = "Dendritic Moduel (DM)",    
                                  labels = paste0('DM', dendro_cluster_order)))) 

# Annotation (bottom) - 2
box_eigencell_DM = box_eigencell_DM %>% 
  mutate(`Cell count` = ifelse(cell_count < 5, "less than 5", "greater than 5"))
box_eigencell_DM$`Cell count` = factor(box_eigencell_DM$`Cell count`, levels = c('less than 5', 'greater than 5'))
group_annotation1 = box_eigencell_DM[, c('Box','Cell count')] %>% 
  mutate(Box = factor(Box, levels = box_order)) %>%
  arrange(Box) %>% select(`Cell count`)
color_mapping1 = c("less than 5" = "#33ffff", "greater than 5" = "#ff99cc")
ha_bottom1 = HeatmapAnnotation(df = group_annotation1, show_legend = TRUE, 
                               col = list(`Cell count` = color_mapping1), 
                               annotation_name_gp = gpar(fontsize = 0),
                               annotation_legend_param = list(
                                 `Cell count` = list(
                                   title = "Box Size",    
                                   labels = c("less than 5 cells", "greater than 5 cells"))))
ha_bottom_combined = c(ha_bottom, ha_bottom1)

# Annotation (right)
ha_right = HeatmapAnnotation(df = data.frame(Cluster = dendro_cluster_order), 
                             which = "row", 
                             annotation_name_gp = gpar(fontsize = 0),
                             show_legend = FALSE, 
                             col = list(Cluster = color_mapping))

# Create and draw the heatmap
ht = Heatmap(
  heatmap_matrix, 
  name = "Correlation",
  col = col_fun,
  cluster_rows = FALSE,  
  cluster_columns = FALSE, 
  show_column_names = FALSE,  
  show_row_names = FALSE,
  column_title = "",
  row_title = "DMs",
  width = unit(ncol(heatmap_matrix) * 1, "mm"),
  height = unit(nrow(heatmap_matrix) * 15, "mm"),
  bottom_annotation = ha_bottom_combined,
  right_annotation = ha_right,
  
  # Move color legend to the top and make it horizontal
  heatmap_legend_param = list(
    title = "Correlation",  
    legend_direction = "horizontal")
)

### Plot 19: Heatmap - Box membership for each DM
pdf(file = "19_Heatmap_box_membership_by_each_DM_each_Box.pdf", wi = 12, he = 7)
scpp(2.5);
par(mar = c(3,5,3,5));
  draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right")
  grid.text("Correlation between box eigencell and DM eigencell\n(Rows: 7 DMs, Columns: Ordered Boxes)", 
            y = unit(0.92, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()


