# Set work directory
setwd("...")

# source functions
func_dir = '../Functions/'

source(paste0(func_dir, "individualAnalysis-General-010.R"));
source(paste0(func_dir, "GNVFunctions-018-02.R"));
source(paste0(func_dir, "networkFunctions-extras-20.R"));
source(paste0(func_dir, "labelPoints2-01.R"));
source(paste0(func_dir, "heatmap.wg.R"));

# import libraries
library(anRichmentMethods)
library(Cairo)
library(gtable)
library(ggplotify)
library(tsne)
library(umap)
library(vioplot)
library(gplots)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(grid)
library(gridExtra)
library(car)
library(RColorBrewer) 
library(ggtext)
library(randomForest)
library(randomGLM)
library(openxlsx);
library(dplyr)
library(stringr)
library(effsize)
library(nlme)
library(lme4)
library(emmeans)
library(data.table)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(pheatmap)
library(tibble)
library(wesanderson)
library(fmsb)

# create directions
dir.create("Plots", recursive = TRUE);
dir.create("Results", recursive = TRUE);

###############################################################################################
############################ read data ########################################################
###############################################################################################
# # read 2-m data
# data_2m = read.csv('morpho_2248CPneurons_with_community.csv')
# # read 12-m data
# data_12m = read.csv('all_htme_brains_with_registration_1168CPneurons_onlyCP.csv')
# ncol(data_12m) # 37
# # only keep 12-m WT data
# data_12m = subset(data_12m, Genotype == 'WT')
# 
# # drop uncommon columns for both datasets
# col_not_common = setdiff(union(colnames(data_2m), colnames(data_12m)), intersect(colnames(data_2m), colnames(data_12m)))
# data_2m = data_2m[, !(names(data_2m) %in% col_not_common)]
# data_2m$Age =  '2m'
# ncol(data_2m) # 36
# data_12m = data_12m[, !(names(data_12m) %in% col_not_common)]
# data_12m$Age = '12m'
# ncol(data_12m) # 36
# 
# # combine data
# stats1 = rbind(data_2m, data_12m)
# write.csv(stats1, 'morpho_2791CPneurons_for_Aging_analysis.csv',row.names = FALSE)

# read all morphometrics of WT neurons (P56 and 12m, CP neurons only)
stats1 = read.csv('morpho_2791CPneurons_for_Aging_analysis.csv')
table(stats1$Brain)
table(stats1$Age) # 2m: 2248, 12m: 543
table(stats1$Striatal.Subregion) # 473 CPr, 1377 CPi, 941 CPc
table(stats1$Type) # 1459 D1, 1332 D2

# create Brain_ind variable, for example ('TME07-1 (2m,WT))
stats1$Brain_ind = paste0(stats1$Brain, ' (', stats1$Age, ',WT)')
table(stats1$Brain_ind)

# re-order the columns
colnames(stats1)
stats1 = stats1[, c(1:4, 36, 37, 5:35)]
colnames(stats1)

###############################################################################################
############################ Rename, Add columns, Process data ################################
###############################################################################################
stats = list(stats1); 
statNames = c("2m_vs_12m"); 
statNames.pretty = c("");
names(stats) = statNames;

# create covariate variables
stats = lapply(stats, function(stats) 
{
  # Age in subregion (CPr, CPi, CPc)
  rest1 = restrictVariableByCovariateLevels(stats$Age, stats$Striatal.Subregion, 
                                            varName = "Age", covarName = "", nameSep = " in ", check.names = FALSE);
  # Age in type (D1, D2)
  rest2 = restrictVariableByCovariateLevels(stats$Age, stats$Type, 
                                            varName = "Age", covarName = "", nameSep = " in ", check.names = FALSE);
  
  stats = data.frame.ncn(rest1, rest2, stats);
  out = setNames(stats, gsub(" ", "_", names(stats)))
});

# add Age in type in subregions
stats = lapply(stats, function(stats) {
  # Age in D1 in subregions
  rest1 = restrictVariableByCovariateLevels(stats$`Age_in_D1`, stats$Striatal.Subregion, 
                                            varName = "Age_in_D1", covarName = "", nameSep = " in ", check.names = FALSE);
  # Age in D2 in subregions
  rest2 = restrictVariableByCovariateLevels(stats$`Age_in_D2`, stats$Striatal.Subregion, 
                                            varName = "Age_in_D2", covarName = "", nameSep = " in ", check.names = FALSE);
  stats = data.frame.ncn(rest1, rest2, stats);
  out = setNames(stats, gsub(" ", "_", names(stats)))
})

# create a dataset with all  numeric statistics
firstStat = "Bif_ampl_local";
lastStat = "Convexity"
morf0 = lapply(stats, function(stats) setRownames(stats[match(firstStat, names(stats)):match(lastStat, names(stats))], 
                                                  make.unique(spaste(stats$`Brain`, ".", stats$`Reconstruction #`))));
numStats = lapply(morf0, dropConstantColumns);

nSets = length(stats);
nNumStats = sapply(numStats, ncol)
nNumStats # 31
numericCols = list(colnames(numStats[[1]]))

# define levels and ordering
order_subregion = c('CPr','CPi','CPc')
order_age = c('2m', '12m')
order_type = c('D1', 'D2')
order_Brain = c('TME07-1', 'TME08-1', 'TME09-1', 'TME10-1', 'TME10-3', 
                'hTME15-1', 'hTME19-2', 'hTME16-1', 'hTME24-2')
order_Brain_ind = c('TME07-1 (2m,WT)', 'TME08-1 (2m,WT)', 'TME09-1 (2m,WT)', 'TME10-1 (2m,WT)', 'TME10-3 (2m,WT)',
                    'hTME15-1 (12m,WT)', 'hTME19-2 (12m,WT)', 'hTME16-1 (12m,WT)', 'hTME24-2 (12m,WT)')

stats[[1]]$Striatal.Subregion = factor(stats[[1]]$Striatal.Subregion, levels=order_subregion)
stats[[1]]$Age = factor(stats[[1]]$Age, levels=order_age)
stats[[1]]$Type = factor(stats[[1]]$Type, levels=order_type)
stats[[1]]$Brain = factor(stats[[1]]$Brain, levels=order_Brain)
stats[[1]]$Brain_ind = factor(stats[[1]]$Brain_ind, levels=order_Brain_ind)

stats1$Striatal.Subregion = factor(stats1$Striatal.Subregion, levels=order_subregion)
stats1$Age = factor(stats1$Age, levels=order_age)
stats1$Type = factor(stats1$Type, levels=order_type)
stats1$Brain = factor(stats1$Brain, levels=order_Brain)
stats1$Brain_ind = factor(stats1$Brain_ind, levels=order_Brain_ind)

# whether or not the variances are equal for groups in the following anova test
var_equal_global = FALSE
p_correction = TRUE 

colnames(stats[[1]])
non_numericCols = colnames(stats[[1]][, c(1:17)])

###############################################################################################
############################ 1. Frequency Distribution ########################################
###############################################################################################
pdf('Plots/1_Frequency_Distribution.pdf')
  # genotype distribution
  grid.table(data.frame(table(stats1$Age)) %>% rename(Age = Var1, Frequency = Freq), rows = NULL)
  plot.new()
  # Type distribution
  grid.table(data.frame(table(stats1$Type)) %>% rename(Type = Var1, Frequency = Freq), rows = NULL)
  plot.new()
  # genotype and type distribution
  grid.table(data.frame(table(stats1$Age, stats1$Type)) %>% rename(Age = Var1, Type = Var2, Frequency = Freq), rows = NULL)
  plot.new()
  # brain distribution
  grid.table(data.frame(table(stats1$Brain)) %>% rename(Brain = Var1, Frequency = Freq), rows = NULL)
  plot.new()
  # genotype, brain distribution
  grid.table(data.frame(table(stats1$Brain,stats1$Age))
             %>% rename(Brain = Var1, Cell.Type = Var2, Frequency = Freq) %>% arrange(Brain), rows = NULL)
  plot.new()
  # level distribution
  grid.table(data.frame(table(stats1$Striatal.Subregion)) %>% rename(Sub.Region = Var1, Frequency = Freq), rows = NULL)
dev.off()

###############################################################################################
################### 2. Distribution plot of D1 and D2 in each subregion #######################
###############################################################################################
df_summary = stats1 %>%
  group_by(Age, Striatal.Subregion, Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Age, Striatal.Subregion) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup() 
total_cnt = df_summary %>% 
  group_by(Age, Striatal.Subregion)%>% 
  summarise(Total = sum(Count), .groups = 'drop')
df_summary = total_cnt %>% 
  left_join(df_summary, by = c("Age", "Striatal.Subregion"))

pdf(paste0("Plots/2_Type_Distribution_in_Region_by_Age.pdf"), wi = 7, he = 6)
  p = ggplot(df_summary, aes(x = Striatal.Subregion, y = Count, fill = Type, label = sprintf("%.1f%%", Percentage))) +
    geom_bar(stat = "identity", position = "stack") + 
    geom_text(position = position_stack(vjust = 0.5), size = 3, color = "white") + 
    geom_text(aes(y = Total, label=Total), vjust = -0.5, size = 3.5, color = "black") + 
    facet_wrap(~Age) + 
    scale_fill_manual(values = c("D1" = "darkgreen", "D2" = "darkred")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() + 
    labs(title = "Distribution of Type across Striatal Subregions", 
         x = "Striatal Subregion", 
         y = "Count") +
    theme(
      plot.title = element_text(size = 16, face = 'bold', hjust=0.5), 
      axis.title = element_text(size = 14), #
      strip.text = element_text(size = 14), 
      axis.text = element_text(size = 12), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),  
      panel.background = element_blank(),  
      axis.line = element_line(color = "black") 
    )
  print(p)
dev.off()

###############################################################################################
############################ 3. histograms of each variable#######################################
###############################################################################################
pdf(file = "Plots/3_statHistograms.pdf", wi = 5, he = 3);
scpp(2);
for (set in 1:nSets) for (col in 1:nNumStats[set])
  hist(numStats[[set]] [, col], breaks = 20, main = numericCols[[set]][col], 
       xlab = numericCols[[set]][col], ylab = "Frequency");

dev.off();

###############################################################################################
######## 4. correlation heatmap of variables ##################################################
###############################################################################################
pdf(file = "Plots/4_correlationHeatmapOfStats.pdf", wi = 18, he = 13);
for (set in 1:nSets)
  corHeatmapWithDendro(bicor(as.matrix(numStats[[set]]), use = 'p', maxPOutliers = 0.01),
                       mar.main = c(10, 13, 2, 1), main = spaste("Correlations of shape statistics", statNames.pretty[set]),
                       dendroWidth = 2/13.5)
dev.off()

###############################################################################################
######## 5. hierarchical clustering of variables############
###############################################################################################
# Collapse the various stats by correlation
numStats.scaled = lapply(numStats, scale);
tree = lapply(numStats, function(x) hclust(as.dist(1-bicor(x, use = 'p', maxPOutliers = 0.05)), method = "average"))
height = 0.15;
clusters = lapply(tree, cutree, h = height); 

pdf(file = "Plots/5_statClusteringTree.pdf", wi = 8, he = 5);
par(mar = c(4, 1, 2, 13));
  plotDendrogram(tree[[set]], horiz = TRUE, main = spaste("Clustering of neuron shape statistics", statNames.pretty[set]),
                 sub = "", xlab = "");
  abline(v = height, col = "red")
dev.off(); 

###############################################################################################
#################### Create Subsets #####################################################
###############################################################################################
# create adjusted subsets
colnames(stats[[1]])
n_non_numeric_col = 18 # number of non-numeric columns to exclude when doing cbind

# 2. for D1 only
stats[[2]] = subset(stats[[1]], Type == 'D1')
numStats[[2]] = numStats[[1]][stats[[1]]$Type == 'D1',]

# 3. for D2 only
stats[[3]] = subset(stats[[1]], Type == 'D2')
numStats[[3]] = numStats[[1]][stats[[1]]$Type == 'D2',]

# 4. for 2-m only
stats[[4]] = subset(stats[[1]], Age == '2m')
numStats[[4]] = numStats[[1]][stats[[1]]$Age == '2m',]

# 5. for 12-m only
stats[[5]] = subset(stats[[1]], Age == '12m')
numStats[[5]] = numStats[[1]][stats[[1]]$Age == '12m',]

# 6. adjusted for subregion
numStats[[6]] = as.data.frame(empiricalBayesLM(numStats[[1]], removedCovariates = stats[[1]][["Striatal.Subregion"]],
                                               getEBadjustedData = FALSE)$adjustedData.OLS)
stats[[6]] = cbind(stats[[1]][1:n_non_numeric_col], numStats[[6]])

statNames = c('All', # 1
              'D1', # 2
              'D2', # 3
              '2-m', #4
              '12-m', #5
              'adj.For.Subregion' #6
)
names(stats) = names(numStats) = statNames;
nSets=length(stats) 
print(nSets) # 6

###############################################################################################
############################ Output CSV for association ##################################################
###############################################################################################
plotTraits_association = c("Age", 
                           # "Striatal.Subregion", 
                           spaste("Age_in_", order_subregion),
                           spaste("Age_in_", order_type),
                           as.vector(outer(spaste("Age_in_", order_type), spaste('_in_', order_subregion), FUN = paste0))
)

# Function to convert specified columns to factors with common levels
convert_factors = function(df, cols, levels) {
  df[cols] = lapply(df[cols], factor, levels = levels)
  return(df)
}
# Apply the function to each dataframe in the list
stats = lapply(stats, convert_factors, cols = plotTraits_association[!plotTraits_association == 'Striatal.Subregion'], levels = order_age)
stats = lapply(stats, convert_factors, cols = c('Striatal.Subregion'), levels = order_subregion)

# create a new loop to generate results.csv based on anova test
total_num_trials = nNumStats 

nSets = 1 

test_list = list()
for (set in 1:nSets){
  data_list = list()
  # for each sub-group 
  for(cat in plotTraits_association){
    
    # create placeholders for shape statistics and test statistics
    morphometric_features = c()
    statistic_list = c()
    p_value_list = c()
    p_value_adjusted_list = c()
    
    # for each shape statistic
    for(col in colnames(numStats[[set]])){
      morphometric_features = c(morphometric_features,col)
      
      # use anova test
      formula = as.formula(paste(col, "~ Age"))
      test_result =  t.test(formula, data = stats[[set]])
      statistic = test_result$statistic
      p_value = test_result$p.value
      p_value_adjusted = p.adjust(p_value, method = "bonferroni", n = total_num_trials)
      
      statistic_list = c(statistic_list, statistic)
      p_value_list = c(p_value_list, p_value)
      p_value_adjusted_list = c(p_value_adjusted_list, p_value_adjusted)
    }
    
    # convert to df
    test_df_sub = data.frame(morphometric_features, statistic_list, p_value_list, p_value_adjusted_list)
    colnames(test_df_sub) = c('ShapeParameter', paste0('statistic for ', cat), paste0('p for ', cat), paste0('p.adjusted for ', cat))
    data_list = append(data_list, list(test_df_sub))
  }
  # generate a big sheet for each set from data_list
  test_list[[set]] = Reduce(function(x, y) merge(x, y, by = "ShapeParameter", all.x = TRUE), data_list)
  # rename rownames
  rownames(test_list[[set]]) = test_list[[set]]$ShapeParameter
  
  write.csv(test_list[[set]], spaste("Results/associationOfShapeParameters-", statNames[set], ".csv"), row.names = FALSE)
}

###############################################################################################
######### 6. Heatmap of Normalized features across Striatal.Subregion ########################
###############################################################################################
# normalize features to its median by brain
normalized_df = stats[[1]] %>% group_by(Brain_ind) %>%
  mutate(across(all_of(numericCols[[1]]), ~ . / median(., na.rm = TRUE)))

# calculate mean in region by brain
mean_df = normalized_df %>% group_by(Brain_ind, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup()

CPr_CPi = mean_df %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPc = mean_df %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPi_CPc = mean_df %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))

# calculate difference in means by region and by brain, for CPr vs. CPi, CPr vs. CPc, CPi vs. CPc
calculate_differences = function(data, group1, group2) {
  data_group1 = data %>% filter(Striatal.Subregion == group1) %>% select(-Striatal.Subregion)
  data_group2 = data %>% filter(Striatal.Subregion == group2) %>% select(-Striatal.Subregion)
  diff_df = data_group1
  diff_df[,-1] = data_group1[,-1] - data_group2[,-1] 
  diff_df
}
CPr_CPi_diff = calculate_differences(mean_df, "CPr", "CPi")
CPr_CPc_diff = calculate_differences(mean_df, "CPr", "CPc")
CPi_CPc_diff = calculate_differences(mean_df, "CPi", "CPc")

# reshape data for heatmap
reshape_for_heatmap = function(data) {
  data = data %>%
    pivot_longer(cols = -Brain_ind, names_to = "Feature", values_to = "Value") %>%
    pivot_wider(names_from = Brain_ind, values_from = Value) %>% as.data.frame()
  
  rownames(data) = data$Feature
  data = data[, -1]
}

CPr_CPi_heatmap_data = reshape_for_heatmap(CPr_CPi_diff)
CPr_CPc_heatmap_data = reshape_for_heatmap(CPr_CPc_diff)
CPi_CPc_heatmap_data = reshape_for_heatmap(CPi_CPc_diff)

# create a combined data for uniform color bar range
min_color = min(rbind(CPr_CPi_heatmap_data, CPr_CPc_heatmap_data, CPi_CPc_heatmap_data))
max_color = max(rbind(CPr_CPi_heatmap_data, CPr_CPc_heatmap_data, CPi_CPc_heatmap_data))

# prepare subset data for pairwise subregion comparison within each brain
CPr_CPi_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPc_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPi_CPc_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))

# reorder heatmap rows by 3 main groups 
# define 3 groups of features and assign colors
features_angle = c('Bif_ampl_local', 'Bif_ampl_remote', 'Bif_tilt_local', 'Bif_tilt_remote', 'Bif_torque_local', 'Bif_torque_remote', 'Centripetal_Bias')
features_length = c('Depth', 'Height', 'Width', 'Length', 'Sum_EucDistance', 'Sum_PathDistance', 'Max_EucDistance', 'Max_PathDistance', 'ABEL_All', 'ABEL_Internal', 'ABEL_Terminal', 'BAPL_All', 'BAPL_Internal', 'BAPL_Terminal')
features_complexity = c('N_bifs', 'N_branch', 'N_tips', 'N_stems', 'Branch_Order', 'Fractal_Dim', 'Partition_asymmetry', 'Terminal_degree', 'Balancing_Factor', 'Convexity')

feature_reorder = data.frame(Feature = c(features_angle, features_length, features_complexity),
                             group = rep(1:3, times = c(length(features_angle), length(features_length), length(features_complexity))))

# function to generate the heatmap
draw_heatmap = function(data_for_heatmap, data_subset, group_name) {
  
  # generate p-value matrix for pairwise subregion comparison per brain
  p_values = list()
  for(Brain_ind in order_Brain_ind) {
    brain_data = data_subset[data_subset$Brain_ind == Brain_ind, ]
    unique_region = unique(as.character(data_subset$Striatal.Subregion))
    
    # Loop through each feature and perform a t-test
    p_values[[Brain_ind]] = sapply(numericCols[[1]], function(feature) {
      test =  wilcox.test(brain_data[brain_data$Striatal.Subregion == unique_region[1], feature],
                          brain_data[brain_data$Striatal.Subregion == unique_region[2], feature],
                          exact = FALSE)
      test$p.value
    })
  }
  
  adjusted_p_values = lapply(p_values, function(p) p.adjust(p, method = "bonferroni", n = nNumStats))
  adjusted_p_values_df = as.data.frame(adjusted_p_values)
  rownames(adjusted_p_values_df) = numericCols[[1]]
  colnames(adjusted_p_values_df) = order_Brain_ind
  significance_df = apply(adjusted_p_values_df, 2, function(x) {
    ifelse(x < 0.001, "***",
           ifelse(x < 0.01, "**",
                  ifelse(x < 0.05, "*", "")))
  })
  
  data_for_heatmap = round(data_for_heatmap, 2)
  # sort data matrix by re-ordered features
  data_for_heatmap = data_for_heatmap[match(feature_reorder$Feature, rownames(data_for_heatmap)), ]
  # sort p_significance matrix by re-ordered features
  significance_df = significance_df[match(feature_reorder$Feature, rownames(significance_df)), ]
  
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
  blue_white_red = colorRamp2(c(min_color, 0, max_color), 
                              c("blue", "white", "red"))
  
  # create data matrix for text display
  text_matrix = matrix(nrow = nrow(data_for_heatmap), ncol = ncol(data_for_heatmap))
  for(i in 1:nrow(data_for_heatmap)) {
    for(j in 1:ncol(data_for_heatmap)) {
      text_matrix[i, j] = paste0(data_for_heatmap[i, j], significance_df[i, j])
    }
  }
  
  # generate heatmap
  ht = Heatmap(as.matrix(data_for_heatmap), 
               name = " ", 
               col = blue_white_red,
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               # cluster_columns = dend_reordered, 
               show_row_names = TRUE, row_title = 'Morphometrics',
               show_column_names = TRUE, column_title = 'Brains',
               width = unit(4, 'inch'), 
               height = unit(5, 'inch'),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 # Use text from text_matrix
                 grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 8))
               },
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
               heatmap_legend_param = list(legend_direction = "horizontal",
                                           title_position = 'topcenter'),
               right_annotation = right_annotation
  ) 
  draw(ht, heatmap_legend_side = 'top')
  grid.text(paste0('Difference of means of\nnormalized morphometrics: ', group_name),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
}

# Draw heatmaps
pdf(file = "Plots/6_heatmap - normalized_feature_region_pairwise_comparison.pdf", wi = 7, he = 8)
scpp(2.5);
par(mar = c(3,5,3,5));

  draw_heatmap(CPr_CPi_heatmap_data, CPr_CPi_sub, "CPr vs. CPi")
  draw_heatmap(CPr_CPc_heatmap_data, CPr_CPc_sub, "CPr vs. CPc")
  draw_heatmap(CPi_CPc_heatmap_data, CPi_CPc_sub, "CPi vs. CPc")

dev.off()

###############################################################################################
####################### 7. Above heatmap further divided by D1D2 ##############################
###############################################################################################
# calculate mean in region by brain
mean_df_d1d2 = normalized_df %>% group_by(Brain_ind, Type, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup()

CPr_CPi_d1d2 = mean_df_d1d2 %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPc_d1d2 = mean_df_d1d2 %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPi_CPc_d1d2 = mean_df_d1d2 %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))

# calculate difference in means by region and by brain, for CPr vs. CPi, CPr vs. CPc, CPi vs. CPc
calculate_differences_d1d2 = function(data, group1, group2) {
  data_group1 = data %>% filter(Striatal.Subregion == group1) %>% select(-Striatal.Subregion)
  data_group2 = data %>% filter(Striatal.Subregion == group2) %>% select(-Striatal.Subregion)
  diff_df = data_group1
  diff_df[,-c(1,2)] = data_group1[,-c(1,2)] - data_group2[,-c(1,2)] 
  diff_df
}
CPr_CPi_diff_d1d2 = calculate_differences_d1d2(mean_df_d1d2, "CPr", "CPi")
CPr_CPc_diff_d1d2 = calculate_differences_d1d2(mean_df_d1d2, "CPr", "CPc")
CPi_CPc_diff_d1d2 = calculate_differences_d1d2(mean_df_d1d2, "CPi", "CPc")

# reshape data for heatmap
reshape_for_heatmap_d1d2 = function(data) {
  data = data %>%
    pivot_longer(cols = -c(Brain_ind,Type), names_to = "Feature", values_to = "Value") %>%
    pivot_wider(names_from = c(Brain_ind,Type), values_from = Value) %>% as.data.frame()
  
  rownames(data) = data$Feature
  data = data[, -1]
}

CPr_CPi_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPr_CPi_diff_d1d2)
CPr_CPc_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPr_CPc_diff_d1d2)
CPi_CPc_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPi_CPc_diff_d1d2)

# prepare subset data for pairwise subregion comparison within each brain
CPr_CPi_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPi_sub$brain_type = paste0(CPr_CPi_sub$Brain, '_', CPr_CPi_sub$Type)
CPr_CPc_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPr_CPc_sub$brain_type = paste0(CPr_CPc_sub$Brain, '_', CPr_CPc_sub$Type)
CPi_CPc_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))
CPi_CPc_sub$brain_type = paste0(CPi_CPc_sub$Brain, '_', CPi_CPc_sub$Type)

# function to create the heatmap
draw_heatmap_d1d2 = function(data_for_heatmap, data_subset, group_name) {
  
  # create brain_type_age variable
  data_subset$brain_type_age = paste(data_subset$Brain_ind, data_subset$Type, sep = " (") %>%
    sapply(function(x) paste0(x, ")"))
  # re-order brain_type_age by type and genotype
  unique_brain_type_age = unique(data_subset$brain_type_age)
  unique_brain_type_age = unique_brain_type_age[order(-grepl("\\(2m,WT\\)", unique_brain_type_age), unique_brain_type_age)]
  
  # generate p-value matrix for pairwise subregion comparison per brain
  p_values = list()
  for(brain in unique_brain_type_age) {
    brain_data = data_subset[data_subset$brain_type_age == brain, ]
    unique_region = unique(as.character(data_subset$Striatal.Subregion))
    
    # Loop through each feature and perform a t-test
    p_values[[brain]] = sapply(numericCols[[1]], function(feature) {
      test =  wilcox.test(brain_data[brain_data$Striatal.Subregion == unique_region[1], feature],
                          brain_data[brain_data$Striatal.Subregion == unique_region[2], feature],
                          exact = FALSE)
      test$p.value
    })
  }
  
  adjusted_p_values = lapply(p_values, function(p) p.adjust(p, method = "bonferroni", n = nNumStats))
  adjusted_p_values_df = as.data.frame(adjusted_p_values)
  rownames(adjusted_p_values_df) = numericCols[[1]]
  colnames(adjusted_p_values_df) = unique(data_subset$brain_type)
  significance_df = apply(adjusted_p_values_df, 2, function(x) {
    ifelse(x < 0.001, "***",
           ifelse(x < 0.01, "**",
                  ifelse(x < 0.05, "*", "")))
  })
  
  data_for_heatmap = round(data_for_heatmap, 2)
  # sort data matrix by re-ordered features
  data_for_heatmap = data_for_heatmap[match(feature_reorder$Feature, rownames(data_for_heatmap)), ]
  # sort p_significance matrix by re-ordered features
  significance_df = significance_df[match(feature_reorder$Feature, rownames(significance_df)), ]
  
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
  blue_white_red = colorRamp2(c(min_color, 0, max_color), 
                              c("blue", "white", "red"))
  
  # create data matrix for text display
  text_matrix = matrix(nrow = nrow(data_for_heatmap), ncol = ncol(data_for_heatmap))
  for(i in 1:nrow(data_for_heatmap)) {
    for(j in 1:ncol(data_for_heatmap)) {
      text_matrix[i, j] = paste0(data_for_heatmap[i, j], significance_df[i, j])
    }
  }
  
  # generate heatmap
  ht = Heatmap(as.matrix(data_for_heatmap), 
               name = " ", 
               col = blue_white_red,
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE, row_title = 'Morphometrics',
               show_column_names = TRUE, column_title = 'Brains',
               width = unit(6, 'inch'), 
               height = unit(6, 'inch'),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 # Use text from text_matrix
                 grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 8))
               },
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
               heatmap_legend_param = list(legend_direction = "horizontal",
                                           title_position = 'topcenter'),
               right_annotation = right_annotation
  ) 
  draw(ht, heatmap_legend_side = 'top')
  grid.text(paste0('Difference of means of\nnormalized morphometrics: ', group_name),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
}

# Draw heatmaps
pdf(file = "Plots/7_heatmap - normalized_feature_region_pairwise_comparison_byD1D2.pdf", wi = 10, he = 9)
scpp(2.5);
par(mar = c(3,5,3,5));

  draw_heatmap_d1d2(CPr_CPi_heatmap_data_d1d2, CPr_CPi_sub, "CPr vs. CPi")
  draw_heatmap_d1d2(CPr_CPc_heatmap_data_d1d2, CPr_CPc_sub, "CPr vs. CPc")
  draw_heatmap_d1d2(CPi_CPc_heatmap_data_d1d2, CPi_CPc_sub, "CPi vs. CPc")

dev.off()

##############################################################################################
############################ 8. Distribution of Brain in CPr, CPi and CPc ########################
###############################################################################################
stacked_bar_distribution = function(data, by_var, stack_var, set_name, scale_fill_color){
  # overall proportion
  overall_proportions = data %>%
    count(!!sym(stack_var)) %>%
    mutate(!!sym(by_var) := "All", 
           Proportion = n / sum(n))
  
  # cluster proportions
  cluster_totals = data %>%
    group_by(!!sym(by_var)) %>%
    summarise(total_per_cluster = n(), .groups = 'drop')
  
  cluster_proportions = data %>%
    group_by(!!sym(by_var), !!sym(stack_var)) %>%
    summarise(n = n(), .groups = 'drop')
  
  cluster_proportions = cluster_proportions %>%
    left_join(cluster_totals, by = by_var) %>%
    mutate(Proportion = n / total_per_cluster)
  if (!is.factor(cluster_proportions[, by_var])) {
    cluster_proportions[, by_var] = cluster_proportions[, by_var]
  }
  
  combined_data = bind_rows(
    overall_proportions,
    cluster_proportions
  ) %>%
    arrange(!!sym(by_var))
  
  # place 'All' group at the first
  combined_data[, by_var] = factor(combined_data[, by_var], levels = c('All', names(table(data[, by_var]))))
  
  # group counts
  grp_cnts = combined_data %>% group_by(!!sym(by_var)) %>% summarize(n = sum(n))
  grp_cnt_strings = paste(grp_cnts[[by_var]], " (", grp_cnts$n, ")", sep = "")
  
  # chi-square
  pval_chi = c("") 
  expected_prop = prop.table(table(data[, stack_var]))
  for(grp in names(table(data[, by_var]))){
    ct = chisq.test(table(data[data[, by_var] == grp, stack_var]), p = expected_prop)
    print(paste0('group', grp, ' - ', ct$p.value))
    pval_chi = append(pval_chi, ifelse(ct$p.value <= 0.05, '*', ''))
  }
  asterisks = rep(pval_chi, each=length(table(data[, stack_var])))

  combined_data = combined_data[order(combined_data[, by_var]), ]
  combined_data$Asterisks = asterisks
  combined_data$Proportion = round(combined_data$Proportion,3)
  
  # Create the stacked bar plot, display percentage on it
  p = ggplot(combined_data, aes(x = get(by_var), y = Proportion, fill = get(stack_var), label = Asterisks)) +
    geom_bar(stat = "identity", width = 0.5) +
    (function() {
      if(is.null(scale_fill_color)){
        return(scale_fill_brewer(palette = 'Set1'))
      } else {
        return(scale_fill_manual(values = scale_fill_color))
      }
    })() +
    labs(
      title = paste0("Proportion of ", stack_var, " by ", by_var, "\n(", set_name, ")"),
      x = by_var,
      y = "Proportion",
      fill = stack_var
    ) +
    scale_x_discrete(labels = grp_cnt_strings)+
    # Add text labels on the bars
    geom_text(aes(label = scales::percent(Proportion)), 
              position = position_stack(vjust = 0.5), 
              size = 3, 
              color = "black") +
    geom_text(aes(y = 1, label = Asterisks), vjust = 0, color = "red") +    
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() + 
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
  print(p)
}

pdf(paste0("Plots/8_Distribution_Plots_brain_level.pdf"), wi = 7, he = 7)
  stacked_bar_distribution(stats[[1]], 'Striatal.Subregion', 'Brain_ind', 'All', NULL)
dev.off()

###############################################################################################
##########  9. Boxplot distribution by subregions  ##############################
########## (2-m, 2-m D1, 2-m D2, 12-m, 12-m D1, 12-m D2)  ######################################
###############################################################################################
# Function to perform overall mixed model fitting and return significance of P-values ofdifferent morphometric features vs. Striatal.Subregion
get_overall_pvalue = function(data, cat, variate) {
  
  formula = as.formula(paste(cat, " ~ ", variate))
  test_result = aov(formula, data)
  test_summary = summary(test_result)
  p_value = test_summary[[1]]["Pr(>F)"][[1]][1]
  
  if(p_correction){
    p_value = p.adjust(p_value, method = "bonferroni", n = nNumStats)
  }
  
  if(p_value <= 0.001){
    '***'
  }else if(p_value > 0.001 && p_value <= 0.01){
    '**'
  }else if(p_value > 0.01 && p_value <= 0.05){
    '*'
  }else if(p_value > 0.05){
    'NS.'
  }
}

# Function to calculate unique y positions for braces
calculate_y_positions = function(test_result, cat, data_combined, y_limits) {
  Q1 = quantile(data_combined[, cat], 0.25)
  significant_pairs = subset(test_result, p.value < 0.05)[, 'contrast']
  y_positions = length(significant_pairs)
  if (y_positions > 0) {
    y_positions = seq(from = 1.02 * y_limits[4], 
                      by = (y_limits[2] - 1.02 * y_limits[4])*0.2, length.out = y_positions)
  }
  return(y_positions)
}

# Function to define outliers of a column and return lower&upper y limits exlcuding the outliers
get_y_lims = function(data, cat){
  # Calculate quartiles and IQR
  Q1 = quantile(data[, cat], 0.25)
  Q3 = quantile(data[, cat], 0.75)
  IQR = Q3 - Q1
  
  lower_limit_IQR = max(Q1 - 1.5 * IQR, min(data[, cat]))
  upper_limit_IQR = min(Q3 + 1.5 * IQR, max(data[, cat]))
  
  lower_limit = min(data[, cat], na.rm = TRUE)
  upper_limit = max(data[, cat], na.rm = TRUE) + Q1
  
  return(c(lower_limit, upper_limit, lower_limit_IQR, upper_limit_IQR))
}

# Function to generate level comparison plots with (WT and by D1/D2, Q140 and by D1/D2)
level_comparison = function(subset_list, 
                            subset_name_list,
                            set_name){
  
  data_combined = do.call(rbind, subset_list)
  data_combined$Subset = rep(subset_name_list, 
                             times = unlist(lapply(subset_list, nrow)))
  # counts of each subset
  subset_counts = data_combined %>% 
    group_by(Subset) %>% 
    summarise(Count = n())
  data_combined = merge(data_combined, subset_counts, by = "Subset")

  # for each L-measure stat
  for(cat in names(numStats[[1]])){
    
    # set y-axis maximum
    y_limits = get_y_lims(data_combined, cat)
    
    # count of each subregion in each subset
    counts_data = data_combined %>%
      group_by(Subset, Striatal.Subregion) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(PosX = Striatal.Subregion, PosY = y_limits[1])
    
    # Get p-values for each subset
    p_values = sapply(split(data_combined, data_combined$Subset), function(x) get_overall_pvalue(x, cat, 'Striatal.Subregion'))
    
    # Prepare annotation data
    annotation_data = data.frame(
      Subset = names(p_values),
      # Label = sprintf("Overall p = %.1e", p_values),  # Format the p-values
      Label = paste0("Overall P value: ", p_values),
      x = 2,
      y = y_limits[2]
    )
    
    # pairwise-contrast
    results_list = lapply(split(data_combined, data_combined$Subset), function(subset_data){
      test_result = pairwise.t.test(subset_data[, cat], subset_data[, 'Striatal.Subregion'], p.adjust.method = "BH")
      p_values_df  = as.data.frame(as.table(test_result$p.value))
      p_values_df = na.omit(p_values_df)
      names(p_values_df) = c("Group1", "Group2", "p.value")
      p_values_df = p_values_df %>%
        mutate(contrast = paste(Group1, Group2, sep = " - "))
      
      significant_comparisons = subset(p_values_df, p.value < 0.05)
      significant_comparisons
    }
    )
    
    # Create the box plot
    p = ggplot(data_combined, aes(x = Striatal.Subregion, y = get(cat), fill = Striatal.Subregion)) +
      geom_boxplot(outlier.shape = NA, # Hides outliers for a cleaner look
                   color = "black",    
                   size = 0.5) +        
      facet_wrap(~ Subset, scales = "free") +
      ggtitle(paste0(cat, ' by Striatal.Subregion, (', set_name, ')')) + 
      ylab(cat) + 
      ylim(y_limits[1], y_limits[2]) +
      theme_minimal(base_size = 14) + 
      theme(
        legend.position = "bottom",    
        legend.title = element_blank(), 
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),         
        axis.title = element_text(size = rel(1.2)),        
        strip.text = element_text(size = rel(1.2)),        
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        plot.background = element_rect(fill = "white", color = NA) 
      ) +
      scale_fill_manual(values = c("gray30", "gray60", "gray90"))
    # Add count annotations
    p = p + geom_text(data = counts_data, aes(x = PosX, y = PosY, label = paste("n =", Count)), 
                      position = position_dodge(width = 0.75), 
                      vjust = 1.2, size = 3, color = "black")
    # Add the annotations to the plot (overall p-values for ANOVA)
    p = p + geom_text(data = annotation_data, aes(label = Label, x = x, y = y),
                      hjust = 0.5, vjust = 0,  
                      inherit.aes = FALSE, size = 3) 
    # Add pairwise test annotations
    for (subset_name in names(results_list)) {
      test_result = results_list[[subset_name]]
      y_positions = calculate_y_positions(test_result, cat, data_combined, y_limits)
      
      if (length(y_positions) > 0 && y_positions[1] != 0) {
        significant_subsets = subset(test_result, p.value < 0.05)
        significant_pairs = significant_subsets[, 'contrast']
        p_vals = significant_subsets[, 'p.value']
        
        for (i in 1:length(significant_pairs)) {
          pair = significant_pairs[i]
          group1 = unlist(strsplit(pair, " - "))[1]
          group2 = unlist(strsplit(pair, " - "))[2]
          
          p_value = p_vals[i]
          
          if(p_value > 0.01 && p_value <= 0.05){
            annotations = '*'
          }else if(p_value > 0.001 && p_value <= 0.01){
            annotations = '**'
          }else if (p_value <= 0.001){
            annotations = '***'
          }
          
          p = p + geom_signif(comparisons = list(c(group1, group2)), map_signif_level = FALSE,
                              y_position = y_positions[i], tip_length = 0.02, textsize = 5,
                              annotations = annotations,
                              data = subset(data_combined, Subset == subset_name))
        }
      }
    }
    
    print(p)
  }
}

subset_list = list(stats[[4]], subset(stats[[4]], Type == 'D1'), subset(stats[[4]], Type == 'D2'),
                   stats[[5]], subset(stats[[5]], Type == 'D1'), subset(stats[[5]], Type == 'D2'))
subset_name_list = c('2-m', '2-m,D1', '2-m,D2', '12-m', '12-m,D1', '12-m,D2')
set_name = 'All'

pdf(spaste("Plots/9_Boxplots-by-level.", set_name, ".pdf"), wi = 9, he = 9)
scpp(2.5);
par(mar = c(6.3, 3, 2.5, 1));
  level_comparison(subset_list, subset_name_list, set_name)
dev.off()

###############################################################################################
############################ 10. Blue & Red Heatmap of 12m vs. P56  ########################
###############################################################################################
# function to generate blue/red color heatmap for D1,D2 comparison
rowwise_color_gradient = function(coef_matrix) {
  color_matrix = matrix(NA, nrow = nrow(coef_matrix), ncol = ncol(coef_matrix),
                        dimnames = dimnames(coef_matrix))
  
  # Define color gradients
  blue_palette = colorRampPalette(c("#FFC2C2", "#FF6565"))
  red_palette = colorRampPalette(c("#C2C2FF", "#6565FF"))
  # Apply row-wise
  for (i in 1:nrow(coef_matrix)) {
    row_values = coef_matrix[i, ]
    max_abs_value = max(abs(row_values))
    # Normalize values in the row between 0 and 1 for both positive and negative
    positive_normalized = row_values[row_values > 0] / max_abs_value
    negative_normalized = -row_values[row_values < 0] / max_abs_value
    # Apply color mapping
    color_matrix[i, row_values > 0] = blue_palette(length(positive_normalized))[rank(positive_normalized)]
    color_matrix[i, row_values < 0] = red_palette(length(negative_normalized))[rank(negative_normalized)]
  }
  return(color_matrix)
}

# generate heatmap, using P56 as baseline
colored_heatmap_by_age = function(data, cat_features, set_name, cat_name){
  
  data_matrix = matrix(NA, nrow = length(numericCols[[1]]))
  text_matrix = matrix(NA, nrow = length(numericCols[[1]]), ncol = 1)
  coef_matrix = matrix(NA, nrow = length(numericCols[[1]]), ncol = 1)
  pval_matrix = matrix(NA, nrow = length(numericCols[[1]]), ncol = 1)
  
  for(cat in cat_features){
    results = data.frame(Feature = numericCols[[1]],
                         Coefficient = NA, 
                         P_Value = NA)
    
    for (feature in results$Feature) {
      
      sub_data = data[!is.na(data[, cat]), ]
      sub_data[, cat] = factor(sub_data[, cat], levels = order_age)
      column_means = colMeans(sub_data[sub_data$Age == '2m', numericCols[[1]]], na.rm = TRUE)
      results$Baseline_mean = column_means
      
      # t-test
      formula = as.formula(paste(feature, " ~ ", 'Age'))
      test_result = t.test(formula, sub_data)
      p_value = test_result$p.value
      coef = mean(sub_data[sub_data$Age == '12m', feature]) - mean(sub_data[sub_data$Age == '2m', feature])
      
      if(p_correction){
        p_value = p.adjust(p_value, method = 'bonferroni', n = nNumStats)
      }
      
      results[results$Feature == feature, "Coefficient"] = coef
      results[results$Feature == feature, "P_Value"] = p_value
      results$Fold_change = round(abs(results$Coefficient) / results$Baseline_mean,1) 
    }
    
    results$P_Value = signif(results$P_Value, 2)
    # account for P-value = 1
    results[results$P_Value == 1, 'P_Value'] = 0.95 
    
    results$neg_log10p = -log10(results$P_Value)
    results$Coefficient = signif(results$Coefficient,3)
    
    data_matrix = cbind(data_matrix, results$neg_log10p)
    # color_matrix = cbind(color_matrix, results$Color)
    text_matrix = cbind(text_matrix, results$Coefficient)
    coef_matrix = cbind(coef_matrix, results$Coefficient)
    pval_matrix = cbind(pval_matrix, results$P_Value)
  }
  data_matrix = data_matrix[,-1]
  # color_matrix = color_matrix[,-1]
  text_matrix = text_matrix[,-1]
  coef_matrix = coef_matrix[,-1]
  pval_matrix = pval_matrix[,-1]
  
  # add * based on -log10P, if p<0.05 then -log10P should > 1.301
  text_matrix = matrix(as.character(text_matrix), nrow=nrow(text_matrix))
  text_matrix[data_matrix > 1.301] = paste0(text_matrix[data_matrix > 1.301], '*')
  text_matrix[data_matrix > 2.301] = paste0(text_matrix[data_matrix > 2.301], '*')
  text_matrix[data_matrix > 3.301] = paste0(text_matrix[data_matrix > 3.301], '*')
  
  # set row and column names of the matrix
  colLabs = gsub("Age", "12-m vs. 2-m", cat_features)
  rowLabs = results$Feature
  rownames(data_matrix) = rownames(text_matrix) = rownames(coef_matrix) = rownames(pval_matrix) = rowLabs
  colnames(data_matrix) = colnames(text_matrix) = colnames(coef_matrix) = colnames(pval_matrix) = colLabs
  
  # save text_matrix
  text_only_asterisk = apply(text_matrix, c(1, 2), function(x) gsub("[^*]", "", x))
  write.csv(text_only_asterisk,  file = paste0('12m-vs-2m-', cat_name, '.csv'), row.names = TRUE)
  
  color_matrix = rowwise_color_gradient(coef_matrix)
  color_matrix[data_matrix<1.301 & coef_matrix > 0] = "#FFF0F0"
  color_matrix[data_matrix<1.301 & coef_matrix < 0] = "#F0F0FF"
      
  data_matrix = data_matrix[feature_reorder$Feature, ]
  color_matrix = color_matrix[feature_reorder$Feature, ]
  text_matrix = text_matrix[feature_reorder$Feature, ]
  coef_matrix = coef_matrix[feature_reorder$Feature, ]
  
  # to map cell colors using color_matrix
  custom_color_fun = function(value, matrix = coef_matrix, colors = color_matrix) {
    color_index = match(value, matrix)
    return(colors[color_index])
  }
  
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
  if(length(cat_features) <= 6){
    wi = 3
    hi = 6
  }else if(length(cat_features) > 6){
    wi = 7
    hi = 6
  }
  
  ht = Heatmap(coef_matrix, 
               col = custom_color_fun,  
               name = " ",  
               show_row_names = TRUE, 
               show_column_names = TRUE,  
               cluster_rows = FALSE,  
               cluster_columns = FALSE,  
               width = unit(wi, 'inch'), 
               height = unit(hi, 'inch'),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 9))
               },
               right_annotation = right_annotation,
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
               show_heatmap_legend = FALSE
  )
  
  draw(ht)
  grid.text(paste0("Coefficient and Significance Heatmap (12-m vs. 2-m), ", set_name),
            x = unit(0.5, "npc"),
            y = unit(0.98, "npc"),  
            gp = gpar(fontsize = 12, fontface = 'bold'))
  
  # Define the viewport to place the legend (adjust x and y as needed)
  viewport = viewport(x = 0.5, y = 0.94, width = 1, height = 0.5)
  pushViewport(viewport)
  
  # Add the custom legend
  grid.rect(gp = gpar(fill = "red"), x = 0.2, width = 0.06, height = 0.03)
  grid.text("12-m higher than 2-m", x = 0.24, just = "left", gp = gpar(fontsize = 8))
  grid.rect(gp = gpar(fill = "blue"), x = 0.5, width = 0.06, height = 0.03)
  grid.text("12-m lower than 2-m", x = 0.54, just = "left", gp = gpar(fontsize = 8))
  
  # Reset viewport
  popViewport()
}

pdf(file = spaste("Plots/10_Heatmap-Association-12m.vs.2m-Blue-Red-",'All', ".pdf"), wi = 9, he = 9);
scpp(2.5);
par(mar = c(3,5,3,5));
  # 1). All, by region
  colored_heatmap_by_age(stats[[1]], c('Age', 'Age_in_CPr', 'Age_in_CPi', 'Age_in_CPc'), 'All', 'Age_in_Region') 
  # 2) All, by D1D2
  colored_heatmap_by_age(stats[[1]], c('Age', 'Age_in_D1', 'Age_in_D2'), 'All', 'Age_in_Type') 
  # 3) All, by Type in Regions
  colored_heatmap_by_age(stats[[1]], c('Age', 
                                       'Age_in_D1_in_CPr', 'Age_in_D2_in_CPr',
                                       'Age_in_D1_in_CPi', 'Age_in_D2_in_CPi',
                                       'Age_in_D1_in_CPc', 'Age_in_D2_in_CPc'), 'All', 'Age_in_Type_in_Region') 
dev.off()

###############################################################################################
################# 11. Box Plot 12-m vs. 2-m Group Comparison (overall, by subregion) #########
###############################################################################################
comparison_var1_across_var2 = function(data, compare_var, across_var, set_name){
  
  all_data = data
  all_data[, across_var] = 'All'
  combined_data = rbind(data, all_data)
  combined_data$Age = factor(combined_data$Age, levels = order_age)
  if(across_var == 'Striatal.Subregion'){
    combined_data$Striatal.Subregion = factor(combined_data$Striatal.Subregion, levels = c('All', order_subregion))
  }
  
  for(cat in numericCols[[1]]){
    # count within each subgroup
    counts_data = combined_data %>%
      group_by(!!sym(across_var), !!sym(compare_var)) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(PosX = !!sym(compare_var), PosY = 0.98*min(combined_data[, cat], na.rm = TRUE))
    
    # Get p-values for each subset (between compare_var)
    p_values = sapply(split(combined_data, combined_data[, across_var]), function(x) get_overall_pvalue(x, cat, 'Age'))
    
    # Prepare annotation data
    annotation_data = data.frame(
      Striatal.Subregion = names(p_values),
      Label = paste0("Overall P: ", p_values),
      x = 1.5,
      y = max(combined_data[, cat], na.rm = TRUE)
    )
    
    if(across_var == 'Striatal.Subregion'){
      annotation_data$Striatal.Subregion = factor(annotation_data$Striatal.Subregion, levels = c('All', order_subregion))
    }
    
    p = ggplot(data = combined_data, aes(x = !!sym(compare_var), y = !!sym(cat), fill = !!sym(compare_var))) +
      geom_boxplot(position = position_dodge(width = 0.6),
                   outlier.shape = NA, 
                   size = 0.5) +
      facet_grid(. ~ get(across_var), scales = "free_x", space = "free_x") +
      theme_minimal() +
      theme(
        legend.position = "bottom",    
        legend.title = element_blank(), 
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text( hjust = 1),         
        axis.title = element_text(size = rel(1.2)),        
        strip.text = element_text(size = rel(1.2)),        
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        plot.background = element_rect(fill = "white", color = NA) 
      ) +
      labs(title = paste0(cat, " by ", compare_var, " and ", across_var, "\n", set_name),
           x = compare_var,
           y = cat) + 
      if(compare_var == 'Type'){
        scale_fill_manual(values = c("D1" = "#8B0000", "D2" = "darkgreen")) 
      }else if(compare_var == 'Genotype'){
        scale_fill_manual(values = c("WT" = "darkgreen", "Q140" = "red")) 
      }else{
        scale_fill_brewer(palette = "Dark2") 
      } 
    # Add count annotations
    p = p + geom_text(data = counts_data, aes(x = PosX, y = PosY, label = paste("n =", Count)), 
                      position = position_dodge(width = 0.75), 
                      vjust = 1.2, size = 3, color = "black")
    # Add the annotations to the plot (overall p-values for ANOVA)
    p = p + geom_text(data = annotation_data, aes(label = Label, x = x, y = y),
                      hjust = 0.5, vjust = 0,  # Adjust horizontal and vertical position
                      inherit.aes = FALSE, size = 3) 
    print(p)
  }
}

# 12-m/2-m comparison, across levels
pdf("Plots/11_Boxplots-Age.by.Striatal.Subregion.All.pdf", width = 7, height = 5)
scpp(2);
  comparison_var1_across_var2(stats[[1]], 'Age', 'Striatal.Subregion', 'All')
dev.off()

###############################################################################################
################# 12. Box Plot 12-m vs. 2-m  (overall, by level and D1D2) ##################
###############################################################################################
create_all_cat = function(x, across_var){
  all = copy(x)
  all[, across_var] = 'All'
  combined_data = rbind(x, all)
  combined_data$Age = factor(combined_data$Age, levels = order_age)
  if(across_var == 'Striatal.Subregion'){
    combined_data$Striatal.Subregion = factor(combined_data$Striatal.Subregion, levels = c('All', order_subregion))
  }
  return (combined_data)
}

comparison_var1_across_var2_all = function(subset_list,
                                           subset_name_list,
                                           compare_var,
                                           across_var){
  
  subset_list_with_all = lapply(subset_list, function(x)create_all_cat(x, across_var))
  data_combined = do.call(rbind, subset_list_with_all)
  data_combined$Subset = rep(subset_name_list, 
                             times = unlist(lapply(subset_list, function(x) 2*nrow(x))))
  data_combined[, across_var] = factor(data_combined[, across_var], levels = c('All', order_subregion))
  
  for(cat in numericCols[[1]]){
    # count within each subgroup
    counts_data = data_combined %>%
      group_by(Subset, !!sym(across_var), !!sym(compare_var)) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(PosX = !!sym(compare_var), PosY = 0.98*min(data_combined[, cat], na.rm = TRUE))
    
    # Get p-values for each subset (between compare_var)
    combined_factor = interaction(data_combined$Subset, data_combined[[across_var]])
    p_values = sapply(split(data_combined,combined_factor), function(x) get_overall_pvalue(x, cat, 'Age'))
    
    # Prepare annotation data
    annotation_data = data.frame(
      Subset_level = names(p_values),
      Label = paste0("Overall P: ", p_values),
      x = 1.5,
      y = max(data_combined[, cat], na.rm = TRUE)
    )
    annotation_data$Subset = sapply(strsplit(annotation_data$Subset_level, split = "\\."), `[`, 1)
    annotation_data[, across_var] = sapply(strsplit(annotation_data$Subset_level, split = "\\."), `[`, 2)
    annotation_data$Subset = factor(annotation_data$Subset, levels = subset_name_list)
    annotation_data[, across_var] = factor(annotation_data[, across_var], levels = c('All', order_subregion))
    
    p = ggplot(data = data_combined, aes(x = !!sym(compare_var), y = !!sym(cat), fill = !!sym(compare_var))) +
      geom_boxplot(position = position_dodge(width = 0.6),
                   outlier.shape = NA, 
                   size = 0.5) +
      facet_grid( Subset ~ get(across_var), scales = "free_x", space = "free_x") +
      theme_minimal() +
      theme(
        legend.position = "bottom",    
        legend.title = element_blank(), 
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text.x = element_text( hjust = 1),         
        axis.title = element_text(size = rel(1.2)),        
        strip.text = element_text(size = rel(1.2)),        
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        plot.background = element_rect(fill = "white", color = NA) 
      ) +
      labs(title = paste0(cat, " by ", compare_var, " and ", across_var),
           x = compare_var,
           y = cat) + 
      if(compare_var == 'Type'){
        scale_fill_manual(values = c("D1" = "#8B0000", "D2" = "#00008B")) 
      }else if(compare_var == 'Genotype'){
        scale_fill_manual(values = c("WT" = "darkgreen", "Q140" = "darkorange")) 
      }else{
        scale_fill_brewer(palette = "Dark2") 
      } 
    # Add count annotations
    p = p + geom_text(data = counts_data, aes(x = PosX, y = PosY, label = paste("n =", Count)), 
                      position = position_dodge(width = 0.75), 
                      vjust = 1.2, size = 3, color = "black")
    p = p + geom_text(data = annotation_data, aes(label = Label, x = x, y = y),
                      hjust = 0.5, vjust = 0,  # Adjust horizontal and vertical position
                      inherit.aes = FALSE, size = 3) 
    print(p)
  }
}


subset_list = list(stats[[1]], stats[[2]], stats[[3]]) # All, D1, D2
subset_name_list = c('All Neurons', 'D1', 'D2')

pdf(paste0("Plots/12_Boxplots-Age.by.Striatal.Subregion.byD1D2.pdf"), width = 7, height = 10)
scpp(2);
  comparison_var1_across_var2_all(subset_list, subset_name_list, 'Age', 'Striatal.Subregion')
dev.off()

###############################################################################################
##########  13. PCA, UMAP, cell clustering to compare 12m and 2m, All, D1 and D2 #############
###############################################################################################
# start cell clustering
pdf(file = spaste("Plots/13_Cell_Clustering_PCA_UMAP.pdf"), wi = 9, he = 9);
scpp(2.5);
par(mar = c(3,5,3,5));
  # 'All', 'D1', 'D2'
  for(set in c(1,2,3)){
    # 1) PCA
    pca_result = prcomp(scale(numStats[[set]]), center = FALSE, scale. = FALSE)
    explained_variance_ratio = round(pca_result$sdev^2 / sum(pca_result$sdev^2),2)
    # the n_PC to keep that explains > 90% variance
    n_keep = which(cumsum(explained_variance_ratio) > 0.90)[1]
    
    pc_data = data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], 
                         PC3 = pca_result$x[,3], PC4 = pca_result$x[,4], 
                         Age = stats[[set]]$Age)
    # PC1 vs. PC2
    p1 = ggplot(pc_data, aes(x = PC1, y = PC2, color = Age)) +
      geom_point() +
      theme_minimal() +
      ggtitle(paste0("PCA: PC1 vs PC2 colored by Age, ", names(stats)[[set]])) +
      xlab(paste0("PC 1 (explained variance: ", explained_variance_ratio[1], ')')) +
      ylab(paste0("PC 2 (explained variance: ", explained_variance_ratio[2], ')')) + 
      theme_minimal() + 
      theme(legend.position = "top",
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size = 14),  
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),  
            axis.text.y = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 
    print(p1)
    
    # PC3 vs. PC4
    p2 = ggplot(pc_data, aes(x = PC3, y = PC4, color = Age)) +
      geom_point() +
      theme_minimal() +
      ggtitle(paste0("PCA: PC3 vs PC4 colored by Age, ", names(stats)[[set]])) +
      xlab(paste0("PC 3 (explained variance: ", explained_variance_ratio[3], ')')) +
      ylab(paste0("PC 4 (explained variance: ", explained_variance_ratio[4], ')')) +
      theme_minimal() + 
      theme(legend.position = "top",
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size = 14),  
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),  
            axis.text.y = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(p2)
    
    # loadings for PC1
    loadings = pca_result$rotation
    loadings_pc1 = loadings[, 1]
    loadings_df = data.frame(Variable = rownames(loadings), Loading = loadings_pc1)
    loadings_df = loadings_df[order(abs(loadings_df$Loading), decreasing = TRUE), ]
    p3 = ggplot(loadings_df, aes(x = reorder(Variable, Loading), y = Loading, fill = Loading > 0)) +
      geom_bar(stat = "identity") +
      coord_flip() + #
      labs(x = "Variable", y = "Loading on PC1", title = paste0("Loadings for PC1, ", names(stats)[[set]])) +
      theme_minimal() +
      scale_fill_manual(name = "Direction", values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("Positive", "Negative")) +
      theme(legend.position = "top") +
      geom_hline(yintercept = -0.3, linetype = "dashed", color = "black") + 
      geom_hline(yintercept = 0.3, linetype = "dashed", color = "black")  
    print(p3)
    
    # 2) UMAP
    umap_result = umap(scale(numStats[[set]]))
    # Prepare the dataframe for plotting
    df_umap = as.data.frame(umap_result$layout)
    colnames(df_umap) = c("UMAP1", "UMAP2")
    df_umap$Age = stats[[set]]$Age
    
    # Plot UMAP1 vs UMAP2, colored by 'Age'
    p4 = ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = Age)) +
      geom_point() +
      theme_minimal() +
      ggtitle(paste0("UMAP Visualization Colored by Age, ", names(stats)[[set]])) + 
      theme_minimal() + 
      theme(legend.position = "top",
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size = 14),  
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 12),  
            axis.text.y = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    print(p4)
    
    # 3) Clustering
    cell_dist = list(dist(pca_result$x[, 1:n_keep]))
    cellTree = lapply(cell_dist, hclust, method = "a");
    # cluster table of all cells
    if(names(stats)[[set]] == 'D1'){
      deepsplit = 2.5
      minclustersize = 300
    }else if(names(stats)[[set]] == 'D2'){
      deepsplit = 3
      minclustersize = 200
    }else if(names(stats)[[set]] == 'All'){
      deepsplit = 3
      minclustersize = 400
    }
    cellClusters = mymapply(cutreeDynamic, cellTree, distM = lapply(cell_dist, as.matrix),
                            MoreArgs = list(deepSplit = deepsplit, 
                                            minClusterSize = minclustersize, cutHeight = NULL))
    print(table(cellClusters))
    sil_score = mean(silhouette(cellClusters[[1]], cell_dist[[1]])[, 'sil_width'])
    print(paste0("Silouette score: ", sil_score)) 
    
    cellStatColors = lapply(list(scale(numStats[[set]])), function(x) 
      numbers2colors(x, lim = c(-quantile(abs(x), prob = 0.99), quantile(abs(x), prob = 0.99)), commonLim = TRUE))
    
    # distribution of categorical variables 
    cellClass = lapply(list(stats[[set]]), function(st) cbind(st[c("Age")]));
    
    cellClassColorData = lapply(cellClass, function(x) dataFrame2colors(x[multiGrepv(c("Age"), names(x))]));
    
    plotDendroAndColors(cellTree[[1]],
                        cbind(labels2colors(cellClusters[[1]], colorSeq = c('brown','gold',"darkgreen", "navy", "purple")), 
                              cellClassColorData[[1]]$colors,
                              cellStatColors[[1]]),
                        c("Cell clusters", cellClassColorData[[1]]$legend, colnames(numStats[[set]])), 
                        rowText = cellClusters[[1]],
                        autoColorHeight = FALSE, colorHeight = 0.65,
                        textPositions = 1,
                        rowWidths = c(1, 2, rep(1, length(cellClassColorData[[1]]$legend)), rep(1, ncol(numStats[[set]]))),
                        rowTextAlignment = "center",
                        marAll = c(0.1, 22, 2, 0.2), 
                        main = spaste("Cell clustering, ", names(stats)[[set]]),
                        dendroLabels = FALSE,
                        addGuide = TRUE)
    
    cell_cluster_df = data.frame(Cluster = cellClusters[[1]])
    cell_cluster_df_result = cbind(stats[[set]], cell_cluster_df)
    cell_cluster_df_result$Cluster = as.character(cell_cluster_df_result$Cluster)
    table(cell_cluster_df_result$Cluster)
    
    # stacked Age in cell Clusters
    stacked_bar_distribution(subset(cell_cluster_df_result, Cluster != 0), 'Cluster', 'Age', names(stats)[[set]], c("2m" = "turquoise", "12m" = "blue"))
  
}
dev.off()

###############################################################################################
######### 14. radar plot show pct% change of mean 12-m compared to P56 in CPr/i/c #############
###############################################################################################
angle_features = c('Bif_ampl_local', 'Bif_ampl_remote', 'Bif_tilt_local', 'Bif_tilt_remote',
                   'Bif_torque_local', 'Bif_torque_remote', 'Centripetal_Bias')
size_features = c('Length', 'Sum_EucDistance', 'Sum_PathDistance', 'Max_PathDistance', 'Max_EucDistance',
                  'ABEL_All', 'BAPL_All', 'ABEL_Terminal', 'BAPL_Terminal', 'ABEL_Internal', 'BAPL_Internal',
                  'Height', 'Width', 'Depth')
complexity_features = c('N_stems', 'N_branch', 'N_bifs', 'N_tips', 'Terminal_degree',
                        'Branch_Order', 'Fractal_Dim', 'Partition_asymmetry',
                        'Balancing_Factor', 'Convexity')
all_features = c(angle_features, size_features, complexity_features)

percent_change = function(new, old) {
  pct_diff = (new - old) / old
  return (round(pct_diff * 100, 1))
}

create_dup = function(df){
  df_dup = copy(df)
  df_dup$Type = 'All'
  df_comb = rbind(df, df_dup)
  df_comb$Type = factor(df_comb$Type, levels = c('All', 'D1', 'D2'))
  return(df_comb)
}

stats_dup = create_dup(stats[[1]])

mean_df = stats_dup %>% group_by(Type, Age, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup() %>% as.data.frame()

pct_changes = mean_df %>% group_by(Type, Striatal.Subregion) %>%
  summarise(across(numericCols[[1]], ~ percent_change(.x[2], .x[1])))

print(min(pct_changes[, numericCols[[1]]]))
plot_min = -40
print(max(pct_changes[, numericCols[[1]]]))
plot_max = 20

pretty_num_col_names_for_radar = c('Bif_ampl\nlocal', 'Bif_ampl\nremote', 'Bif tilt\nlocal', 'Bif tilt\nremote', 
                                   'Bif torque\nlocal', 'Bif torque\nremote','Centripetal\nBias',
                                   'Length','Sum Euc\nDistance', 'Sum Path\nDistance', 
                                   'Max Euc\nDistance', 'Max Path\nDistance', 'ABEL\nAll', 'BAPL\nAll',
                                   'ABEL\nTerminal', 'BAPL\nTerminal', 'ABEL\nInternal', 'BAPL\nInternal',
                                   'Height', 'Width', 'Depth',  
                                   'N stems', 'N branch', 'N bifs', 'N tips',  'Terminal\ndegree',
                                   'Branch\nOrder', 'Fractal\nDim',
                                   'Partition\nasymmetry', 'Balancing\nFactor', 'Convexity'
)

create_radar_pct_change = function(df, grp, plot_min, plot_max){
  df = df[, all_features]
  # add max and min value
  df = rbind(rep(plot_max,length(numericCols[[1]])) , 
             rep(plot_min,length(numericCols[[1]])) , df)
  # prettify the feature names
  colnames(df) = pretty_num_col_names_for_radar
  # Define colors for All, D1 and D2
  colors_border = c('blue', 'green', 'red')
  
  # Create the radar chart
  radarchart(df, axistype = 4, seg= 6,
             # Custom polygon
             pcol = colors_border,plwd = 3, plty = 1,
             pfcol = NULL,
             
             # Custom the grid
             cglcol = "gray", cglty = 1, axislabcol = "black", 
             caxislabels = c('-40%', '-30%', '-20%\n', '-10%\n', '0\n', '10%\n', '20%\n'), 
             cglwd = 0.5,
             
             # Custom labels
             vlcex = 0.8)
  
  n_vertices = length(numericCols[[1]]) * 3
  angles = seq(0, 2 * pi, length.out = n_vertices)
  radius_0 = 0.71
  x_coords = radius_0 * cos(angles)
  y_coords = radius_0 * sin(angles)
  lines(x_coords, y_coords, col = "black", lwd = 4, lty = 2)
  
  title(main = paste0("Percent (%) change of means - 12-m vs. P56\n", grp), col.main = "black", font.main = 4)
  legend(x = 1.2, y = 1, legend = c('All', 'D1', 'D2'), bty = "n", pch = 20, col = colors_border, 
         text.col = "black", cex = 1.2, pt.cex = 3)
  
}

pdf(file = spaste("Plots/14_Radar_CPr.CPi.CPc_pct_change_12m_2m_byD1D2.pdf"), wi = 13, he = 10);
scpp(2.5);
par(mar = c(3,5,3,5));
  create_radar_pct_change(subset(pct_changes, Striatal.Subregion == 'CPr'), 'CPr', plot_min, plot_max)
  create_radar_pct_change(subset(pct_changes, Striatal.Subregion == 'CPi'), 'CPi', plot_min, plot_max)
  create_radar_pct_change(subset(pct_changes, Striatal.Subregion == 'CPc'), 'CPc', plot_min, plot_max)
dev.off()

###############################################################################################
##########  15. Same radar plot but all subregions combined ################################
###############################################################################################
mean_df = stats_dup %>% group_by(Type, Age) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup() %>% as.data.frame()

pct_changes = mean_df %>% group_by(Type) %>%
  summarise(across(numericCols[[1]], ~ percent_change(.x[2], .x[1])))

print(min(pct_changes[, numericCols[[1]]]))
plot_min = -40
print(max(pct_changes[, numericCols[[1]]]))
plot_max = 20

pdf(file = spaste("Plots/15_Radar_pct_change_12m_2m_byD1D2_all_regions.pdf"), wi = 13, he = 10);
scpp(2.5);
par(mar = c(3,5,3,5));
  create_radar_pct_change(pct_changes, 'All', plot_min, plot_max)
dev.off()

