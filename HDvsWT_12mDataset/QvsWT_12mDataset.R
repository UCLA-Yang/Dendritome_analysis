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
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
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
library(car)
library(RColorBrewer) 
library(ggtext)
library(randomForest)
library(randomGLM)
library(openxlsx);
library(stringr)
library(effsize)
library(nlme)
library(lme4)
library(emmeans)
library(data.table)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(wesanderson)
library(paletteer)
library(fmsb)

# create directions
dir.create("Plots", recursive = TRUE);
dir.create("Results", recursive = TRUE);

###############################################################################################
############################ read data #######################################
###############################################################################################
stats1 = read.csv('all_htme_brains_with_registration_1168CPneurons_onlyCP.csv',check.names = FALSE)
colnames(stats1)
dim(stats1) # 1167 neurons
table(stats1$Genotype) # 543 WT, 624 Q140

###############################################################################################
############################ Rename, Add columns, Process data #######################################
###############################################################################################
stats = list(stats1); # convert data to list
statNames = c("Q140_WT"); 
statNames.pretty = c("");
names(stats) = statNames;

# create covariate variables 
stats = lapply(stats, function(stats) 
{
  # genotpe in subregions (3 anatomical levels, CPr, CPi, CPc)
  rest1 = restrictVariableByCovariateLevels(stats$Genotype, stats$Striatal.Subregion, 
                                            varName = "Genotype", covarName = "", nameSep = " in ", check.names = FALSE);
  # genotype in type
  rest2 = restrictVariableByCovariateLevels(stats$Genotype, stats$Type, 
                                            varName = "Genotype", covarName = "", nameSep = " in ", check.names = FALSE);
  
  stats = data.frame.ncn(rest1, rest2, stats);
  out = setNames(stats, gsub(" ", "_", names(stats)))
});

stats = lapply(stats, function(stats) {
  # genotype in D1 in regions
  rest1 = restrictVariableByCovariateLevels(stats$`Genotype_in_D1`, stats$Striatal.Subregion, 
                                            varName = "Genotype_in_D1", covarName = "", nameSep = " in ", check.names = FALSE);
  # genotype in D2 in regions
  rest2 = restrictVariableByCovariateLevels(stats$`Genotype_in_D2`, stats$Striatal.Subregion, 
                                            varName = "Genotype_in_D2", covarName = "", nameSep = " in ", check.names = FALSE);
  stats = data.frame.ncn(rest1, rest2, stats);
  out = setNames(stats, gsub(" ", "_", names(stats)))
})

# create a dataset with all numeric statistics
firstStat = "Bif_ampl_local";
lastStat = "Convexity"
morf0 = lapply(stats, function(stats) setRownames(stats[match(firstStat, names(stats)):match(lastStat, names(stats))], 
                                                  make.unique(spaste(stats$`Brain`, ".", stats$`Reconstruction #`))));
numStats = lapply(morf0, dropConstantColumns);

nSets = length(stats);
numericCols = lapply(numStats, colnames);
nNumStats = sapply(numStats, ncol) 
nNumStats[[1]] # 31 numeric features

# define levels of genotype, subregion, and type
order_type = c('D1', 'D2')
order_subregion = c('CPr','CPi','CPc')
order_genotype = c('WT', 'Q140')

stats[[1]]$Striatal.Subregion = factor(stats[[1]]$Striatal.Subregion, levels=order_subregion)
stats[[1]]$Genotype = factor(stats[[1]]$Genotype, levels=order_genotype)

stats1$Striatal.Subregion = factor(stats1$Striatal.Subregion, levels=order_subregion)
stats1$Genotype = factor(stats1$Genotype, levels=order_genotype)

# define global variables
var_equal_global = TRUE
p_correction = TRUE

###############################################################################################
############################ 1. Frequency Distribution Table ##################################
###############################################################################################
pdf('Plots/1_Frequency_Distribution.pdf')
  # genotype distribution
  grid.table(data.frame(table(stats1$Genotype)) %>% rename(Genotype = Var1, Frequency = Freq), rows = NULL)
  plot.new()
  # Type distribution
  grid.table(data.frame(table(stats1$Type)) %>% rename(Type = Var1, Frequency = Freq), rows = NULL)
  plot.new()
  # genotype and type distribution
  grid.table(data.frame(table(stats1$Genotype, stats1$Type)) %>% rename(Genotype = Var1, Type = Var2, Frequency = Freq), rows = NULL)
  plot.new()
  # brain distribution
  grid.table(data.frame(table(stats1$Brain)) %>% rename(Brain = Var1, Frequency = Freq), rows = NULL)
  plot.new()
  # genotype, brain distribution
  grid.table(data.frame(table(stats1$Brain,stats1$Genotype))
             %>% rename(Brain = Var1, Cell.Type = Var2, Frequency = Freq) %>% arrange(Brain), rows = NULL)
  plot.new()
  # level distribution
  grid.table(data.frame(table(stats1$Striatal.Subregion)) %>% rename(Sub.Region = Var1, Frequency = Freq), rows = NULL)
dev.off()

###############################################################################################
############################ Distribution Plots ###############################################
###### 2. Distribution plot of D1 and D2 in each level/region, facet Q and WT ######################
###############################################################################################
df_summary = stats1 %>%
  group_by(Genotype, Striatal.Subregion, Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Genotype, Striatal.Subregion) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup() 
total_cnt = df_summary %>% 
  group_by(Genotype, Striatal.Subregion)%>% 
  summarise(Total = sum(Count), .groups = 'drop')
df_summary = total_cnt %>% 
  left_join(df_summary, by = c("Genotype", "Striatal.Subregion"))

pdf(paste0("Plots/2_Type_Distribution_in_Region_by_Genotype.pdf"), wi = 7, he = 6)
  p = ggplot(df_summary, aes(x = Striatal.Subregion, y = Count, fill = Type, label = sprintf("%.1f%%", Percentage))) +
    geom_bar(stat = "identity", position = "stack") + 
    geom_text(position = position_stack(vjust = 0.5), size = 3, color = "white") + 
    geom_text(aes(y = Total, label=Total), vjust = -0.5, size = 3.5, color = "black") + 
    facet_wrap(~Genotype) + 
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
######################### 4. Hierarchical clustering of variables ################################
###############################################################################################
# Collapse the various stats by correlation
numStats.scaled = lapply(numStats, scale);
tree = lapply(numStats, function(x) hclust(as.dist(1-bicor(x, use = 'p', maxPOutliers = 0.05)), method = "average"))
height = 0.15;
clusters = lapply(tree, cutree, h = height); 

pdf(file = "Plots/4_statClusteringTree.pdf", wi = 8, he = 5);
  par(mar = c(4, 1, 2, 13));
  plotDendrogram(tree[[set]], horiz = TRUE, main = spaste("Clustering of neuron shape statistics", statNames.pretty[set]),
                 sub = "", xlab = "");
  abline(v = height, col = "red")
dev.off(); 

###############################################################################################
############################ Create Subsets ###################################################
###############################################################################################
# create adjusted subsets
colnames(stats[[1]])
n_non_numeric_col = 16 # number of non-numeric columns to exclude when doing cbind

# 2. for D1 only
stats[[2]] = subset(stats[[1]], Type == 'D1')
numStats[[2]] = numStats[[1]][stats[[1]]$Type == 'D1',]

# 3. for D2 only
stats[[3]] = subset(stats[[1]], Type == 'D2')
numStats[[3]] = numStats[[1]][stats[[1]]$Type == 'D2',]

# 4. for WT only
stats[[4]] = subset(stats[[1]], Genotype == 'WT')
numStats[[4]] = numStats[[1]][stats[[1]]$Genotype == 'WT',]

# 5. for Q140 only
stats[[5]] = subset(stats[[1]], Genotype == 'Q140')
numStats[[5]] = numStats[[1]][stats[[1]]$Genotype == 'Q140',]

# 6. adjusted for subregion
numStats[[6]] = as.data.frame(empiricalBayesLM(numStats[[1]], removedCovariates = stats[[1]][["Striatal.Subregion"]],
                                               getEBadjustedData = FALSE)$adjustedData.OLS)
stats[[6]] = cbind(stats[[1]][1:n_non_numeric_col], numStats[[6]])

statNames = c('All', # 1
              'D1', # 2
              'D2', # 3
              'WT', #4
              'Q140', #5
              'adj.For.Subregion' #6
)
names(stats) = names(numStats) = statNames;
nSets=length(stats) 
print(nSets) # 6

###############################################################################################
############################ Output CSV for association #######################################
###############################################################################################
plotTraits_association = c("Genotype", 
                           "Striatal.Subregion", 
                           spaste("Genotype_in_", order_subregion),
                           spaste("Genotype_in_", order_type),
                           as.vector(outer(spaste("Genotype_in_", order_type), spaste('_in_', order_subregion), FUN = spaste))
)

# Function to convert specified columns to factors with common levels
convert_factors = function(df, cols, levels) {
  df[cols] = lapply(df[cols], factor, levels = levels)
  return(df)
}
# Apply the function to each dataframe in the list
stats = lapply(stats, convert_factors, cols = plotTraits_association[!plotTraits_association == 'Striatal.Subregion'], levels = order_genotype)
stats = lapply(stats, convert_factors, cols = c('Striatal.Subregion'), levels = order_subregion)

# To test for Genotype differences in All, D1, D2
nSets = 1 

total_num_trials = nNumStats 
test_list = list()
for (set in 1:nSets){
  # print(set)
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
      
      # anova test
      formula = as.formula(paste(col, '~', cat))
      test_result =  aov(formula, data = stats[[set]])
      test_result = summary(test_result)
      statistic = test_result[[1]]$`F value`[[1]][1] 
      p_value = test_result[[1]]$`Pr(>F)`[[1]][1]
      
      if(p_correction == TRUE){
        p_value_adjusted = p.adjust(p_value, method = "bonferroni", n = total_num_trials)
      }
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
  rownames(test_list[[set]]) = test_list[[set]]$ShapeParameter
  
  write.csv(test_list[[set]], spaste("Results/associationOfShapeParameters-", statNames[set], ".csv"), row.names = FALSE)
}

###############################################################################################
######### 5. Heatmap of Normalized Features -subregion pairwise comparison #######################
######### (All WT, All Q140, followed by each individual brain in Q and WT) ###################
###############################################################################################
# normalize features to its median by brain
normalized_df = stats[[1]] %>% group_by(Brain) %>%
  mutate(across(all_of(numericCols[[1]]), ~ . / median(., na.rm = TRUE)))

# set levels for brains
order_brains = c('hTME15-1', 'hTME19-2', 'hTME16-1', 'hTME24-2',
                 'hTME15-2', 'hTME18-1', 'hTME20-1')
normalized_df$Brain = factor(normalized_df$Brain, levels = order_brains)

# create brain_geno variable
normalized_df$brain_geno = paste(normalized_df$Brain, normalized_df$Genotype, sep = " (") %>%
  sapply(function(x) paste0(x, ")"))

# re-order brain_geno by genotype
unique_brain_geno = c(paste0(c('hTME15-1', 'hTME19-2', 'hTME16-1', 'hTME24-2'), rep(' (WT)',4)), 
                      paste0(c('hTME15-2', 'hTME18-1', 'hTME20-1'), rep(' (Q140)',3)))

# calculate mean in region for All WT, All Q140 and by individual brain
mean_df_brain = normalized_df %>% group_by(brain_geno, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup()

mean_df_all_Q_WT = normalized_df %>% group_by(Genotype, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup() %>% mutate(Genotype = if_else(Genotype == 'WT', 'All WT', 'All Q140')) %>%
  rename(brain_geno = Genotype)

mean_df = bind_rows(mean_df_all_Q_WT, mean_df_brain)

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
    pivot_longer(cols = -brain_geno, names_to = "Feature", values_to = "Value") %>%
    pivot_wider(names_from = brain_geno, values_from = Value) %>% as.data.frame()
  
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
stats[[1]]$brain_geno = paste(stats[[1]]$Brain, stats[[1]]$Genotype, sep = " (") %>%
  sapply(function(x) paste0(x, ")"))

# Function to add All WT and All Q140 groups
create_dup_brain_geno = function(df){
  df_dup = copy(df)
  df_dup$brain_geno = ifelse(df_dup$Genotype == 'WT', 'All WT', 'All Q140')
  df_comb = rbind(df, df_dup)
  return(df_comb)
}

stats_with_dup_brain_geno = create_dup_brain_geno(stats[[1]])
unique_brain_geno_all = c('All WT', 'All Q140', unique_brain_geno)
stats_with_dup_brain_geno$brain_geno = factor(stats_with_dup_brain_geno$brain_geno, unique_brain_geno_all)

# create subsets for region-pairs
CPr_CPi_sub = stats_with_dup_brain_geno %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPc_sub = stats_with_dup_brain_geno %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPi_CPc_sub = stats_with_dup_brain_geno %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))

# prettify brain names
replace_substrings = function(v) {
  # 12m WT
  v = gsub("hTME15-1\\s*\\(WT\\)", "12m_WT#1", v)
  v = gsub("hTME19-2\\s*\\(WT\\)", "12m_WT#2", v)
  v = gsub("hTME16-1\\s*\\(WT\\)", "12m_WT#3", v)
  v = gsub("hTME24-2\\s*\\(WT\\)", "12m_WT#4", v)
  
  # 12m Q140
  v = gsub("hTME15-2\\s*\\(Q140\\)", "12m_Q140#1", v)
  v = gsub("hTME18-1\\s*\\(Q140\\)", "12m_Q140#2", v)
  v = gsub("hTME20-1\\s*\\(Q140\\)", "12m_Q140#3", v)
  
  return(v)
}

# reorder heatmap rows by 3 main groups 
# define 3 groups of features and assign colors
features_angle = c('Bif_ampl_local', 'Bif_ampl_remote', 'Bif_tilt_local', 'Bif_tilt_remote', 'Bif_torque_local', 'Bif_torque_remote', 'Centripetal_Bias')
features_length = c('Depth', 'Height', 'Width', 'Length', 'Sum_EucDistance', 'Sum_PathDistance', 'Max_EucDistance', 'Max_PathDistance', 'ABEL_All', 'ABEL_Internal', 'ABEL_Terminal', 'BAPL_All', 'BAPL_Internal', 'BAPL_Terminal')
features_complexity = c('N_bifs', 'N_branch', 'N_tips', 'N_stems', 'Branch_Order', 'Fractal_Dim', 'Partition_asymmetry', 'Terminal_degree', 'Balancing_Factor', 'Convexity')

feature_reorder = data.frame(Feature = c(features_angle, features_length, features_complexity),
                             group = rep(1:3, times = c(length(features_angle), length(features_length), length(features_complexity))))


draw_heatmap = function(data_for_heatmap, data_subset, group_name) {
  
  # generate p-value matrix for pairwise subregion comparison per brain
  p_values = list()
  statistic_values = list()
  
  for(geno in unique_brain_geno_all) {
    sub_data = data_subset[data_subset$brain_geno == geno, ]
    unique_region = unique(as.character(data_subset$Striatal.Subregion))
    
    # Loop through each feature and perform a test
    p_values[[geno]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(sub_data[sub_data$Striatal.Subregion == unique_region[1], feature],
                     sub_data[sub_data$Striatal.Subregion == unique_region[2], feature])
      test$p.value
    })
    
    statistic_values[[geno]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(sub_data[sub_data$Striatal.Subregion == unique_region[1], feature],
                     sub_data[sub_data$Striatal.Subregion == unique_region[2], feature])
      test$statistic
    })
  }
  
  df_pval = as.data.frame(p_values, col.names = paste0('p.value_', unique_brain_geno_all))
  df_statistic = as.data.frame(statistic_values, col.names = paste0('statistics_', unique_brain_geno_all))
  
  adjusted_p_values = lapply(p_values, function(p) p.adjust(p, method = "bonferroni", n = nNumStats))
  adjusted_p_values_df = as.data.frame(adjusted_p_values, col.names = paste0('p.value_adjusted_', unique_brain_geno_all))
  
  rownames(adjusted_p_values_df) = rownames(df_pval) = rownames(df_statistic) = numericCols[[1]]
  
  # write data to a csv file
  df_output = cbind(df_statistic, df_pval, adjusted_p_values_df)
  ordered_columns = as.vector(outer(
    c('statistics_', 'p.value_', 'p.value_adjusted_'),
    gsub('\\)', '.', gsub(' \\(', '..', gsub('-', '.', unique_brain_geno_all))),
    paste0
  ))
  ordered_columns = gsub(' ', '.', ordered_columns)
  
  df_output = df_output[, ordered_columns]
  write.csv(df_output, paste0('Results/', group_name, '_statistic_output_by_individual_brain.csv'), row.names = TRUE)
  
  colnames(adjusted_p_values_df) = unique_brain_geno_all
  
  significance_df = apply(adjusted_p_values_df, 2, function(x) {
    ifelse(x < 0.001, "***",
           ifelse(x < 0.01, "**",
                  ifelse(x < 0.05, "*", "")))
  })
  
  # sort data matrix by re-ordered features
  data_for_heatmap = data_for_heatmap[unique_brain_geno_all]
  data_for_heatmap = round(data_for_heatmap, 2)
  data_for_heatmap = data_for_heatmap[match(feature_reorder$Feature, rownames(data_for_heatmap)), ]
  # sort p_significance matrix by re-ordered features
  significance_df = significance_df[match(feature_reorder$Feature, rownames(significance_df)), ]
  
  # add annotation
  annotation_df = data.frame(Cluster = feature_reorder[,'group'])
  rownames(annotation_df) = feature_reorder$Feature
  
  palette_colors = c("darkorange", "lightseagreen", "mediumorchid")
  annotation_color_mapping = setNames(palette_colors, c('1', '2', '3'))
  ### 
  
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
  
  # use prettified names for brains
  colnames(data_for_heatmap) = colnames(significance_df) = replace_substrings(colnames(data_for_heatmap))
  
  # generate heatmap
  ht = Heatmap(as.matrix(data_for_heatmap), 
               name = " ", 
               col = blue_white_red,
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE, row_title = 'Morphometrics',
               show_column_names = TRUE, column_title = 'Brain',
               width = unit(5, 'inch'), 
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
pdf(file = "Plots/5_heatmap - normalized_feature_region_pairwise_comparison_by_brain_with_all_pretty-name_new_feature_grouping.pdf", wi = 9, he = 8)
scpp(2.5);
par(mar = c(3,5,3,5));
  draw_heatmap(CPr_CPi_heatmap_data, CPr_CPi_sub, "CPr vs. CPi")
  draw_heatmap(CPr_CPc_heatmap_data, CPr_CPc_sub, "CPr vs. CPc")
  draw_heatmap(CPi_CPc_heatmap_data, CPi_CPc_sub, "CPi vs. CPc")
dev.off()

###############################################################################################
######### 6. Heatmap of Normalized Features, subregion pairwise comparison ######################
######### The above heatmap further divided by D1/D2 ##########################################
###############################################################################################
# calculate mean in region by Genotype, and with All WT and all Q140 group
mean_df_d1d2 = normalized_df %>% group_by(brain_geno, Type, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup()

mean_df_d1d2_with_all = normalized_df %>% group_by(Genotype, Type, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup() %>% mutate(Genotype = case_when(
    Genotype == 'WT'  ~ 'All WT',
    Genotype == 'Q140'  ~ 'All Q140'
  )) %>% rename(brain_geno = Genotype)

mean_df_d1d2 = bind_rows(mean_df_d1d2, mean_df_d1d2_with_all)

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
    pivot_longer(cols = -c(brain_geno,Type), names_to = "Feature", values_to = "Value") %>%
    pivot_wider(names_from = c(brain_geno,Type), values_from = Value) %>% as.data.frame()
  
  rownames(data) = data$Feature
  data = data[, -1]
}

CPr_CPi_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPr_CPi_diff_d1d2)
CPr_CPc_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPr_CPc_diff_d1d2)
CPi_CPc_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPi_CPc_diff_d1d2)

# prepare subset data for pairwise subregion comparison within each brain
CPr_CPi_sub = stats_with_dup_brain_geno %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPi_sub$geno_type = paste0(CPr_CPi_sub$brain_geno, '_', CPr_CPi_sub$Type)
CPr_CPc_sub = stats_with_dup_brain_geno %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPr_CPc_sub$geno_type = paste0(CPr_CPc_sub$brain_geno, '_', CPr_CPc_sub$Type)
CPi_CPc_sub = stats_with_dup_brain_geno %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))
CPi_CPc_sub$geno_type = paste0(CPi_CPc_sub$brain_geno, '_', CPi_CPc_sub$Type)

# define the desired order of the heatmap columns
all_combs = expand.grid(unique_brain_geno_all, c('_D1', '_D2'))
all_combs$Var1 = factor(all_combs$Var1, levels = c('All WT', 'All Q140', unique_brain_geno))
all_combs = all_combs %>% arrange(Var1, Var2)
order_geno_type = paste0(all_combs$Var1, all_combs$Var2)

# Function to draw heatmap
draw_heatmap_d1d2 = function(data_for_heatmap, data_subset, group_name) {
  
  data_for_heatmap = data_for_heatmap[, order_geno_type]
  
  # generate p-value matrix for pairwise subregion comparison per brain
  p_values = list()
  statistic_values = list()
  
  for(geno_type in order_geno_type) {
    sub_data = data_subset[data_subset$geno_type == geno_type, ]
    unique_region = unique(as.character(data_subset$Striatal.Subregion))
    
    # Loop through each feature and perform a t-test
    p_values[[geno_type]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(sub_data[sub_data$Striatal.Subregion == unique_region[1], feature],
                     sub_data[sub_data$Striatal.Subregion == unique_region[2], feature])
      test$p.value
    })
    
    statistic_values[[geno_type]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(sub_data[sub_data$Striatal.Subregion == unique_region[1], feature],
                     sub_data[sub_data$Striatal.Subregion == unique_region[2], feature])
      test$statistic
    })
  }
  
  df_pval = as.data.frame(p_values, col.names = paste0('p.value_', order_geno_type))
  df_statistic = as.data.frame(statistic_values, col.names = paste0('statistics_', order_geno_type))
  
  adjusted_p_values = lapply(p_values, function(p) p.adjust(p, method = "bonferroni", n = nNumStats))
  adjusted_p_values_df = as.data.frame(adjusted_p_values, col.names = paste0('p.value_adjusted_', order_geno_type))
  
  rownames(adjusted_p_values_df) = rownames(df_pval) = rownames(df_statistic) = numericCols[[1]]
  
  # write data to a csv file
  df_output = cbind(df_statistic, df_pval, adjusted_p_values_df)
  ordered_columns = as.vector(outer(
    c('statistics_', 'p.value_', 'p.value_adjusted_'),
    gsub('\\)', '.', gsub(' \\(', '..', gsub('-', '.', order_geno_type))),
    paste0
  ))
  ordered_columns = gsub(' ', '.', ordered_columns)
  df_output = df_output[, ordered_columns]
  write.csv(df_output, paste0('Results/', group_name, '_statistic_output_by_individual_brain_by_D1D2.csv'), row.names = TRUE)
  
  colnames(adjusted_p_values_df) = order_geno_type
  significance_df = apply(adjusted_p_values_df, 2, function(x) {
    ifelse(x < 0.001, "***",
           ifelse(x < 0.01, "**",
                  ifelse(x < 0.05, "*", "")))
  })
  
  # sort data matrix by re-ordered features
  data_for_heatmap = data_for_heatmap[, order_geno_type]
  data_for_heatmap = round(data_for_heatmap, 2)
  data_for_heatmap = data_for_heatmap[match(feature_reorder$Feature, rownames(data_for_heatmap)), ]
  # sort p_significance matrix by re-ordered features
  significance_df = significance_df[match(feature_reorder$Feature, rownames(significance_df)), ]
  
  # add annotation
  annotation_df = data.frame(Cluster = feature_reorder[,'group'])
  rownames(annotation_df) = feature_reorder$Feature
  
  palette_colors = c("darkorange", "lightseagreen", "mediumorchid")
  annotation_color_mapping = setNames(palette_colors, c('1', '2', '3'))
  ### 
  
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
  
  # prettify names of brains
  colnames(data_for_heatmap) = colnames(significance_df) = replace_substrings(colnames(data_for_heatmap))
  
  # generate heatmap
  ht = Heatmap(as.matrix(data_for_heatmap), 
               name = " ", 
               col = blue_white_red,
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE, row_title = 'Morphometrics',
               show_column_names = TRUE, column_title = 'Genotype_Type',
               width = unit(7, 'inch'), 
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
pdf(file = "Plots/6_heatmap - normalized_feature_region_pairwise_comparison_by_brain_byD1D2_with_all_pretty-name_new_feature_grouping.pdf", wi = 11, he = 9)
scpp(2.5);
par(mar = c(3,5,3,5));
  draw_heatmap_d1d2(CPr_CPi_heatmap_data_d1d2, CPr_CPi_sub, "CPr vs. CPi")
  draw_heatmap_d1d2(CPr_CPc_heatmap_data_d1d2, CPr_CPc_sub, "CPr vs. CPc")
  draw_heatmap_d1d2(CPi_CPc_heatmap_data_d1d2, CPi_CPc_sub, "CPi vs. CPc")
dev.off()

###############################################################################################
################## 7. Blue & Red Heatmap - Q140 vs. WT on Morphometrics ########################
###############################################################################################
# function to normalize the colors row-wise (feature-wise)
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

# function to generate blue/red color heatmap for Q140 compared to WT (baseline)
colored_heatmap_by_genotype = function(data, cat_features, set_name, cat_name){
  data_matrix = matrix(NA, nrow = nNumStats)
  text_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  coef_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  
  for(cat in cat_features){
    results = data.frame(Feature = colnames(numStats[[1]]),
                         Coefficient = NA, 
                         P_Value = NA)
    
    for (feature in results$Feature) {
      
      sub_data = data[!is.na(data[, cat]), ]
      sub_data[, cat] = factor(sub_data[, cat], levels = c('WT', 'Q140'))
      column_means = colMeans(sub_data[sub_data$Genotype == 'WT', numericCols[[1]]], na.rm = TRUE)
      results$Baseline_mean = column_means
      
      # use t-test
      formula = as.formula(paste(feature, " ~ ", 'Genotype'))
      t_test_result = t.test(formula, data = sub_data, var.equal = var_equal_global)
      p_value = t_test_result$p.value
      coef = mean(sub_data[sub_data$Genotype == 'Q140', feature]) - mean(sub_data[sub_data$Genotype == 'WT', feature])
      
      if(p_correction){
        p_value = p.adjust(p_value, method = 'bonferroni', n = nNumStats)
      }
      
      results[results$Feature == feature, "Coefficient"] = coef
      results[results$Feature == feature, "P_Value"] = p_value
      results$Fold_change = round(abs(results$Coefficient) / results$Baseline_mean,2)
    }
    
    results$P_Value = signif(results$P_Value, 2)
    # account for P-value = 1
    results[results$P_Value == 1, 'P_Value'] = 0.95 
    
    results$neg_log10p = -log10(results$P_Value)
    results$Coefficient = signif(results$Coefficient,3)
    
    data_matrix = cbind(data_matrix, results$neg_log10p)
    # conditionally round coefficients, for Fractal_Dim, Terminal_degree, Partition_asymmetry, Balancing_Factor, Convexity: 3 decimal places
    # for other features: 1 decimal places 
    small_features = c('Fractal_Dim', 'Terminal_degree', 'Partition_asymmetry', 'Balancing_Factor', 'Convexity')
    results_coef = results$Coefficient
    
    results_coef[!(results$Feature %in% small_features)] = round(results_coef[!(results$Feature %in% small_features)], 2)
    results_coef[results$Feature %in% small_features] = sprintf("%.1e", results_coef[results$Feature %in% small_features])
    
    # color_matrix = cbind(color_matrix, results$Color)
    text_matrix = text_matrix = cbind(text_matrix, results_coef)
    coef_matrix = cbind(coef_matrix, results$Coefficient)
  }
  
  data_matrix = data_matrix[,-1]
  # color_matrix = color_matrix[,-1]
  text_matrix = text_matrix[,-1]
  coef_matrix = coef_matrix[,-1]
  
  # add * based on -log10P, if p<0.05 then -log10P should > 1.301
  text_matrix = matrix(as.character(text_matrix), nrow=nrow(text_matrix))
  text_matrix[data_matrix > 1.301] = paste0(text_matrix[data_matrix > 1.301], '*')
  text_matrix[data_matrix > 2.301] = paste0(text_matrix[data_matrix > 2.301], '*')
  text_matrix[data_matrix > 3.301] = paste0(text_matrix[data_matrix > 3.301], '*')
  
  # set row and column names of the matrix
  colLabs = gsub("Genotype", "Q140 vs. WT", cat_features)
  rowLabs = results$Feature
  rownames(data_matrix) = rownames(text_matrix) = rownames(coef_matrix) = rowLabs
  colnames(data_matrix) = colnames(text_matrix) = colnames(coef_matrix) = colLabs
  
  # save text_matrix
  text_only_asterisk = apply(text_matrix, c(1, 2), function(x) gsub("[^*]", "", x))
  write.csv(text_only_asterisk,  file = paste0('Q-vs-WT-', cat_name, '.csv'), row.names = TRUE)
  
  # adjust shade based on significance
  color_matrix = rowwise_color_gradient(coef_matrix)
  color_matrix[data_matrix<1.301 & coef_matrix > 0] = "#FFF0F0"
  color_matrix[data_matrix<1.301 & coef_matrix < 0] = "#F0F0FF"
  
  # reorder and group features
  data_matrix = data_matrix[feature_reorder$Feature, ]
  color_matrix = color_matrix[feature_reorder$Feature, ]
  text_matrix = text_matrix[feature_reorder$Feature, ]
  coef_matrix = coef_matrix[feature_reorder$Feature, ]
    
  # to map cell colors using color_matrix for the heatmap
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
                 grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 8))
               },
               right_annotation = right_annotation,
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
               show_heatmap_legend = FALSE
  )
  
  draw(ht)
  grid.text(paste0("Mean Difference and Significance Heatmap (Q140 vs. WT), ", set_name),
            x = unit(0.5, "npc"),
            y = unit(0.98, "npc"),  
            gp = gpar(fontsize = 12, fontface = 'bold'))
  
  # Define the viewport to place the legend 
  viewport = viewport(x = 0.5, y = 0.94, width = 1, height = 0.5)
  pushViewport(viewport)
  
  # Add the custom legend
  grid.rect(gp = gpar(fill = "red"), x = 0.2, width = 0.06, height = 0.03)
  grid.text("Q140 higher than WT", x = 0.24, just = "left", gp = gpar(fontsize = 8))
  grid.rect(gp = gpar(fill = "blue"), x = 0.5, width = 0.06, height = 0.03)
  grid.text("Q140 lower than WT", x = 0.54, just = "left", gp = gpar(fontsize = 8))
  
  # Reset viewport
  popViewport()
}

pdf(file = spaste("Plots/7_Heatmap-Association-QvsWT-Blue-Red-",'All_ttest', ".pdf"), wi = 9, he = 9);
scpp(2.5);
par(mar = c(3,5,3,5));
  # 1). All, by region
  colored_heatmap_by_genotype(stats[[1]], c('Genotype', 'Genotype_in_CPr', 'Genotype_in_CPi', 'Genotype_in_CPc'), 'All', 'Genotype_in_Region') 
  # 2) All, by D1D2
  colored_heatmap_by_genotype(stats[[1]], c('Genotype', 'Genotype_in_D1', 'Genotype_in_D2'), 'All', 'Genotype_in_Type') 
  # 3) All, by Type in Regions
  colored_heatmap_by_genotype(stats[[1]], c('Genotype', 
                                            'Genotype_in_D1_in_CPr', 'Genotype_in_D2_in_CPr',
                                            'Genotype_in_D1_in_CPi', 'Genotype_in_D2_in_CPi',
                                            'Genotype_in_D1_in_CPc', 'Genotype_in_D2_in_CPc'), 'All', 'Genotype_in_Type_in_Region') 
dev.off()

###############################################################################################
############ 8. D1 vs. D2 differences, overall and by region, in Q compared to WT ################
###############################################################################################
# create variables
stats_type_geno = lapply(stats, function(stats) 
{
  # type in Genotype
  rest1 = restrictVariableByCovariateLevels(stats$Type, stats$Genotype, 
                                            varName = "Type", covarName = "", nameSep = " in ", check.names = FALSE);
  stats = data.frame.ncn(rest1,  stats);
  out = setNames(stats, gsub(" ", "_", names(stats)))
});

stats_type_geno = lapply(stats_type_geno, function(stats) 
{
  # type in WT in striatal subregions
  rest1 = restrictVariableByCovariateLevels(stats$Type_in_WT, stats$Striatal.Subregion, 
                                            varName = "Type_in_WT", covarName = "", nameSep = " in ", check.names = FALSE);
  # type in Q140 in striatal subregions
  rest2 = restrictVariableByCovariateLevels(stats$Type_in_Q140, stats$Striatal.Subregion, 
                                            varName = "Type_in_Q140", covarName = "", nameSep = " in ", check.names = FALSE);
  stats = data.frame.ncn(rest1, rest2, stats);
  out = setNames(stats, gsub(" ", "_", names(stats)))
});

# function to generate the heatmap showing D2 vs. D1 using D1 as the baseline
colored_heatmap_by_type_by_genotype = function(data, cat_features, set_name){
  data_matrix = matrix(NA, nrow = nNumStats)
  text_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  coef_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  
  for(cat in cat_features){
    print(cat)
    results = data.frame(Feature = colnames(numStats[[1]]),
                         Coefficient = NA, 
                         P_Value = NA)
    
    for (feature in results$Feature) {
      
      sub_data = data[!is.na(data[, cat]), ]
      sub_data[, cat] = factor(sub_data[, cat], levels = c('D1', 'D2'))
      column_means = colMeans(sub_data[sub_data$Type == 'D1', numericCols[[1]]], na.rm = TRUE)
      results$Baseline_mean = column_means
      
      # t-test
      formula = as.formula(paste(feature, " ~ ", cat))
      test_result = t.test(formula, sub_data)
      p_value = test_result$p.value
      coef = mean(sub_data[sub_data$Type == 'D2', feature]) - mean(sub_data[sub_data$Type == 'D1', feature])
      
      if(p_correction){
        p_value = p.adjust(p_value, method = 'bonferroni', n = nNumStats)
      }
      
      results[results$Feature == feature, "Coefficient"] = coef
      results[results$Feature == feature, "P_Value"] = p_value
      results$Fold_change = round(abs(results$Coefficient) / results$Baseline_mean,1)
    }
    
    results$P_Value = signif(results$P_Value, 2)
    results[results$P_Value == 1, 'P_Value'] = 0.95 
    
    results$neg_log10p = -log10(results$P_Value)
    results$Coefficient = signif(results$Coefficient,3)
    
    data_matrix = cbind(data_matrix, results$neg_log10p)
    text_matrix = cbind(text_matrix, round(results$Coefficient,2))
    coef_matrix = cbind(coef_matrix, results$Coefficient)
  }
  
  data_matrix = data_matrix[,-1]
  text_matrix = text_matrix[,-1]
  coef_matrix = coef_matrix[,-1]
  
  # add * based on -log10P, if p<0.05 then -log10P should > 1.301
  text_matrix = matrix(as.character(text_matrix), nrow=nrow(text_matrix))
  text_matrix[data_matrix > 1.301] = paste0(text_matrix[data_matrix > 1.301], '*')
  text_matrix[data_matrix > 2.301] = paste0(text_matrix[data_matrix > 2.301], '*')
  text_matrix[data_matrix > 3.301] = paste0(text_matrix[data_matrix > 3.301], '*')
  
  # set row and column names of the matrix
  colLabs = gsub("Type", "D2 vs. D1", cat_features)
  rowLabs = results$Feature
  rownames(data_matrix) = rownames(text_matrix) = rownames(coef_matrix) = rowLabs
  colnames(data_matrix) = colnames(text_matrix) = colnames(coef_matrix) = colLabs
  
  # ade based on difference but not p-value
  color_matrix = rowwise_color_gradient(coef_matrix)
  color_matrix[data_matrix<1.301 & coef_matrix > 0] = "#FFF0F0"
  color_matrix[data_matrix<1.301 & coef_matrix < 0] = "#F0F0FF"
      
  # to map cell colors using color_matrix
  custom_color_fun = function(value, matrix = coef_matrix, colors = color_matrix) {
    color_index = match(value, matrix)
    return(colors[color_index])
  }
  
  # group and reorder features
  data_matrix = data_matrix[feature_reorder$Feature, ]
  color_matrix = color_matrix[feature_reorder$Feature, ]
  text_matrix = text_matrix[feature_reorder$Feature, ]
  coef_matrix = coef_matrix[feature_reorder$Feature, ]
  
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
  grid.text(paste0("Coefficient and Significance Heatmap (D2 vs. D1), ", set_name),
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

cat_features = c('Type_in_WT', 'Type_in_Q140', 
                 'Type_in_WT_in_CPr', 'Type_in_Q140_in_CPr',
                 'Type_in_WT_in_CPi', 'Type_in_Q140_in_CPi',
                 'Type_in_WT_in_CPc', 'Type_in_Q140_in_CPc')

pdf(file = spaste("Plots/8_Heatmap-Association-Type-Blue-Red-D1vsD2_by_Genotype",'', ".pdf"), 
    wi = 10, he = 9);
scpp(2.5);
par(mar = c(3,5,3,5));
  # 1). only D2vsD1 in Q and WT
  colored_heatmap_by_type_by_genotype(stats_type_geno[[1]], c('Type_in_WT', 'Type_in_Q140'), 'All ') 
  # 2). D2vsD1 in Q and WT by region
  colored_heatmap_by_type_by_genotype(stats_type_geno[[1]], cat_features, 'All ') 
dev.off()

###############################################################################################
##########  9. Boxplot, comparison between striatal subregions ##############################
##########  WT, WT D1, WT D2, Q140, Q140 D1, Q140 D2  #########################################
###############################################################################################

# Function to calculate overall significance of P-values of different morphometric features vs. Striatal.Subregion
get_overall_pvalue = function(data, cat, variate) {
  formula = as.formula(paste(cat, " ~ ", variate))
  test = aov(formula, data)
  test_summary = summary(test)
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
  
  # for each feature
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

# box plot distribution showing CPr, CPi and CPc
# in Q140, Q140 D1, Q140 D2, WT, WT D1, WT D2
subset_list = list(stats[[4]], subset(stats[[4]], Type == 'D1'), subset(stats[[4]], Type == 'D2'),
                   stats[[5]], subset(stats[[5]], Type == 'D1'), subset(stats[[5]], Type == 'D2'))
subset_name_list = c('WT', 'WT,D1', 'WT,D2', 'Q140', 'Q140,D1', 'Q140,D2')
set_name = 'All'

pdf(spaste("Plots/9_Boxplots-by-level.", set_name, ".pdf"), wi = 9, he = 9)
scpp(2.5);
par(mar = c(6.3, 3, 2.5, 1));
  level_comparison(subset_list, subset_name_list, set_name)
dev.off()

###############################################################################################
################## 10. Box plot:  Q vs. WT Group Comparison (overall, by level) ###################
###############################################################################################
comparison_var1_across_var2 = function(data, compare_var, across_var, set_name){
  all_data = data
  all_data[, across_var] = 'All'
  combined_data = rbind(data, all_data)
  combined_data$Genotype = factor(combined_data$Genotype, levels = order_genotype)
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
    p_values = sapply(split(combined_data, combined_data[, across_var]), function(x) get_overall_pvalue(x, cat, 'Genotype'))
    
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
        scale_fill_manual(values = c("D1" = "darkgreen", "D2" = "darkred")) 
      }else if(compare_var == 'Genotype'){
        scale_fill_manual(values = c("WT" = "darkblue", "Q140" = "red")) 
      }else{
        scale_fill_brewer(palette = "Dark2") 
      } 
    # Add count annotations
    p = p + geom_text(data = counts_data, aes(x = PosX, y = PosY, label = paste("n =", Count)), 
                      position = position_dodge(width = 0.75), 
                      vjust = 1.2, size = 3, color = "black")
    # Add the annotations to the plot (overall p-values for ANOVA)
    p = p + geom_text(data = annotation_data, aes(label = Label, x = x, y = y),
                      hjust = 0.5, vjust = 0,  
                      inherit.aes = FALSE, size = 3) 
    print(p)
  }
}

# Q/WT comparison, across levels
pdf("Plots/10_Boxplots-Genotype.by.level.All.pdf", width = 7, height = 5)
scpp(2);
  comparison_var1_across_var2(stats[[1]], 'Genotype', 'Striatal.Subregion', 'All')
dev.off()

##########################################################################################################
#################### 11. Box plot: feature distribution by brain grouped by genotype #########################
##########################################################################################################
pdf(paste0("Plots/11_Boxplot", ".by.Brain", ".pdf"), width = 5, height = 4)
scpp(2);
for(cat in numericCols[[1]]){
  p = ggplot(data = stats[[1]], aes(x = Genotype, y = !!sym(cat), fill = Brain)) +
    geom_boxplot(position = position_dodge(width = 0.8),
                 outlier.shape = NA, 
                 size = 0.5) +
    theme_minimal() +
    theme(
      legend.position = "bottom",    
      legend.title = element_blank(), 
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text( face = 'bold', hjust = 1),         
      axis.title = element_text(size = rel(1.2)),        
      strip.text = element_text(size = rel(1.2)),        
      panel.border = element_rect(color = "black", fill = NA, size = 1), 
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(), 
      plot.background = element_rect(fill = "white", color = NA) 
    ) 
  print(p)
}
dev.off()

###############################################################################################
############################ 12. Distribution of Brain in CPr, CPi and CPc ########################
###############################################################################################
# define a function to generate stacked bar plot of level distribution within each cluster
stacked_bar_distribution = function(data, by_var, stack_var, set_name){
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
  
  combined_data = bind_rows(
    overall_proportions,
    cluster_proportions
  ) %>%
    arrange(!!sym(by_var))
  
  if(by_var == 'Striatal.Subregion'){
    combined_data$Striatal.Subregion = factor(combined_data$Striatal.Subregion, levels = c('All', order_subregion))
  }
  
  # group counts
  grp_cnts = combined_data %>% group_by(!!sym(by_var)) %>% summarize(n = sum(n))
  grp_cnt_strings = paste(grp_cnts[[by_var]], " (", grp_cnts$n, ")", sep = "")
  
  # chi-square
  # Chi-square test
  pval_chi = c("") # placeholder for 'All' group
  # get the expected proportion of by_var
  expected_prop = prop.table(table(data[, stack_var]))
  # chi-square goodness of fit test for each cluster
  for(grp in names(table(data[, by_var]))){
    ct = chisq.test(table(data[data[, by_var] == grp, stack_var]), p = expected_prop)
    print(paste0(grp, ct$p.value))
    pval_chi = append(pval_chi, ifelse(ct$p.value <= 0.05, '*', ''))
  }
  combined_data$Asterisks = rep(pval_chi, each=length(table(data[, stack_var])))
  
  # Create the stacked bar plot, display percentage on it
  ggplot(combined_data, aes(x = get(by_var), y = Proportion, fill = get(stack_var), label = Asterisks)) +
    geom_bar(stat = "identity", width = 0.5) +
    # scale_fill_manual(values = wes_palette("Darjeeling1")) +
    scale_fill_manual(values = as.character(paletteer_d("ggthemes::Classic_Blue_Red_12"))) + 
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
              size = 4, 
              color = "black") +
    geom_text(aes(y = 1, label = Asterisks), vjust = 0, color = "red") +    
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() + 
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = 14),  # Adjust x-axis label size
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),  # Adjust x-axis tick label size
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
  
}

pdf(paste0("Plots/12_Distribution_Plots_brain_level.pdf"), wi = 7, he = 7)
  stacked_bar_distribution(stats[[1]], 'Striatal.Subregion', 'Brain', 'All')
dev.off()

###############################################################################################
################ 13. Radar plot: pct% change of mean Q140 vs. WT, in CPr/i/c ######################
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

# function to calculate percentage
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

mean_df = stats_dup %>% group_by(Type, Genotype, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup() %>% as.data.frame()

pct_changes = mean_df %>% group_by(Type, Striatal.Subregion) %>%
  summarise(across(numericCols[[1]], ~ percent_change(.x[2], .x[1])))

print(min(pct_changes[, numericCols[[1]]]))
plot_min = -20
print(max(pct_changes[, numericCols[[1]]]))
plot_max = 30

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
  radarchart(df, axistype = 4, seg= 5,
             # Custom polygon
             pcol = colors_border,plwd = 3, plty = 1,
             pfcol = NULL,
             
             # Custom the grid
             cglcol = "grey", cglty = 1, axislabcol = "black", 
             caxislabels = c('-20%\n', '-10%\n', '0\n', '10%\n', '20%\n', '30%\n'), 
             cglwd = 0.5,
             
             # Custom labels
             vlcex = 0.8)
  
  n_vertices = length(numericCols[[1]]) * 2
  angles = seq(0, 2 * pi, length.out = n_vertices)
  radius_0 = 0.5
  x_coords = radius_0 * cos(angles)
  y_coords = radius_0 * sin(angles)
  lines(x_coords, y_coords, col = "black", lwd = 4, lty = 2)
  
  title(main = paste0("Percent (%) change of means - Q140 vs. WT\n", grp), col.main = "black", font.main = 4)
  legend(x = 1.2, y = 1, legend = c('All', 'D1', 'D2'), bty = "n", pch = 20, col = colors_border, 
         text.col = "black", cex = 1.2, pt.cex = 3)
  
}

pdf(file = spaste("Plots/13_Radar_CPr.CPi.CPc_pct_change_byD1D2.pdf"), wi = 13, he = 10);
scpp(2.5);
par(mar = c(3,5,3,5));
  create_radar_pct_change(subset(pct_changes, Striatal.Subregion == 'CPr'), 'CPr', plot_min, plot_max)
  create_radar_pct_change(subset(pct_changes, Striatal.Subregion == 'CPi'), 'CPi', plot_min, plot_max)
  create_radar_pct_change(subset(pct_changes, Striatal.Subregion == 'CPc'), 'CPc', plot_min, plot_max)
dev.off()

###############################################################################################
################ 14. Same Radar plot, but for all regions combined  ########################
###############################################################################################
mean_df = stats_dup %>% group_by(Type, Genotype) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup() %>% as.data.frame()

pct_changes = mean_df %>% group_by(Type) %>%
  summarise(across(numericCols[[1]], ~ percent_change(.x[2], .x[1])))

print(min(pct_changes[, numericCols[[1]]]))
plot_min = -20
print(max(pct_changes[, numericCols[[1]]]))
plot_max = 30

pdf(file = spaste("Plots/14_Radar_pct_change_byD1D2_All_regions.pdf"), wi = 13, he = 10);
scpp(2.5);
par(mar = c(3,5,3,5));
  create_radar_pct_change(pct_changes, 'All', plot_min, plot_max)
dev.off()

# ###############################################################################################
# ########### Write the significance to Excel to for plot annotation ###########################
# ###############################################################################################
# # get significance comparison, CPr, CPi and CPc separately
# CPr_sub = subset(stats_dup, Striatal.Subregion == 'CPr')
# CPi_sub = subset(stats_dup, Striatal.Subregion == 'CPi')
# CPc_sub = subset(stats_dup, Striatal.Subregion == 'CPc')
# 
# Q_WT_comparison_significance_df_byD1D2 = function(df, grp_name){
# 
#   p_values = list()
# 
#   for(t in c('All', 'D1', 'D2')){
#     df_sub = subset(df, Type == t)
#     p_values[[t]] = sapply(numericCols[[1]], function(feature) {
#       formula = as.formula(paste(feature, "~ Genotype"))
#       test = t.test(formula, data = df_sub, var.equal = TRUE)
#       test$p.value
#     })
#   }
# 
#   adjusted_p_values = lapply(p_values, function(p) p.adjust(p, method = "bonferroni", n = nNumStats))
#   adjusted_p_values_df = as.data.frame(adjusted_p_values)
#   rownames(adjusted_p_values_df) = numericCols[[1]]
#   colnames(adjusted_p_values_df) = c('All', 'D1', 'D2')
#   significance_df = apply(adjusted_p_values_df, 2, function(x) {
#     ifelse(x < 0.001, "***",
#            ifelse(x < 0.01, "**",
#                   ifelse(x < 0.05, "*", "")))
#   })
# 
#   return(significance_df)
# }
# 
# # Create a new workbook
# wb = createWorkbook()
#   # Add worksheets and write dataframes into them
#   addWorksheet(wb, 'CPr')
#   writeData(wb, sheet = 'CPr', Q_WT_comparison_significance_df_byD1D2(CPr_sub, 'CPr'), rowNames = TRUE)
#   addWorksheet(wb, 'CPi')
#   writeData(wb, sheet = 'CPi', Q_WT_comparison_significance_df_byD1D2(CPi_sub, 'CPi'), rowNames = TRUE)
#   addWorksheet(wb, 'CPc')
#   writeData(wb, sheet = 'CPc', Q_WT_comparison_significance_df_byD1D2(CPc_sub, 'CPc'), rowNames = TRUE)
# # Save the workbook
# saveWorkbook(wb, file = "Radar_CPr.CPi.CPc_pct_change_significance_df_byD1D2.xlsx", overwrite = TRUE)

