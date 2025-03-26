# set work directory
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
library(openxlsx)
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
library(tidyr)
library(dplyr)
library(stringr)
library(effsize)
library(nlme)
library(lme4)
library(emmeans)
library(circlize)
library(ComplexHeatmap)
library(tibble)

# create directions
dir.create("Plots", recursive = TRUE);
dir.create("Results", recursive = TRUE);

###############################################################################################
############################ read data #######################################
###############################################################################################
data_dir = '../morpho_1871CPneurons_excluded30.csv'
stats1 = read.csv(data_dir,check.names = FALSE)
colnames(stats1)
dim(stats1) # 1871

table(stats1$Type) # 1017 D1, 854 D2
table(stats1$Striatal.Subregion) # 279 CPr, 951 CPi, 641 CPc
table(stats1$Striatal.Community)

###############################################################################################
############################ Manipulation #####################################################
###############################################################################################
stats = list(stats1); 
statNames = c("D1D2")
names(stats) = statNames;

stats = lapply(stats, function(stats) 
{ 
  # Type in subregion (CPr, CPi, CPc)
  rest1 = restrictVariableByCovariateLevels(stats$Type, stats$Striatal.Subregion, 
                                            varName = "Type", covarName = "", nameSep = " in ", check.names = FALSE);
  # Type in sex
  rest2 = restrictVariableByCovariateLevels(stats$Type, stats$Sex, 
                                            varName = "Type", covarName = "", nameSep = " in ", check.names = FALSE);
  # Sex in subregion
  rest3 = restrictVariableByCovariateLevels(stats$Sex, stats$`Striatal.Subregion`, 
                                            varName = "Sex", covarName = "", nameSep = " in ", check.names = FALSE)
  # Sex in Type
  rest4 = restrictVariableByCovariateLevels(stats$Sex, stats$`Type`, 
                                            varName = "Sex", covarName = "", nameSep = " in ", check.names = FALSE)
  # Type in communities
  rest5 = restrictVariableByCovariateLevels(stats$Type, stats$Striatal.Community, 
                                            varName = "Type", covarName = "", nameSep = " in ", check.names = FALSE);
  
  stats = data.frame.ncn(rest1,rest2, rest3,rest4, rest5, stats);
  stats = setNames(stats, gsub(" ", "_", names(stats)))
  return(stats)
});

# datasets with numeric features only
firstStat = "Bif_ampl_local";
lastStat = "Convexity"
morf0 = lapply(stats, function(stats) setRownames(stats[match(firstStat, names(stats)):match(lastStat, names(stats))], 
                                                  make.unique(spaste(stats$`Brain`, ".", stats$`Reconstruction #`))));
numStats = lapply(morf0, dropConstantColumns);

nSets = length(stats);
numericCols = lapply(numStats, colnames);
nNumStats = sapply(numStats, ncol) 
print(nNumStats) # 31

# define levels of type, subregion and communities
order_type = c('D1', 'D2')
order_subregion = c('CPr','CPi','CPc')
order_community = c("CPr_m", 'CPr_imd', "CPr_imv", "CPr_l",
                    "CPi_dm", "CPi_dl", "CPi_vm", "CPi_vl",
                    "CPc_d", "CPc_i", "CPc_v")

# factorize type, subregion, community
factorize_columns = function(data, type_levels, subregion_levels, community_levels) {
  data$Type = factor(data$Type, levels = type_levels)
  data$Striatal.Subregion = factor(data$Striatal.Subregion, levels = subregion_levels)
  data$Striatal.Community = factor(data$Striatal.Community, levels = community_levels)
  return(data)
}

stats[[1]] = factorize_columns(stats[[1]], order_type, order_subregion, order_community)
stats1 = factorize_columns(stats1, order_type, order_subregion, order_community)

# equal_group_variance
var_equal_global = TRUE
# apply p-value correction
p_correction = TRUE 

# define continuous variable
continuous_var = colnames(numStats[[1]])

# number of independent trials for morphometric comparison
total_num_trials = nNumStats

###############################################################################################
############################ 1. Frequency Distribution  #################
###############################################################################################
# 1. Plot: frequency table
pdf('Plots/1_Frequency_Distribution.pdf')

  # Function to create and plot tables
  plot_table = function(data, var_names) {
    grid.table(data %>% rename(!!!var_names), rows = NULL)
    plot.new()
  }
  
  # Type frequency
  type_data = data.frame(table(stats1$Type)) 
  plot_table(type_data, c(Cell.Type = "Var1", Frequency = "Freq"))
  
  # Sex frequency
  sex_data = data.frame(table(stats1$Sex))
  plot_table(sex_data, c(Sex = "Var1", Frequency = "Freq"))
  
  # Brain frequency
  brain_data = data.frame(table(stats1$Brain))
  plot_table(brain_data, c(Brain = "Var1", Frequency = "Freq"))
  
  # Subregion frequency
  subregion_data = data.frame(table(stats1$Striatal.Subregion))
  plot_table(subregion_data, c(Sub.Region = "Var1", Frequency = "Freq"))
  
  # Type, sex, brain frequency
  type_sex_brain_data = data.frame(table(stats1$Brain, stats1$Type, stats1$Sex)) %>% 
    rename(Brain = Var1, Cell.Type = Var2, Sex = Var3, Frequency = Freq) %>% 
    arrange(Brain)
  plot_table(type_sex_brain_data, c(Brain = "Brain", Cell.Type = "Cell.Type", Sex = "Sex", Frequency = "Frequency"))

dev.off()

# 2. Plot: Distribution plot of D1 D2 in each subregion
# function: Add labels and theme adjustments
add_label_bar_plot = function(p, var1, var2) {
  p = p + 
    ylab('Count') +
    theme_classic() +
    theme(
      text = element_text(size = 15),
      axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
      plot.title = element_text(hjust = 0.5, vjust = 5),
      plot.margin = ggplot2::margin(1, 0.5, 0.5, 0.5, 'cm')  # t, r, b, l
    )
  
  xlab_title = switch(var1,
                      'Striatal.Subregion' = 'Striatal Subregions',
                      'Striatal.Communities' = 'Striatal Communities',
                      'Sex' = 'Sex',
                      var1)
  plot_title = paste0("Distribution of ", var2, " in ", xlab_title)
  p = p + xlab(xlab_title) + ggtitle(plot_title)
  
  return(p)
}

# function:stacked bar plot (with D1D2 proportion and total number displayed) in each subregion
distribution_plot = function(data, var1, var2) {
  
  # create dataset for annotate the bar plot
  a = table(data[, var1], data[, var2])
  b = prop.table(a, margin = 1)
  
  a = as.data.frame(a)
  colnames(a) = c(var1, var2, 'freq')
  
  b = as.data.frame(b)
  colnames(b) = c(var1, var2, 'prop')
  
  c = merge(a, b, by = c(var1, var2))
  
  # Determine totals for each group based on var1
  totals = c %>% group_by(.data[[var1]]) %>% summarize(total = sum(freq)) %>% as.data.frame()
  
  # Create base plot
  if(length(unique(data[, var1])) > 3){
    wi = 0.8
    font_size = 3
  }else{
    wi = 0.6
    font_size = 5
  }
  
  p = ggplot(c, aes(fill = eval(as.name(var2)), y = freq, x = .data[[var1]])) + 
    geom_bar(position = "stack", stat = "identity", width = wi) + 
    geom_text(size = font_size, position = position_stack(vjust = 0.5), label = paste0(round(c$prop * 100, 1), '%'), color = 'white') +
    geom_text(aes(.data[[var1]], total, label = total, fill = NULL), data = totals, vjust = -0.5, size = 5) +
    guides(fill = guide_legend(title = var2)) +
    scale_y_continuous(expand = c(0, 0))
  
  # Add custom scale_fill based on var2
  p = p + if (var2 == 'Type') {
    scale_fill_manual(values = c('darkgreen', 'darkred'))
  } else {
    scale_fill_brewer(palette = "Dark2")
  }
  
  # Adjust plot titles
  p = switch(var1,
             'Striatal.Subregion' = p + scale_y_continuous(limits = c(0, 1300)) + scale_x_discrete(),
             'Striatal.Community' = p + scale_x_discrete(expand = c(0.05, 0)) + 
               scale_y_continuous(limits = c(0, max(totals$total) + 50)),
             'Sex' = p + scale_y_continuous(limits = c(0, 1500)) + scale_x_discrete())
  
  p = add_label_bar_plot(p, var1, var2)
  print(p)
}

pdf(paste0("Plots/2_Distribution_Plots.pdf"), wi = 7, he = 7)
  distribution_plot(stats1, 'Striatal.Subregion', 'Type')
  distribution_plot(stats1, 'Striatal.Community', 'Type')
  distribution_plot(stats1, 'Sex','Type')
  distribution_plot(stats1, 'Striatal.Subregion','Sex')
dev.off()

# Plot: for paper, create singple plot for D1D2 distribution in community 
pdf(paste0("Plots/Type_distribution_in_community.pdf"), wi = 12, he = 7)
distribution_plot(stats1, 'Striatal.Community', 'Type')
dev.off()

###############################################################################################
############################ 4. histograms of each variable #######################################
###############################################################################################
# 4. Plot: feature histogram distribution
pdf(file = "Plots/4_statHistograms.pdf", wi = 5, he = 3);
scpp(2);
for (set in 1:nSets) for (col in 1:nNumStats[set])
  hist(numStats[[set]] [, col], breaks = 20, main = numericCols[[set]][col], 
       xlab = numericCols[[set]][col], ylab = "Frequency");

dev.off();

###############################################################################################
######## 5. correlation heatmap of variables  ########################
###############################################################################################
# 5. Plot: feature correlation heatmap
pdf(file = "Plots/5_correlationHeatmapOfStats.pdf", wi = 18, he = 13);
  for (set in 1:nSets)
    corHeatmapWithDendro(bicor(as.matrix(numStats[[set]]), use = 'p', maxPOutliers = 0.01),
                         mar.main = c(10, 13, 2, 1), main = "Correlations of shape statistics",
                         dendroWidth = 2/13.5)
dev.off()

###############################################################################################
#################### 6. hierarchical clustering of variables ########################
###############################################################################################
# Collapse the various stats by correlation
numStats.scaled = lapply(numStats, scale);
tree = lapply(numStats, function(x) hclust(as.dist(1-bicor(x, use = 'p', maxPOutliers = 0.05)), method = "average"))
height = 0.15;
clusters = lapply(tree, cutree, h = height); 

# 6. Plot: feature clustering dendrogram with tree cutting
pdf(file = "Plots/6_statClusteringTree.pdf", wi = 8, he = 5);
par(mar = c(4, 1, 2, 13));
  plotDendrogram(tree[[set]], horiz = TRUE, main = "Clustering of neuron shape statistics",
                 sub = "", xlab = "");
  abline(v = height, col = "red")
dev.off(); 

###############################################################################################
############ 7.8. PCA and UMAP plot, by 2 type, by 3 region, and by 11 communities #################
###############################################################################################
# PCA
pca_result = prcomp(stats1[, numericCols[[1]]], scale. = TRUE)
pca_data = data.frame(PC1 = pca_result$x[,1], 
                       PC2 = pca_result$x[,2], 
                       Type = stats1$Type,
                       Striatal_Subregion = stats1$Striatal.Subregion,
                       Striatal_Community = stats1$Striatal.Community)

# 7. Plot: PCA plots
pdf(file = "Plots/7_PCA_plots.pdf", wi = 7, he = 6);
for(cat in c('Type', 'Striatal_Subregion', 'Striatal_Community')){
  if(cat %in% c('Striatal_Subregion', 'Striatal_Community')){
    opacity = 0.6
  }else{
    opacity = 0.7
  }
  p = ggplot(pca_data, aes(x = PC1, y = PC2, color = !!sym(cat))) +
    geom_point(size = 1.5, alpha = opacity) +
    labs(title = paste0("PCA, colored by ", cat)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  
      axis.title = element_text(color = "black"),  
      axis.text = element_text(color = "black"),   
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Add conditional scale for the color
  if(cat == 'Type'){
    p = p + scale_color_manual(values = c("D1" = "darkgreen", "D2" = "darkred"))
  }else if(cat == 'Striatal_Subregion'){
    p = p + scale_color_manual(values = c("CPr" = "red", "CPi" = "green", 'CPc' = 'blue'))
  }else if (cat == 'Striatal_Community'){
    p = p + scale_color_manual(values = c("CPr_m" = "brown3", "CPr_imd" = "brown2", 'CPr_imv' = 'brown1', 'CPr_l' = 'brown',
                                          'CPi_dm' = 'chartreuse3', 'CPi_dl' = 'chartreuse2', 'CPi_vm' = 'chartreuse1', 'CPi_vl' = 'chartreuse',
                                          'CPc_d' = 'cadetblue3', 'CPc_i' = 'cadetblue2', 'CPc_v' = 'cadetblue1'))
  }
  print(p)
}
dev.off()

# same PCA plot but show one level in color at a time 
pdf(file = "Plots/7.1_PCA_plots_one_level_a_time.pdf", wi = 7, he = 6);
for(cat in c('Type', 'Striatal_Subregion')){
  if(cat %in% c('Striatal_Subregion')){
    opacity = 0.9
    levels_cat = order_subregion
  }else{
    opacity = 0.9
    levels_cat = order_type
  }
  
  if(cat == 'Type'){
    for(level in levels_cat){
      color_mapping = ifelse(pca_data[[cat]] == level, 
                             if(level == 'D1') 'red' else 'green',  
                             'grey')
       p = ggplot(pca_data, aes(x = PC1, y = PC2)) +
         geom_point(aes(color = color_mapping), size = 1.5, alpha = opacity) + 
         labs(title = paste0("PCA, highlighting ", level)) +
         theme_minimal() +
         theme(
           panel.grid = element_blank(),  
           axis.title = element_text(color = "black"),  
           axis.text = element_text(color = "black"),   
           axis.line = element_line(color = "black"),
           plot.title = element_text(hjust = 0.5)
         ) +
         scale_color_identity()
       print(p)
    }
  }
  else if(cat == 'Striatal_Subregion'){
      highlight_levels = c('CPr', 'CPi', 'CPc')
      highlight_colors = c('red', 'green', 'blue')
      
      # Loop through each level of 'highlight_levels' and create the corresponding plot
      for(i in seq_along(highlight_levels)) {
        color_mapping = ifelse(pca_data[[cat]] == 'CPr', 
                               if(i == 1) 'red' else 'grey',
                               ifelse(pca_data[[cat]] == 'CPi',
                                      if(i == 2) 'green' else 'grey',
                                      ifelse(pca_data[[cat]] == 'CPc',
                                             if(i == 3) 'blue' else 'grey',
                                             'grey')))
                                             
        # Generate the plot
        p = ggplot(pca_data, aes(x = PC1, y = PC2)) +
          geom_point(aes(color = color_mapping), size = 1.5, alpha = opacity) +  
          labs(title = paste0("PCA, highlighting ", highlight_levels[i])) +
          theme_minimal() +
          theme(
            panel.grid = element_blank(),  
            axis.title = element_text(color = "black"),  
            axis.text = element_text(color = "black"),   
            axis.line = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5)
          ) +
          scale_color_identity()  
        
        print(p)
      }
    }
}
dev.off()

# UMAP
umap_result = umap(stats1[, numericCols[[1]]] %>% scale())
umap_data = data.frame(UMAP1 = umap_result$layout[,1], 
                        UMAP2 = umap_result$layout[,2], 
                        Type = stats1$Type,
                        Striatal_Subregion = stats1$Striatal.Subregion,
                        Striatal_Community = stats1$Striatal.Community)

# 8. Plot: UMAP plots
pdf(file = "Plots/8_UMAP_plots.pdf", wi = 7, he = 6);
for(cat in c('Type', 'Striatal_Subregion', 'Striatal_Community')){
  if(cat %in% c('Striatal_Subregion', 'Striatal_Community')){
    opacity = 0.9
  }else{
    opacity = 0.7
  }
  p = ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = !!sym(cat))) +
    geom_point(size = 1.5, alpha = opacity) +
    labs(title = paste0("UMAP, colored by ", cat)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  
      axis.title = element_text(color = "black"),  
      axis.text = element_text(color = "black"),   
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Add conditional scale for the color
  if(cat == 'Type'){
    p = p + scale_color_manual(values = c("D1" = "darkgreen", "D2" = "darkred"))
  }else if(cat == 'Striatal_Subregion'){
    p = p + scale_color_manual(values = c("CPr" = "red", "CPi" = "green", 'CPc' = 'blue'))
  }else if (cat == 'Striatal_Community'){
    p = p + scale_color_manual(values = c("CPr_m" = "brown3", "CPr_imd" = "brown2", 'CPr_imv' = 'brown1', 'CPr_l' = 'brown',
                                          'CPi_dm' = 'chartreuse3', 'CPi_dl' = 'chartreuse2', 'CPi_vm' = 'chartreuse1', 'CPi_vl' = 'chartreuse',
                                          'CPc_d' = 'cadetblue3', 'CPc_i' = 'cadetblue2', 'CPc_v' = 'cadetblue1'))
  }
  print(p)
}
dev.off()

# same UMAP plot but show each level in color at a time
pdf(file = "Plots/8.1_UMAP_plots_one_level_a_time.pdf", wi = 7, he = 6);
for(cat in c('Type', 'Striatal_Subregion')){
  if(cat %in% c('Striatal_Subregion')){
    opacity = 0.9
    levels_cat = order_subregion
  }else{
    opacity = 0.9
    levels_cat = order_type
  }
  
  if(cat == 'Type'){
    for(level in levels_cat){
      color_mapping = ifelse(umap_data[[cat]] == level, 
                             if(level == 'D1') 'red' else 'green',  
                             'grey')
                             p = ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
                               geom_point(aes(color = color_mapping), size = 1.5, alpha = opacity) + 
                               labs(title = paste0("UMAP, highlighting ", level)) +
                               theme_minimal() +
                               theme(
                                 panel.grid = element_blank(),  
                                 axis.title = element_text(color = "black"),  
                                 axis.text = element_text(color = "black"),   
                                 axis.line = element_line(color = "black"),
                                 plot.title = element_text(hjust = 0.5)
                               ) +
                               scale_color_identity()
                             print(p)
    }
  }
  else if(cat == 'Striatal_Subregion'){
    highlight_levels = c('CPr', 'CPi', 'CPc')
    highlight_colors = c('red', 'green', 'blue')
    
    # Loop through each level of 'highlight_levels' and create the corresponding plot
    for(i in seq_along(highlight_levels)) {
      color_mapping = ifelse(umap_data[[cat]] == 'CPr', 
                             if(i == 1) 'red' else 'grey',
                             ifelse(umap_data[[cat]] == 'CPi',
                                    if(i == 2) 'green' else 'grey',
                                    ifelse(umap_data[[cat]] == 'CPc',
                                           if(i == 3) 'blue' else 'grey',
                                           'grey')))
                                           
      # Generate the plot
      p = ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(color = color_mapping), size = 1.5, alpha = opacity) +  
        labs(title = paste0("UMAP, highlighting ", highlight_levels[i])) +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),  
          axis.title = element_text(color = "black"),  
          axis.text = element_text(color = "black"),   
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5)
        ) +
        scale_color_identity()  
      
      print(p)
    }
  }
}
dev.off()

###############################################################################################
#################### Create Subsets and brain adjusted data ####################################
###############################################################################################
# create different subsets
colnames(stats[[1]])
# number of non-numeric columns
n_non_numeric_col = 27 

# subset 2. for female only
stats[[2]] = subset(stats[[1]], Sex == 'Female')
numStats[[2]] = numStats[[1]][stats[[1]]$Sex == 'Female',]

# subset 3.  for male only
stats[[3]] = subset(stats[[1]], Sex == 'Male')
numStats[[3]] = numStats[[1]][stats[[1]]$Sex == 'Male',]

# subset 4.  for D1 only
stats[[4]] = subset(stats[[1]], Type == 'D1')
numStats[[4]] = numStats[[1]][stats[[1]]$Type == 'D1',]

# subset 5. for D2 only
stats[[5]] = subset(stats[[1]], Type == 'D2')
numStats[[5]] = numStats[[1]][stats[[1]]$Type == 'D2',]

# subset 6. adjusted for subregion
numStats[[6]] = as.data.frame(empiricalBayesLM(numStats[[1]], removedCovariates = stats[[1]][["Striatal.Subregion"]],
                                               getEBadjustedData = FALSE)$adjustedData.OLS)
stats[[6]] = cbind(stats[[1]][1:n_non_numeric_col], numStats[[6]])

statNames = c('All', # 1
              'Female', # 2
              'Male', # 3
              'D1', # 4
              'D2', # 5
              'adj.For.Subregion' # 6
)
names(stats) = names(numStats) = statNames;
nSets=length(stats) 
print(nSets) # 6

# # calculate brain-adjusted data and save
# brain_adj_stats1 = as.data.frame(empiricalBayesLM(numStats[[1]], removedCovariates = stats[[1]][["Brain"]],
#                                                   getEBadjustedData = FALSE)$adjustedData.OLS)
# brain_adj_stats1 = cbind(stats[[1]][11:n_non_numeric_col], brain_adj_stats1)
# write.csv(brain_adj_stats1, 'morhpmetrics_brain_adjusted.csv', row.names = FALSE)

###############################################################################################
########## Output CSV for association (Type effect, Striatal Subregion effect) ################
###############################################################################################
plotTraits_association = c("Type", 
                           "Striatal.Subregion", 
                           spaste("Type_in_", order_subregion),
                           spaste("Type_in_", order_community)
)

# Function to factorize column and assign levels
convert_factors = function(df, cols, levels) {
  df[cols] = lapply(df[cols], factor, levels = levels)
  return(df)
}
stats = lapply(stats, convert_factors, cols = plotTraits_association[!plotTraits_association == 'Striatal.Subregion'], levels = order_type)

nSets = 6

test_list = list()
for (set in 1:nSets){
  
  # do not test for Striatal subregion effect in set 6
  # do not test for Type effect in set 4,5
  if(set == 6){
    plotTraits_association_test = plotTraits_association[plotTraits_association != 'Striatal.Subregion']
  }else if(set %in% c(4,5)){
    plotTraits_association_test = 'Striatal.Subregion'
  }else{
    plotTraits_association_test = plotTraits_association
  }
  
  data_list = list()
  # for each sub-group 
  for(cat in plotTraits_association_test){

    # create placeholders for shape statistics and test statistics
    statistic_list = c()
    p_value_list = c()
    p_value_adjusted_list = c()
    
    # for each shape statistic
    for(col in numericCols[[1]]){
      # use ANOVA 
      formula = as.formula(paste0(col, " ~ ", cat))
      test = aov(formula, data = stats[[set]])
      test_summary = summary(test)
      statistic = test_summary[[1]]["F value"][[1]][1]
      p_value = test_summary[[1]]["Pr(>F)"][[1]][1]
      p_value_adjusted = p.adjust(p_value, method = "bonferroni", n = total_num_trials)
  
      statistic_list = c(statistic_list, statistic)
      p_value_list = c(p_value_list, p_value)
      p_value_adjusted_list = c(p_value_adjusted_list, p_value_adjusted)
    }
    
    # convert to df
    test_df_sub = data.frame(numericCols[[1]], statistic_list, p_value_list, p_value_adjusted_list)
    colnames(test_df_sub) = c('Morphometric', paste0('statistic for ', cat), 
                              paste0('p for ', cat), paste0('p.adjusted for ', cat))
    data_list = append(data_list, list(test_df_sub))
  }
  
  # generate a big sheet for each set from data_list
  test_list[[set]] = Reduce(function(x, y) merge(x, y, by = "Morphometric", all.x = TRUE), data_list)
 
  write.csv(test_list[[set]], spaste("Results/Association_with_morphometrics", statNames[set], ".csv"), row.names = FALSE)
}

###############################################################################################
################  9. Heatmap of Normalized features by Striatal.Community ####################
###############################################################################################
# normalize features to its median by brain
normalized_df = stats[[1]] %>% group_by(Brain) %>%
  mutate(across(all_of(numericCols[[1]]), ~ . / median(., na.rm = TRUE)))

mean_all = normalized_df %>% group_by(Striatal.Community) %>% summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE))

# reorder heatmap rows by 3 main groups not the finer subgroups
# define 3 groups of features and assign colors
features_angle = c('Bif_ampl_local', 'Bif_ampl_remote', 'Bif_tilt_local', 'Bif_tilt_remote', 'Bif_torque_local', 'Bif_torque_remote', 'Centripetal_Bias')
features_length = c('Depth', 'Height', 'Width', 'Length', 'Sum_EucDistance', 'Sum_PathDistance', 'Max_EucDistance', 'Max_PathDistance', 'ABEL_All', 'ABEL_Internal', 'ABEL_Terminal', 'BAPL_All', 'BAPL_Internal', 'BAPL_Terminal')
features_complexity = c('N_bifs', 'N_branch', 'N_tips', 'N_stems', 'Branch_Order', 'Fractal_Dim', 'Partition_asymmetry', 'Terminal_degree', 'Balancing_Factor', 'Convexity')

feature_reorder = data.frame(Feature = c(features_angle, features_length, features_complexity),
                             group = rep(1:3, times = c(length(features_angle), length(features_length), length(features_complexity))))


# 9. Plot: Normalized L-measure stats distribution across Striatal.Community
dist_avg_normalized_feat_by_striatal_community = function(){
  min_color = min(mean_all[, numericCols[[1]]])
  max_color = max(mean_all[, numericCols[[1]]])
  blue_white_red = colorRamp2(c(min_color, 1, max_color), 
                              c("blue", "white", "red"))
  
  data_for_heatmap = t(mean_all[, numericCols[[1]]])
  colnames(data_for_heatmap) = order_community
  
  # reorder morphometric features
  # re-order feature by 3 main  groups
  # sort data matrix by re-ordered features
  data_for_heatmap = data_for_heatmap[match(feature_reorder$Feature, rownames(data_for_heatmap)), ]
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
                              
  ht = Heatmap(as.matrix(data_for_heatmap), 
          name = " ", 
          col = blue_white_red,
          cluster_rows = FALSE, 
          cluster_columns = FALSE, 
          show_row_names = TRUE, row_title = 'Morphometrics',
          show_column_names = TRUE, column_title = 'Striatal.Communities',
          width = unit(3, 'inch'), 
          height = unit(5, 'inch'),
          row_names_gp = gpar(fontsize= 8),
          column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
          heatmap_legend_param = list(legend_direction = "horizontal",
                                      title_position = 'topcenter'),
          right_annotation = right_annotation
  ) 
  draw(ht, heatmap_legend_side = 'top')
  grid.text(paste0('Distribution of Averaged Normalized Features: All CP neurons\n(excluded surface neurons)'),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
}

pdf(file = "Plots/9_heatmap - avg_normalized_feature_by_striatal_community.pdf", wi = 5.5, he = 7.5)
scpp(2.5);
par(mar = c(3,5,3,5));
  dist_avg_normalized_feat_by_striatal_community()
dev.off()

###############################################################################################
####### 9.1  Heatmap of Normalized features, All, 3 Subregion, 11 Communities ##############
###############################################################################################
# normalize features to its median by brain
normalized_df = stats[[1]] %>% group_by(Brain) %>%
  mutate(across(all_of(numericCols[[1]]), ~ . / median(., na.rm = TRUE)))

a = normalized_df %>% ungroup() %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  mutate(group = 'All')

b = normalized_df %>% group_by(Striatal.Subregion) %>% 
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  rename('group' = 'Striatal.Subregion')

c = normalized_df %>% group_by(Striatal.Community) %>% 
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  rename('group' = 'Striatal.Community')

mean_all = rbind(a, b, c)
mean_all$group = factor(mean_all$group, levels = c('All', order_subregion, order_community))

# 9.1: Plot: Normalized L-measure stats distribution across All, 3 subregions, 11 communities
dist_avg_normalized_feat_by_all_subregion_community = function(){
  min_color = min(mean_all[, numericCols[[1]]])
  max_color = max(mean_all[, numericCols[[1]]])
  blue_white_red = colorRamp2(c(min_color, 1, max_color), 
                              c("blue", "white", "red"))
  
  data_for_heatmap = t(mean_all[, numericCols[[1]]])
  colnames(data_for_heatmap) = levels(mean_all$group)
  
  # re-order feature by 3 main  groups
  # sort data matrix by re-ordered features
  data_for_heatmap = data_for_heatmap[match(feature_reorder$Feature, rownames(data_for_heatmap)), ]
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
  
  ht = Heatmap(as.matrix(data_for_heatmap), 
               name = " ", 
               col = blue_white_red,
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE, row_title = 'Morphometrics',
               show_column_names = TRUE, column_title = 'Striatal.Communities',
               width = unit(3, 'inch'), 
               height = unit(5, 'inch'),
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
               heatmap_legend_param = list(legend_direction = "horizontal",
                                           title_position = 'topcenter'),
               right_annotation = right_annotation
  ) 
  draw(ht, heatmap_legend_side = 'top')
  grid.text(paste0('Distribution of Averaged Normalized Features: All CP Neurons\n(excluded surface neurons)'),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
}

pdf(file = "Plots/9.1_heatmap - avg_normalized_feature_by_All_subregion_community.pdf", wi = 7, he = 7.5)
scpp(2.5);
par(mar = c(3,5,3,5));
  dist_avg_normalized_feat_by_all_subregion_community()
dev.off()

###############################################################################################
########### 10. Heatmap of Normalized feature CPr/i/c pairwise comparison by brain  ##############
###############################################################################################
# prettify brain names
replace_substrings = function(v) {
  # P56
  v = gsub("TME07-1", "P56_WT#1", v)
  v = gsub("TME08-1", "P56_WT#2", v)
  v = gsub("TME09-1", "P56_WT#3", v)
  v = gsub("TME10-1", "P56_WT#4", v)
  v = gsub("TME10-3", "P56_WT#5", v)
  return(v)
}

# calculate mean in region by brain
mean_df = normalized_df %>% group_by(Brain, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup()
# add the mean in region for all brains combined
mean_all = normalized_df %>% group_by(Striatal.Subregion) %>% summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE))
mean_all$Brain = 'All'
mean_df = rbind(mean_df, mean_all)

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
    pivot_longer(cols = -Brain, names_to = "Feature", values_to = "Value") %>%
    pivot_wider(names_from = Brain, values_from = Value) %>% as.data.frame()
  
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
# replicate data to create Brain = All observations
create_dup = function(data){
  data_dup = data.frame(data)
  orig_brain_levels = unique(data_dup$Brain)
  data_dup$Brain = 'All'
  data = rbind(data, data_dup)
  data$Brain = factor(data$Brain, levels = c('All', orig_brain_levels))
  return(data)
}

CPr_CPi_sub = normalized_df %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPi_sub = create_dup(CPr_CPi_sub)
CPr_CPc_sub = normalized_df %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPr_CPc_sub = create_dup(CPr_CPc_sub)
CPi_CPc_sub = normalized_df %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))
CPi_CPc_sub = create_dup(CPi_CPc_sub)

draw_heatmap = function(data_for_heatmap, data_subset, group_name) {
  
  # re-order the brains according to levels of data_subset
  data_for_heatmap = data_for_heatmap[, levels(data_subset$Brain)]
  # generate p-value matrix for pairwise subregion comparison per brain
  p_values = list()
  statistic_values = list()
  
  for(brain in levels(data_subset$Brain)) {
    print(brain)
    brain_data = data_subset[data_subset$Brain == brain, ] %>% as.data.frame()
    unique_region = unique(as.character(data_subset$Striatal.Subregion))
    
    # Loop through each feature and perform a t-test
    p_values[[brain]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(brain_data[brain_data$Striatal.Subregion == unique_region[1], feature],
                          brain_data[brain_data$Striatal.Subregion == unique_region[2], feature])
      test$p.value
    })
    statistic_values[[brain]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(brain_data[brain_data$Striatal.Subregion == unique_region[1], feature],
                     brain_data[brain_data$Striatal.Subregion == unique_region[2], feature])
      test$statistic
    })
  }
  
  df_pval = as.data.frame(p_values, col.names = paste0('p.value_', levels(data_subset$Brain)))
  df_statistic = as.data.frame(statistic_values, col.names = paste0('statistics_', levels(data_subset$Brain)))
  
  adjusted_p_values = lapply(p_values, function(p) p.adjust(p, method = "bonferroni", n = nNumStats))
  adjusted_p_values_df = as.data.frame(adjusted_p_values, col.names = paste0('p.value_adjusted_', levels(data_subset$Brain)))
  
  rownames(adjusted_p_values_df) = rownames(df_pval) = rownames(df_statistic) = numericCols[[1]]
  
  # write data to a csv file
  df_output = cbind(df_statistic, df_pval, adjusted_p_values_df)
  ordered_columns = as.vector(outer(c('statistics_', 'p.value_', 'p.value_adjusted_'), gsub('-', '.', levels(data_subset$Brain)), paste0))
  df_output = df_output[, ordered_columns]
  write.csv(df_output, paste0('Results/', group_name, '_statistic_output_by_individual_brain.csv'), row.names = TRUE)
  
  colnames(adjusted_p_values_df) = levels(data_subset$Brain)
  
  significance_df = apply(adjusted_p_values_df, 2, function(x) {
    ifelse(x < 0.001, "***",
           ifelse(x < 0.01, "**",
                  ifelse(x < 0.05, "*", "")))
  })
  
  data_for_heatmap = round(data_for_heatmap, 2)
  
  # re-order feature by 3 main  groups
  # sort data matrix by re-ordered features
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
  
  # prettify brain names
  colnames(data_for_heatmap) = colnames(significance_df) = replace_substrings(colnames(data_for_heatmap))
  
  # generate heatmap
  ht = Heatmap(as.matrix(data_for_heatmap), 
               name = " ", 
               col = blue_white_red,
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE, row_title = 'Morphometrics',
               show_column_names = TRUE, column_title = 'Brains',
               width = unit(3, 'inch'), 
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
  grid.text(paste0('Difference of means of\nnormalized morphometrics: ', group_name, '\n(excluded surface neurons)'),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
}

# Draw heatmaps
pdf(file = "Plots/10_heatmap - normalized_feature_region_pairwise_comparison_by_brain-pretty-name.pdf", wi = 5.5, he = 7.5)
scpp(2.5);
par(mar = c(3,5,3,5));

  draw_heatmap(CPr_CPi_heatmap_data, CPr_CPi_sub, "CPr vs. CPi")
  draw_heatmap(CPr_CPc_heatmap_data, CPr_CPc_sub, "CPr vs. CPc")
  draw_heatmap(CPi_CPc_heatmap_data, CPi_CPc_sub, "CPi vs. CPc")

dev.off()

###############################################################################################
###### 11. Heatmap of Normalized feature CPr/i/c pairwise comparison by D1D2, all brains combined  ###
###############################################################################################
# 11. Plot: normalized feature CPr vs. CPi vs. CPc by D1D2, all brains combined
mean_df = normalized_df %>% group_by(Type, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup()
# add the mean in region for all brains combined
mean_all = normalized_df %>% group_by(Striatal.Subregion) %>% summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE))
mean_all$Type = 'All'
mean_df = rbind(mean_df, mean_all)

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
reshape_for_heatmap <- function(data) {
  data = data %>%
    pivot_longer(cols = -Type, names_to = "Feature", values_to = "Value") %>%
    pivot_wider(names_from = Type, values_from = Value) %>% as.data.frame()
  
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
# replicate data to create Brain = All observations
create_dup = function(data){
  data_dup = data.frame(data)
  orig_type_levels = names(table(data_dup$Type))
  data_dup$Type = 'All'
  data = rbind(data, data_dup)
  data$Type = factor(data$Type, levels = c('All', orig_type_levels))
  return(data)
}

CPr_CPi_sub = normalized_df %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPi_sub = create_dup(CPr_CPi_sub)
CPr_CPc_sub = normalized_df %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPr_CPc_sub = create_dup(CPr_CPc_sub)
CPi_CPc_sub = normalized_df %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))
CPi_CPc_sub = create_dup(CPi_CPc_sub)

draw_heatmap = function(data_for_heatmap, data_subset, group_name, only_profound = FALSE) {
  
  # re-order the brains according to levels of data_subset
  data_for_heatmap = data_for_heatmap[, levels(data_subset$Type)]
  # generate p-value matrix for pairwise subregion comparison per brain
  p_values = list()
  statistic_values = list()
  
  for(t in levels(data_subset$Type)) {
    print(t)
    type_data = data_subset[data_subset$Type == t, ] %>% as.data.frame()
    unique_region = unique(as.character(data_subset$Striatal.Subregion))
    
    # Loop through each feature and perform a t-test
    p_values[[t]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(type_data[type_data$Striatal.Subregion == unique_region[1], feature],
                          type_data[type_data$Striatal.Subregion == unique_region[2], feature],
                          var_equal = var_equal_global)
      test$p.value
    })
    
    statistic_values[[t]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(type_data[type_data$Striatal.Subregion == unique_region[1], feature],
                     type_data[type_data$Striatal.Subregion == unique_region[2], feature],
                     var_equal = var_equal_global)
      test$statistic
    })
  }
  
  df_pval = as.data.frame(p_values, col.names = c('p.value_All', 'p.value_D1', 'p.value_D2'))
  df_statistic = as.data.frame(statistic_values, col.names = c('statistics_All', 'statistics_D1', 'statistics_D2'))
  
  adjusted_p_values = lapply(p_values, function(p) p.adjust(p, method = "bonferroni", n = nNumStats))
  adjusted_p_values_df = as.data.frame(adjusted_p_values, col.names = c('p.value_adjusted_All', 'p.value_adjusted_D1', 'p.value_adjusted_D2'))
  
  rownames(adjusted_p_values_df) = rownames(df_pval) = rownames(df_statistic) = numericCols[[1]]
  
  # write data to a csv file
  df_output = cbind(df_statistic, df_pval, adjusted_p_values_df)
  ordered_columns = as.vector(outer(c('statistics_', 'p.value_', 'p.value_adjusted_'), c('All', 'D1', 'D2'), paste0))
  df_output = df_output[, ordered_columns]
  write.csv(df_output, paste0('Results/', group_name, '_statistic_output.csv'), row.names = TRUE)
  
  colnames(adjusted_p_values_df) = levels(data_subset$Brain)
  significance_df = apply(adjusted_p_values_df, 2, function(x) {
    ifelse(x < 0.001, "***",
           ifelse(x < 0.01, "**",
                  ifelse(x < 0.05, "*", "")))
  })
  
  # re-order feature by 3 main  groups
  # filter to keep only profound features
  if(only_profound){
    feature_reorder_filtered = feature_reorder[feature_reorder$Feature %in% profound_featurs, ]
  }else{
    feature_reorder_filtered = feature_reorder
  }
  
  data_for_heatmap = round(data_for_heatmap, 2)
  
  # sort data matrix by re-ordered features
  data_for_heatmap = data_for_heatmap[match(feature_reorder_filtered$Feature, rownames(data_for_heatmap)), ]
  # sort p_significance matrix by re-ordered features
  significance_df = significance_df[match(feature_reorder_filtered$Feature, rownames(significance_df)), ]
  
  # add annotation
  annotation_df = data.frame(Cluster = feature_reorder_filtered[,'group'])
  rownames(annotation_df) = feature_reorder_filtered$Feature
  
  palette_colors = c("darkorange", "lightseagreen", "mediumorchid")
  annotation_color_mapping = setNames(palette_colors, c('1', '2', '3'))
  ### 
  right_annotation = HeatmapAnnotation(df = annotation_df, 
                                       which = "row", 
                                       annotation_name_side = NULL,
                                       show_legend = FALSE,
                                       col = list(Cluster = annotation_color_mapping))
  
  # create color map for cells
  blue_white_red <- colorRamp2(c(min_color, 0, max_color), 
                               c("blue", "white", "red"))
  
  # create data matrix for text display
  # only show significanc, do not show numbers

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
               show_column_names = TRUE, column_title = 'Type',
               width = unit(2, 'inch'), # original: 3
               height = unit(4, 'inch'), # original: 5
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
  grid.text(paste0('Difference of means of\nnormalized morphometrics: ', group_name, '\n(excluded surface neurons)'),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
}

pdf(file = "Plots/11_heatmap - normalized_feature_region_pairwise_comparison_by_D1D2_all_brains-combined.pdf", wi = 5.5, he = 7.5)
scpp(2.5);
par(mar = c(3,5,3,5));
  draw_heatmap(CPr_CPi_heatmap_data, CPr_CPi_sub, "CPr vs. CPi")
  draw_heatmap(CPr_CPc_heatmap_data, CPr_CPc_sub, "CPr vs. CPc")
  draw_heatmap(CPi_CPc_heatmap_data, CPi_CPc_sub, "CPi vs. CPc")
dev.off()

# if only keep profound features
profound_featurs = c('Bif_ampl_remote', 'Depth', 'Sum_EucDistance', 'Sum_PathDistance',
                     'Length', 'N_bifs', 'N_branch', 'N_tips', 'N_stems', 'Width',
                     'Bif_tilt_rmote', 'ABEL_All', 'BAPL_All', 'ABEL_Terminal', 'BAPL_Terminal',
                     'ABEL_Internal', 'BAPL_Internal', 'Balancing_Factor', 'Centripetal_Bias')

pdf(file = "Plots/11.1_heatmap - normalized_feature_region_pairwise_comparison_by_D1D2_all_brains-combined_profound_features.pdf", wi = 5.5, he = 7.5)
scpp(2.5);
par(mar = c(3,5,3,5));
  draw_heatmap(CPr_CPi_heatmap_data, CPr_CPi_sub, "CPr vs. CPi", TRUE)
  draw_heatmap(CPr_CPc_heatmap_data, CPr_CPc_sub, "CPr vs. CPc", TRUE)
  draw_heatmap(CPi_CPc_heatmap_data, CPi_CPc_sub, "CPi vs. CPc", TRUE)
dev.off()

###############################################################################################
#################### 12. The same heatmap further divided by D1D2  ############################
###############################################################################################
# 12. Plot: normalized feature CPr vs. CPi vs. CPc by each brain and by D1D2
mean_df_d1d2 = normalized_df %>% group_by(Brain, Type, Striatal.Subregion) %>%
  summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE)) %>%
  ungroup()
mean_df_d1d2_all = mean_df_d1d2 %>% group_by(Type, Striatal.Subregion) %>% summarise(across(all_of(numericCols[[1]]), mean, na.rm = TRUE))
mean_df_d1d2_all$Brain = 'All'
mean_df_d1d2 = rbind(mean_df_d1d2, mean_df_d1d2_all)

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
reshape_for_heatmap_d1d2 <- function(data) {
  data = data %>%
    pivot_longer(cols = -c(Brain,Type), names_to = "Feature", values_to = "Value") %>%
    pivot_wider(names_from = c(Brain,Type), values_from = Value) %>% as.data.frame()
  
  rownames(data) = data$Feature
  data = data[, -1]
}

CPr_CPi_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPr_CPi_diff_d1d2)
CPr_CPc_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPr_CPc_diff_d1d2)
CPi_CPc_heatmap_data_d1d2 = reshape_for_heatmap_d1d2(CPi_CPc_diff_d1d2)

# replicate data to create All D1 and All D2 observations
create_dup_d1d2 = function(data){
  data_dup = data.frame(data)
  orig_brain_type_levels = sort(unique(data_dup$brain_type))
  data_dup[data_dup$Type == 'D1', 'brain_type'] = 'All_D1'
  data_dup[data_dup$Type == 'D2', 'brain_type'] = 'All_D2'
  data = rbind(data, data_dup)
  data$brain_type = factor(data$brain_type, levels = c('All_D1', 'All_D2', orig_brain_type_levels))
  return(data)
}

# prepare subset data for pairwise subregion comparison within each brain
CPr_CPi_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPr", "CPi"))
CPr_CPi_sub$brain_type = paste0(CPr_CPi_sub$Brain, '_', CPr_CPi_sub$Type)
CPr_CPi_sub = create_dup_d1d2(CPr_CPi_sub)

CPr_CPc_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPr", "CPc"))
CPr_CPc_sub$brain_type = paste0(CPr_CPc_sub$Brain, '_', CPr_CPc_sub$Type)
CPr_CPc_sub = create_dup_d1d2(CPr_CPc_sub)

CPi_CPc_sub = stats[[1]] %>% filter(Striatal.Subregion %in% c("CPi", "CPc"))
CPi_CPc_sub$brain_type = paste0(CPi_CPc_sub$Brain, '_', CPi_CPc_sub$Type)
CPi_CPc_sub = create_dup_d1d2(CPi_CPc_sub)

draw_heatmap_d1d2 = function(data_for_heatmap, data_subset, group_name) {
  
  # re-order the brains according to levels of data_subset
  data_for_heatmap = data_for_heatmap[, levels(data_subset$brain_type)]
  
  # generate p-value matrix for pairwise subregion comparison per brain
  p_values = list()
  statistic_values = list()
  
  for(brain in levels(data_subset$brain_type)) {
    brain_data = data_subset[data_subset$brain_type == brain, ]
    unique_region = unique(as.character(data_subset$Striatal.Subregion))
    
    # Loop through each feature and perform a t-test
    p_values[[brain]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(brain_data[brain_data$Striatal.Subregion == unique_region[1], feature],
                          brain_data[brain_data$Striatal.Subregion == unique_region[2], feature],
                          var.equal = var_equal_global)
      test$p.value
    })
    
    statistic_values[[brain]] = sapply(numericCols[[1]], function(feature) {
      test =  t.test(brain_data[brain_data$Striatal.Subregion == unique_region[1], feature],
                     brain_data[brain_data$Striatal.Subregion == unique_region[2], feature],
                     var.equal = var_equal_global)
      test$statistic
    })
  }
  
  df_pval = as.data.frame(p_values, col.names = paste0('p.value_', levels(data_subset$brain_type)))
  df_statistic = as.data.frame(statistic_values, col.names = paste0('statistics_', levels(data_subset$brain_type)))
  
  adjusted_p_values = lapply(p_values, function(p) p.adjust(p, method = "bonferroni", n = nNumStats))
  adjusted_p_values_df = as.data.frame(adjusted_p_values, col.names = paste0('p.value_adjusted_', levels(data_subset$brain_type)))
  
  rownames(adjusted_p_values_df) = rownames(df_pval) = rownames(df_statistic) = numericCols[[1]]
  # write data to a csv file
  df_output = cbind(df_statistic, df_pval, adjusted_p_values_df)
  ordered_columns = as.vector(outer(c('statistics_', 'p.value_', 'p.value_adjusted_'), gsub('-', '.', levels(data_subset$brain_type)), paste0))
  df_output = df_output[, ordered_columns]
  write.csv(df_output, paste0('Results/', group_name, 'by_brain_by_D1D2_statistic_output.csv'), row.names = TRUE)
  
  colnames(adjusted_p_values_df) = levels(data_subset$brain_type)
  
  significance_df = apply(adjusted_p_values_df, 2, function(x) {
    ifelse(x < 0.001, "***",
           ifelse(x < 0.01, "**",
                  ifelse(x < 0.05, "*", "")))
  })
  
  # re-order feature by 3 main  groups
  # sort data matrix by re-ordered features
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
  
  # update the column names based on nameing conventions in the paper
  colnames(data_for_heatmap) = colnames(adjusted_p_values_df) = replace_substrings(colnames(data_for_heatmap))
  
  # generate heatmap
  ht = Heatmap(as.matrix(data_for_heatmap), 
               name = " ", 
               col = blue_white_red,
               cluster_rows = FALSE, 
               cluster_columns = FALSE, 
               show_row_names = TRUE, row_title = 'Morphometrics',
               show_column_names = TRUE, column_title = 'Brains',
               width = unit(5, 'inch'), 
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
  grid.text(paste0('Difference of means of\nnormalized morphometrics: ', group_name, '\n(excluded surface neurons)'),
            x = unit(0.5, 'npc'),
            y = unit(0.95, 'npc'),
            gp = gpar(fontsize = 12, fontface = 'bold'))
}

pdf(file = "Plots/12_heatmap - normalized_feature_region_pairwise_comparison_by_brain_by_D1D2_pretty-name1.pdf", wi = 8.5, he = 9)
scpp(2.5);
par(mar = c(3,5,3,5));

  draw_heatmap_d1d2(CPr_CPi_heatmap_data_d1d2, CPr_CPi_sub, "CPr vs. CPi")
  draw_heatmap_d1d2(CPr_CPc_heatmap_data_d1d2, CPr_CPc_sub, "CPr vs. CPc")
  draw_heatmap_d1d2(CPi_CPc_heatmap_data_d1d2, CPi_CPc_sub, "CPi vs. CPc")

dev.off()


###############################################################################################
###################### 13. Blue & Red Heatmap of D2 vs. D1 using t-test ########################
########################## Show mean differences ###############################################
###############################################################################################

# function to create row-wise normalization of colors
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

colored_heatmap_by_type = function(data, cat_features, set_name){
  data_matrix = matrix(NA, nrow = nNumStats)
  text_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  coef_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  
  for(cat in cat_features){
    results = data.frame(Feature = colnames(numStats[[1]]),
                         Coefficient = NA, 
                         P_Value = NA)
    
    for (feature in results$Feature) {
      
      sub_data = data[!is.na(data[, cat]), ]
      sub_data[, cat] = factor(sub_data[, cat], levels = c('D1','D2'))
      column_means = colMeans(sub_data[sub_data$Type == 'D1', numericCols[[1]]], na.rm = TRUE)
      results$Baseline_mean = column_means
      
      formula = as.formula(paste(feature, " ~ ", cat))
      test = t.test(formula, data = sub_data)
      coef = diff(test$estimate)[[1]]
      p_value = test$p.value
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
    results$Coefficient = round(results$Coefficient,3)
    
    data_matrix = cbind(data_matrix, results$neg_log10p)
    
    # conditionally round coefficients, for Fractal_Dim, Terminal_degree, Partition_asymmetry, Balancing_Factor, Convexity: 3 decimal places
    # for other features: 1 decimal places 
    small_features = c('Fractal_Dim', 'Terminal_degree', 'Partition_asymmetry', 'Balancing_Factor', 'Convexity')
    results_coef = results$Coefficient
    
    results_coef[!(results$Feature %in% small_features)] = round(results_coef[!(results$Feature %in% small_features)], 2)
    results_coef[results$Feature %in% small_features] = sprintf("%.1e", results_coef[results$Feature %in% small_features])
    
    text_matrix = cbind(text_matrix, results_coef)
    
    # if only show significance,  do not add coefficient here
    #  text_matrix = cbind(text_matrix, rep('', nNumStats))
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
  colLabs = gsub("Type", "D2 vs. D1", cat_features)
  rowLabs = results$Feature
  rownames(data_matrix) = rownames(text_matrix) = rownames(coef_matrix) = rowLabs
  colnames(data_matrix) = colnames(text_matrix) = colnames(coef_matrix) = colLabs
  
  # shade based on difference but not p-value
  color_matrix = rowwise_color_gradient(coef_matrix)
  color_matrix[data_matrix<1.301 & coef_matrix > 0] = "#FFF0F0"
  color_matrix[data_matrix<1.301 & coef_matrix < 0] = "#F0F0FF"
  
  # re-order feature by 3 main  groups
  # sort data matrix by re-ordered features
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
  
  # define cell size
  if(length(cat_features) <= 5){
    wi = 3
    he = 6
  }else if(length(cat_features) > 5 & length(cat_features) <= 8){
    wi = 9
    he = 6
  }else if(length(cat_features) > 8){
    wi = 12
    he = 6
  }
  
  ht = Heatmap(coef_matrix, 
               col = custom_color_fun,  
               name = " ",  
               show_row_names = TRUE, 
               show_column_names = TRUE,  
               cluster_rows = FALSE,  
               cluster_columns = FALSE,  
               width = unit(wi, 'inch'), 
               height = unit(he, 'inch'),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 9))
               },
               right_annotation = right_annotation,
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
               show_heatmap_legend = FALSE
  )
  
  draw(ht)
  grid.text(paste0("Coefficient and Significance Heatmap (D2 vs. D1), ", set_name, '\n(excluded surface neurons)'),
            x = unit(0.5, "npc"),
            y = unit(0.98, "npc"),  # Adjust this value to move the title up
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
  
  # output matrix of mean differences and a matrix for significance
  return(list(data_matrix = data_matrix, text_matrix = text_matrix, coef_matrix = coef_matrix))
}

# 13. Plot: heatmap of D2 vs. D1, Type and Type in subregion, show mean differences
pdf(file = spaste("Plots/13_Heatmap-Association-Type-Blue-Red-mean_diff", 'All_subregion', ".pdf"), wi = 5, he = 9);
scpp(2.5);
par(mar = c(3,5,3,5));
  cat_features = c('Type', 'Type_in_CPr', 'Type_in_CPi', 'Type_in_CPc')
  colored_heatmap_by_type(stats[[1]], cat_features, 'All') 
dev.off()

# 14. Plot: heatmap of D2 vs. D1, Type and Type in subregions and Type in communities, show mean differences
pdf(file = spaste("Plots/14_Heatmap-Association-Type-Blue-Red-mean_diff", 'All_subregion_community', "_new_feature_grouping.pdf"), wi = 15, he = 9.5);
scpp(2.5);
par(mar = c(3,5,3,5));
  cat_features = c('Type', 'Type_in_CPr', 'Type_in_CPi', 'Type_in_CPc', c(paste0('Type_in_', order_community)))
  colored_heatmap_by_type(stats[[1]], cat_features, 'All') 
dev.off()

# ### output matrix of mean differences and a matrix for significance
# result = colored_heatmap_by_type(stats[[1]], cat_features, 'All')
# data_matrix = result$data_matrix
# text_matrix = result$text_matrix
# coef_matrix = result$coef_matrix
# 
# # remove the numbers from text_matrix
# remove_numbers = function(x) {
#   gsub("[^*]", "", x)  
# }
# asterisk_matrix = apply(text_matrix, c(1,2), remove_numbers)
# 
# # save coef_matrix and asterisk_matrix 
# write.csv(asterisk_matrix, "asterisk_matrix_exclude30.csv")
# write.csv(coef_matrix, "coef_matrix_exclude30.csv")

###############################################################################################
#########################15. Blue & Red Heatmap of D2 vs. D1 using t-test ########################
############################ Show percentage change, D2 compared to D1 ###########################
###############################################################################################
# change the color gradient to based on the values in the whole matrix, not row-wise
rowwise_color_gradient_pct = function(coef_matrix) {
  color_matrix = matrix(NA, nrow = nrow(coef_matrix), ncol = ncol(coef_matrix),
                        dimnames = dimnames(coef_matrix))
  
  # Define color gradients
  red_palette = colorRampPalette(c("#FFC2C2", "#FF6565"))
  blue_palette = colorRampPalette(c("#C2C2FF", "#6565FF"))
  num_colors = 100
  red_colors = red_palette(num_colors)
  blue_colors = blue_palette(num_colors)
  
  # Map negative values to blue colors
  negative_indices = coef_matrix < 0
  if (any(negative_indices)) {
    negative_values = coef_matrix[negative_indices]
    negative_color_indices = cut(abs(negative_values), breaks = num_colors, labels = FALSE)
    color_matrix[negative_indices] = blue_colors[negative_color_indices]
  }
  
  # Map positive values to red colors
  positive_indices = coef_matrix > 0
  if (any(positive_indices)) {
    positive_values = coef_matrix[positive_indices]
    positive_color_indices = cut(positive_values, breaks = num_colors, labels = FALSE)
    color_matrix[positive_indices] = red_colors[positive_color_indices]
  }
  return(color_matrix)
}

colored_heatmap_by_type_pct = function(data, cat_features, set_name){
  data_matrix = matrix(NA, nrow = nNumStats)
  text_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  coef_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  
  for(cat in cat_features){
    results = data.frame(Feature = colnames(numStats[[1]]),
                         Coefficient = NA, 
                         P_Value = NA)
    
    for (feature in results$Feature) {
      
      sub_data = data[!is.na(data[, cat]), ]
      sub_data[, cat] = factor(sub_data[, cat], levels = c('D1','D2'))
      column_means = colMeans(sub_data[sub_data$Type == 'D1', numericCols[[1]]], na.rm = TRUE)
      results$Baseline_mean = column_means
      
      # use t-test
      formula = as.formula(paste(feature, " ~ ", 'Type'))
      test_result = t.test(formula, sub_data)
      coef = diff(test_result$estimate)[[1]]
      p_value = test_result$p.value
      if(p_correction){
        p_value = p.adjust(p_value, method = 'bonferroni', n = nNumStats)
      }
      
      results[results$Feature == feature, "Coefficient"] = coef
      results[results$Feature == feature, "P_Value"] = p_value
      # update fold change to percentage change
      results$Fold_change = results$Coefficient * 100 / results$Baseline_mean
      
    }
    
    results$P_Value = signif(results$P_Value, 2)
    results[results$P_Value == 1, 'P_Value'] = 0.95 # account for P-value = 1
    
    results$neg_log10p = -log10(results$P_Value)
    results$Coefficient = round(results$Coefficient,3)
    
    data_matrix = cbind(data_matrix, results$neg_log10p)
    
    # change text matrix to the string of percentage change 
    text_matrix = cbind(text_matrix, round(results$Fold_change,1))
    
    # change coef_matrix to the value of percentage change
    coef_matrix = cbind(coef_matrix, results$Fold_change)
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
  
  color_matrix = rowwise_color_gradient_pct(coef_matrix)
  color_matrix[data_matrix<1.301 & coef_matrix > 0] = "#FFF0F0"
  color_matrix[data_matrix<1.301 & coef_matrix < 0] = "#F0F0FF"
    
  # re-order feature by 3 main  groups
  # sort data matrix by re-ordered features
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
    
  # define cell size
  if(length(cat_features) <= 5){
    wi = 3
    he = 6
  }else if(length(cat_features) > 5 & length(cat_features) <= 8){
    wi = 9
    he = 6
  }else if(length(cat_features) > 8){
    wi = 12
    he = 6
  }
  
  ht = Heatmap(coef_matrix, 
               col = custom_color_fun,  
               name = " ",  
               show_row_names = TRUE, 
               show_column_names = TRUE,  
               cluster_rows = FALSE,  
               cluster_columns = FALSE,  
               width = unit(wi, 'inch'), 
               height = unit(he, 'inch'),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 9))
               },
               right_annotation = right_annotation,
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
               show_heatmap_legend = FALSE
  )
  
  draw(ht)
  grid.text(paste0("Percentage change and Significance Heatmap (D2 vs. D1), ", set_name, '\n(excluded surface neurons)'),
            x = unit(0.5, "npc"),
            y = unit(0.98, "npc"),  # Adjust this value to move the title up
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

# 15. Plot: heatmap of D2 vs. D1, Type and Type in subregions and Type in communities, show percentage change
pdf(file = spaste("Plots/15_Heatmap-Association-Type-Blue-Red-pct_change", 'All_subregion_community', "_new_feature_grouping.pdf"), wi = 15, he = 9.5);
scpp(2.5);
par(mar = c(3,5,3,5));
  cat_features = c('Type', 'Type_in_CPr', 'Type_in_CPi', 'Type_in_CPc', c(paste0('Type_in_', order_community)))
  colored_heatmap_by_type_pct(stats[[1]], cat_features, 'All') 
dev.off()

###############################################################################################
########## 16. Heatmap, feature region pairwise comparison, show mean differences ########
########## All, by D1D2, all brains combined ##########################################
###############################################################################################
colored_heatmap_by_region_CPr_CPi_CPc = function(data, cat_features, set_name){
  data_matrix = matrix(NA, nrow = nNumStats)
  text_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  coef_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  
  # subset data based on striatal.subregion
  if(set_name == 'CPr vs. CPi'){
    data_sub = subset(data, Striatal.Subregion %in% c('CPr', 'CPi'))
    level_order = c('CPi', 'CPr')
  }else if(set_name == 'CPr vs. CPc'){
    data_sub = subset(data, Striatal.Subregion %in% c('CPr', 'CPc'))
    level_order = c('CPc', 'CPr')
  }else if (set_name == 'CPi vs. CPc'){
    data_sub = subset(data, Striatal.Subregion %in% c('CPi', 'CPc'))
    level_order = c('CPc', 'CPi')
  }
  
  # loop through 'All', 'D1', 'D2'
  for(cat in cat_features){
    results = data.frame(Feature = colnames(numStats[[1]]),
                         Coefficient = NA, 
                         P_Value = NA)
    
    if(cat == 'All'){
      sub_data = data_sub
    }else{
      sub_data = subset(data_sub, Type == cat)
    }
    sub_data[, 'Striatal.Subregion'] = factor(sub_data[, 'Striatal.Subregion'] , levels = level_order)
    
    for (feature in results$Feature) {
      
      column_means = colMeans(sub_data[sub_data$Striatal.Subregion == level_order[1], numericCols[[1]]], na.rm = TRUE)
      results$Baseline_mean = column_means
      
      # change from mixed model to t-test
      formula = as.formula(paste(feature, " ~ ", 'Striatal.Subregion'))
      test_result = t.test(formula, sub_data)
      coef = diff(test_result$estimate)[[1]]
      p_value = test_result$p.value
      if(p_correction){
        p_value = p.adjust(p_value, method = 'bonferroni', n = nNumStats)
      }
      
      results[results$Feature == feature, "Coefficient"] = coef
      results[results$Feature == feature, "P_Value"] = p_value
      results$Fold_change = round(abs(results$Coefficient) / results$Baseline_mean,2)
    }
    
    results$P_Value = signif(results$P_Value, 2)
    results[results$P_Value == 1, 'P_Value'] = 0.95 # account for P-value = 1
    
    results$neg_log10p = -log10(results$P_Value)
    results$Coefficient = signif(results$Coefficient,3)
    
    data_matrix = cbind(data_matrix, results$neg_log10p)
    text_matrix = cbind(text_matrix, results$Coefficient)
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
  colLabs = cat_features
  rowLabs = results$Feature
  rownames(data_matrix) = rownames(text_matrix) = rownames(coef_matrix) = rowLabs
  colnames(data_matrix) = colnames(text_matrix) = colnames(coef_matrix) = colLabs
  
  # shade based on difference 
  color_matrix = rowwise_color_gradient(coef_matrix)
  color_matrix[data_matrix<1.301 & coef_matrix > 0] = "#FFF0F0"
    color_matrix[data_matrix<1.301 & coef_matrix < 0] = "#F0F0FF"
      
    # reorder features by 3 main groups
    # to reorder features and add annotation
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
    ### 
    right_annotation = HeatmapAnnotation(df = annotation_df, 
                                         which = "row", 
                                         annotation_name_side = NULL,
                                         show_legend = FALSE,
                                         col = list(Cluster = annotation_color_mapping))
    
    ht = Heatmap(coef_matrix, 
                 col = custom_color_fun,  
                 name = " ",  
                 show_row_names = TRUE, 
                 show_column_names = TRUE,  
                 cluster_rows = FALSE,  
                 cluster_columns = FALSE,  
                 width = unit(3, 'inch'), 
                 height = unit(6, 'inch'),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 9))
                 },
                 right_annotation = right_annotation,
                 row_names_gp = gpar(fontsize= 8),
                 column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
                 show_heatmap_legend = FALSE
    )
    
    draw(ht)
    grid.text(paste0("Coefficient and Significance Heatmap, ", set_name, '\n(excluded surface neurons)'),
              x = unit(0.5, "npc"),
              y = unit(0.98, "npc"),  
              gp = gpar(fontsize = 12, fontface = 'bold'))
    
    # Define the viewport to place the legend (adjust x and y as needed)
    viewport = viewport(x = 0.5, y = 0.94, width = 1, height = 0.5)
    pushViewport(viewport)
    
    # Add the custom legend
    grid.rect(gp = gpar(fill = "red"), x = 0.2, width = 0.06, height = 0.03)
    grid.text(paste0(level_order[2]," higher than ", level_order[1]), x = 0.24, just = "left", gp = gpar(fontsize = 8))
    grid.rect(gp = gpar(fill = "blue"), x = 0.5, width = 0.06, height = 0.03)
    grid.text(paste0(level_order[2]," lower than ", level_order[1]), x = 0.54, just = "left", gp = gpar(fontsize = 8))
}

# 16. Plot: heatmap of regional pairwise comparison, showing mean differences
pdf(file = "Plots/16_Heatmap-Association-Subregion-Blue-Red-mean_diffAll_D1D2.pdf", wi = 5, he = 8);
scpp(2.5);
par(mar = c(3,5,3,5));
  cat_features = c('All', 'D1', 'D2')
  colored_heatmap_by_region_CPr_CPi_CPc(stats[[1]], cat_features, 'CPr vs. CPi') 
  colored_heatmap_by_region_CPr_CPi_CPc(stats[[1]], cat_features, 'CPr vs. CPc') 
  colored_heatmap_by_region_CPr_CPi_CPc(stats[[1]], cat_features, 'CPi vs. CPc') 
dev.off()

###############################################################################################
########## 17. Heatmap, mean(std) distribution All, subregion and community for all features ########
###############################################################################################
rowwise_color_gradient_raw_distribution = function(coef_matrix) {
  color_matrix = matrix(NA, nrow = nrow(coef_matrix), ncol = ncol(coef_matrix),
                        dimnames = dimnames(coef_matrix))
  # Define color gradients
  palette = colorRampPalette(c("white", "red"))
  # Apply row-wise
  for (i in 1:nrow(coef_matrix)) {
    row_values = coef_matrix[i, ]
    max_abs_value = max(abs(row_values))
    normalized = row_values[row_values > 0] / max_abs_value
    color_matrix[i,] = palette(length(normalized))[rank(normalized)]
  }
  return(color_matrix)
}

colored_heatmap_raw_distribution = function(data, cat_features, set_name){
  data_matrix = matrix(NA, nrow = nNumStats)
  text_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  coef_matrix = matrix(NA, nrow = nNumStats, ncol = 1)
  
  for(cat in cat_features){
    if(cat == 'All'){
      sub_data = data
    }else if(cat %in% order_subregion){
      sub_data = subset(data, Striatal.Subregion == cat)
    }else if(cat %in% order_community){
      sub_data = subset(data, Striatal.Community == cat)
    }
    
    results = data.frame(Feature = colnames(numStats[[1]]),
                         Mean = NA, 
                         std = NA)
    results$Mean = round(apply(sub_data[, numericCols[[1]]], 2, mean),2)
    results$std = round(apply(sub_data[, numericCols[[1]]], 2, sd),2)
    results$Mean_std = paste0(results$Mean, ' (', results$std, ')')
    
    data_matrix = cbind(data_matrix, results$Mean_std)
    text_matrix = cbind(text_matrix, results$Mean_std)
    coef_matrix = cbind(coef_matrix, results$Mean)
  }
  
  data_matrix = data_matrix[,-1]
  text_matrix = text_matrix[,-1]
  coef_matrix = coef_matrix[,-1]
  
  # set row and column names of the matrix
  colLabs = cat_features
  rowLabs = results$Feature
  rownames(data_matrix) = rownames(text_matrix) = rownames(coef_matrix) = rowLabs
  colnames(data_matrix) = colnames(text_matrix) = colnames(coef_matrix) = colLabs
  
  # shade based on mean values
  color_matrix = rowwise_color_gradient_raw_distribution(coef_matrix)
  
  # reorder features by 3 main groups
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
  ###
  right_annotation = HeatmapAnnotation(df = annotation_df, 
                                       which = "row", 
                                       annotation_name_side = NULL,
                                       show_legend = FALSE,
                                       col = list(Cluster = annotation_color_mapping))
  
  # define cell size
  wi=18
  he = 6
  
  ht = Heatmap(coef_matrix, 
               col = custom_color_fun,  
               name = " ",  
               show_row_names = TRUE, 
               show_column_names = TRUE,  
               cluster_rows = FALSE,  
               cluster_columns = FALSE,  
               width = unit(wi, 'inch'), 
               height = unit(he, 'inch'),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 9))
               },
               right_annotation = right_annotation,
               row_names_gp = gpar(fontsize= 8),
               column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
               show_heatmap_legend = FALSE
  )
  
  draw(ht)
  grid.text(paste0("Distribution of mean (std), ", set_name, '\n(excluded surface neurons'),
            x = unit(0.5, "npc"),
            y = unit(0.98, "npc"), 
            gp = gpar(fontsize = 12, fontface = 'bold'))
}

# 17. Plot: heatmap of mean(std) for All, across subregion, and community
pdf(file = spaste("Plots/17_Heatmap-distribution_mean_std.pdf"), wi = 20, he = 8.5);
scpp(2.5);
par(mar = c(3,5,3,5));
  cat_features = c('All', order_subregion, order_community)
  colored_heatmap_raw_distribution(stats[[1]], cat_features, 'All')
dev.off()

###############################################################################################
########## 18. boxplot, comparison between striatal subregions, All and in D1 and D2 ########
###############################################################################################
# function to calculate p-values, use ANOVA
get_overall_pvalue = function(data, cat, variate) {
  
  formula = as.formula(paste(cat, " ~ ", variate))
  test = aov(formula, data = data)
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

# Function to calculate unique y positions for braces (pairwise comparison)
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

# Function to define outliers of a column and return lower&upper y limits excluding the outliers
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

# Function to generate level comparison plots with (All, D1 only and D2 only)
subregion_comparison = function(subset1, subset1_name, 
                            subset2, subset2_name,
                            subset3, subset3_name,
                            set_name){
  
  data_combined = rbind(subset1, subset2, subset3)
  data_combined$Subset = rep(c(subset1_name, subset2_name, subset3_name), 
                             times = c(nrow(subset1), nrow(subset2), nrow(subset3)))
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
      Label = paste0("Overall P value: ", p_values),
      x = 2,
      y = y_limits[2]
    )
    
    # pairwise-contrast
    results_list = lapply(split(data_combined, data_combined$Subset), function(subset_data){
    # Perform ANOVA
    formula = as.formula(paste(cat, " ~ Striatal.Subregion"))
    m1 = aov(formula, data = subset_data, na.action = na.omit)  
    
    # Post-hoc pairwise comparison after ANOVA
    emm = emmeans(m1, ~ Striatal.Subregion)  
    pairwise_comparisons = pairs(emm, adjust = "none")
    
    comparison_results = summary(pairwise_comparisons)
    comparison_results$p.value = p.adjust(comparison_results$p.value, method = 'bonferroni', n = nNumStats)
    significant_comparisons = subset(comparison_results, p.value < 0.05)
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
      theme_minimal(base_size = 14) + # Larger base font size for readability
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
                      hjust = 0.5, vjust = 0,  # Adjust horizontal and vertical position
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

# 18. Plot: Box plot, CPr-CPi-CPc, and by D1D2
pdf(spaste("Plots/18_Boxplots-by-level.pdf"), wi = 9, he = 6)
scpp(2.5);
par(mar = c(6.3, 3, 2.5, 1));
  subregion_comparison(stats[[1]], 'All', stats[[4]], 'D1 only', stats[[5]], 'D2 only', 'All')
dev.off()

###############################################################################################
############### 19. Box plot:  D1 vs. D2  Comparison (overall, by subregion) #########
###############################################################################################
comparison_var1_across_var2 = function(data, compare_var, across_var, set_name){
  
  all_data = data
  all_data[, across_var] = 'All'
  combined_data = rbind(data, all_data)
  combined_data$Type = factor(combined_data$Type, levels = order_type)
  if(across_var == 'Striatal.Subregion'){
    combined_data$Striatal.Subregion = factor(combined_data$Striatal.Subregion, levels = c('All', order_subregion))
  }

  for(cat in numericCols[[1]]){
    # count within each subgroup
    counts_data = combined_data %>%
      group_by(!!sym(across_var), !!sym(compare_var)) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(PosX = !!sym(compare_var), PosY = 0.98*min(combined_data[, cat], na.rm = TRUE))
    
    p_values = sapply(split(combined_data, combined_data[, across_var]), function(x) get_overall_pvalue(x, cat, 'Type'))
    
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

# 19. Plot: Box plot, D1 vs. D2 in All, and by subregion
pdf(paste0("Plots/19_Boxplots-D1vsD2-All_Subregion.pdf"), width = 7, height = 5)
scpp(2);
  comparison_var1_across_var2(stats[[1]], 'Type', 'Striatal.Subregion', 'All')
dev.off()
