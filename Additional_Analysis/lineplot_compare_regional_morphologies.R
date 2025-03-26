# Set work directory
setwd('...')
# source functions
func_dir = '../Functions/'

source(paste0(func_dir, "individualAnalysis-General-010.R"));
source(paste0(func_dir, "GNVFunctions-018-02.R"));
source(paste0(func_dir, "networkFunctions-extras-20.R"));
source(paste0(func_dir, "labelPoints2-01.R"));
source(paste0(func_dir, "heatmap.wg.R"));

library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)

########################## READ AND COMBINE DATA #########################################
# read both 2-m and 12-m data (CP neurons only)
# read 2-m data, CP neurons only
data_2m = read.csv('morpho_2248CPneurons_with_community.csv')
# read 12-m data, CP neurons only
data_12m = read.csv('all_htme_brains_with_registration_1168CPneurons_onlyCP.csv')
# create group indicator
data_2m$group = 'P56 WT'
data_12m$group = ifelse(data_12m$Genotype == 'WT', '12m WT', '12m Q140')

# drop uncommon columns for both datasets
col_not_common = setdiff(union(colnames(data_2m), colnames(data_12m)), intersect(colnames(data_2m), colnames(data_12m)))
data_2m = data_2m[, !(names(data_2m) %in% col_not_common)]
ncol(data_2m)
data_12m = data_12m[, !(names(data_12m) %in% col_not_common)]
ncol(data_12m)

# combine data
data_combined = bind_rows(data_2m, data_12m)
dim(data_combined) # 3415 rows, 36 columns

# define levels and orders
order_type = c('D1', 'D2')
order_subregion = c('CPr','CPi','CPc')
order_group = c('P56 WT','12m WT', '12m Q140')

data_combined$Type = factor(data_combined$Type, levels = order_type)
data_combined$Striatal.Subregion = factor(data_combined$Striatal.Subregion, levels = order_subregion)
data_combined$group = factor(data_combined$group, levels = order_group)

colnames(data_combined)
numeric_features = colnames(data_combined[5:35])
length(numeric_features) # 31

nNumStats = length(numeric_features)

# normalize features to its median by brain
normalized_df = data_combined %>% group_by(Brain) %>%
  mutate(across(all_of(numeric_features), ~ . / median(., na.rm = TRUE)))

# Define the desired order and grouping of the features
desired_feature_order = c("Bif_ampl_local", "Bif_ampl_remote", "Bif_tilt_local", "Bif_tilt_remote", "Bif_torque_local", "Bif_torque_remote", 'Centripetal_Bias',
                          "Length", "Sum_EucDistance", "Sum_PathDistance", "Max_EucDistance", "Max_PathDistance", "ABEL_All", "BAPL_All", "ABEL_Terminal", "BAPL_Terminal", "ABEL_Internal", "BAPL_Internal", "Height", "Width", "Depth",
                          "N_stems", "N_branch", "N_bifs", "N_tips", "Terminal_degree", "Branch_Order", "Fractal_Dim", "Partition_asymmetry", "Balancing_Factor", "Convexity")

#############################################################################################
#### Output statistics for group differences of normalized features (ANOVA, posthoc) ########
######### for p-value annotation only, and full statistics #####################
############################################################################################################
# function for p-value labeling
significance_stars = function(p_values) {
  sapply(p_values, function(p) {
    if (is.na(p)) {
      return(NA)      # Return NA if p is NA
    } else if (p < 0.001) {
      return("***") 
    } else if (p < 0.01) {
      return("**")   
    } else if (p < 0.05) {
      return("*")   
    } else {
      return("")     
    }
  })
}

# specify output column names to generate clean column names
output_columns = c("12m Q140-12m WT_diff",  "12m Q140-12m WT_lwr",   "12m Q140-12m WT_p adj", "12m Q140-12m WT_upr",   "12m Q140-P56 WT_diff", 
                   "12m Q140-P56 WT_lwr",   "12m Q140-P56 WT_p adj", "12m Q140-P56 WT_upr",   "12m WT-P56 WT_diff",    "12m WT-P56 WT_lwr",    
                   "12m WT-P56 WT_p adj",   "12m WT-P56 WT_upr" )

# function to generate and organize result statistics
get_group_difference_anova = function(df, convert_p_value = TRUE){
  result_df = data.frame(Features = desired_feature_order)
  
  adhoc_final_all = list()
  
  for(t in c('D1', 'D2')){
    data_sub = subset(df, Type == t)
    
    statistics_vec = c()
    p_val_vec = c()
    p_val_adjusted_vec = c()
    
    adhoc_test_df_list = list()
    
    for(col in desired_feature_order){
      formula = as.formula(paste0(col, ' ~ ', 'group'))
      test = aov(formula, data = data_sub)
      
      statistic =  summary(test)[[1]]['F value'][[1]][1]
      p_val = summary(test)[[1]]["Pr(>F)"][[1]][1]
      p_val_adjusted = p.adjust(p_val, method = 'bonferroni', n = nNumStats)
      
      statistics_vec = c(statistics_vec, statistic)
      p_val_vec = c(p_val_vec, p_val)
      p_val_adjusted_vec = c(p_val_adjusted_vec, p_val_adjusted)
      
      # only do post-hoc if the anova is significant
      if(p_val_adjusted <= 0.05){
        # ad-hoc analysis, Tukey's HSD
        post_hoc = TukeyHSD(test)
        tukey_result = as.data.frame(post_hoc$group)
        tukey_result$group_comparison = rownames(tukey_result)
        
        tukey_long = tukey_result %>%
          gather(key = "metric", value = "value", -group_comparison)
        
        # Add a column to distinguish each group comparison and metric
        tukey_long$metric_group = paste0(tukey_long$group_comparison, "_", tukey_long$metric)
        
        # Pivot the data into wide format
        tukey_wide = tukey_long %>%
          select(metric_group, value) %>%
          spread(key = metric_group, value = value)
        
        adhoc_test_df_list[[col]] = tukey_wide
      }
    }
    
    if(convert_p_value){
      p_val_vec = significance_stars(p_val_vec)
      p_val_adjusted_vec = significance_stars(p_val_adjusted_vec)
    }
    
    result_df[, paste0(t, '_statistic')] = statistics_vec
    result_df[, paste0(t, '_p')] = p_val_vec
    result_df[, paste0(t, '_p_adjusted')] = p_val_adjusted_vec
    
    # aggregate adhoc data
    adhoc_final = bind_rows(adhoc_test_df_list, .id = 'Morphometrics')
    
    if(nrow(adhoc_final > 0)){
      # add missing morphometrics in the adhoc result
      rownames(adhoc_final) = adhoc_final[, 'Morphometrics']
      adhoc_final = adhoc_final[, -1]
      # further adjust for multi-feature comparisons
      num_to_adjust = nrow(adhoc_final)
      adhoc_final[grepl("p adj", names(adhoc_final))] = adhoc_final[grepl("p adj", names(adhoc_final))] * num_to_adjust
      
      # add missing features
      missing_col = setdiff(desired_feature_order, rownames(adhoc_final))
      num_empty_rows = length(missing_col)
      empty_rows = data.frame(matrix(NA, ncol = ncol(adhoc_final), nrow = num_empty_rows))
      colnames(empty_rows) = colnames(adhoc_final)
      rownames(empty_rows) = missing_col
      
      adhoc_final = rbind(adhoc_final, empty_rows)
      adhoc_final = adhoc_final[match(desired_feature_order, rownames(adhoc_final)), ]
      
      # rename based on D1 or D2
      adhoc_final = adhoc_final %>% rename_with(~ paste0(t, '_', .))
    }else{
      adhoc_final = data.frame(matrix(NA, ncol = length(output_columns), nrow = length(desired_feature_order)))
      colnames(adhoc_final) = output_columns
      adhoc_final = adhoc_final %>% rename_with(~ paste0(t, '_', .))
    }
    
    adhoc_final_all[[t]] = adhoc_final
  }
  
  adhoc_final_all_df = bind_cols(adhoc_final_all)
  
  return(list(result_df = result_df, adhoc_final = adhoc_final_all_df))
}

# 1) with statistic and un-converted p-values
results_cpr1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPr'), FALSE)
results_cpi1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPi'), FALSE)
results_cpc1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPc'), FALSE)

# output adhoc comparison, adhoc with statistics
wb = createWorkbook()
  addWorksheet(wb, "CPr")
  addWorksheet(wb, "CPi")
  addWorksheet(wb, "CPc")
  
  writeData(wb, sheet = 'CPr', results_cpr1$adhoc_final, rowNames = TRUE)
  writeData(wb, sheet = 'CPi', results_cpi1$adhoc_final, rowNames = TRUE)
  writeData(wb, sheet = 'CPc', results_cpc1$adhoc_final, rowNames = TRUE)
saveWorkbook(wb, "differences_in_subregion_by_D1D2_adhoc_comparison.xlsx", overwrite = TRUE)

# 2) adhoc comparison with p-value asterisk labels for plot annotation
input_file = "differences_in_subregion_by_D1D2_adhoc_comparison.xlsx" 
wb = loadWorkbook(input_file)
  sheet_names = getSheetNames(input_file)
  new_wb = createWorkbook()
  
  # for each subregion
  for (sheet in sheet_names) {
    df = read.xlsx(input_file, sheet = sheet)
    colnames(df)[1] = "Morphometric"
    
    # Keep only columns that contain 'p adj' in names
    p_adj_cols = grep("p.adj", colnames(df), value = TRUE)
    
    if (length(p_adj_cols) > 0) {
      df[p_adj_cols] = lapply(df[p_adj_cols], significance_stars)
    }
    addWorksheet(new_wb, sheet) 
    writeData(new_wb, sheet, df[, c('Morphometric', p_adj_cols)], rowNames = FALSE) 
  }
  output_file = "differences_in_subregion_by_D1D2_adhoc_comparison_asterisk.xlsx"  
saveWorkbook(new_wb, output_file, overwrite = TRUE)


## ANOVA significance result, asterisk only
results_cpr1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPr'), TRUE)
results_cpi1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPi'), TRUE)
results_cpc1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPc'), TRUE)

wb = createWorkbook()
addWorksheet(wb, "CPr")
addWorksheet(wb, "CPi")
addWorksheet(wb, "CPc")

writeData(wb, sheet = 'CPr', results_cpr1$result_df, rowNames = TRUE)
writeData(wb, sheet = 'CPi', results_cpi1$result_df, rowNames = TRUE)
writeData(wb, sheet = 'CPc', results_cpc1$result_df, rowNames = TRUE)
saveWorkbook(wb, "significance_by_region_3sets_ANOVA_asterisks.xlsx", overwrite = TRUE)

## ANOVA significance result, original statistics
results_cpr1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPr'), FALSE)
results_cpi1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPi'), FALSE)
results_cpc1 = get_group_difference_anova(subset(normalized_df, Striatal.Subregion == 'CPc'), FALSE)

wb = createWorkbook()
addWorksheet(wb, "CPr")
addWorksheet(wb, "CPi")
addWorksheet(wb, "CPc")

writeData(wb, sheet = 'CPr', results_cpr1$result_df, rowNames = TRUE)
writeData(wb, sheet = 'CPi', results_cpi1$result_df, rowNames = TRUE)
writeData(wb, sheet = 'CPc', results_cpc1$result_df, rowNames = TRUE)
saveWorkbook(wb, "differences_ANOVA_statistics.xlsx", overwrite = TRUE)

############################################################################################################
######### Line plot: normalized feature distribution in CPr, CPi and CPc, by D1 and D2 #####################
############################################################################################################
pdf(file = "lineplot_regional_differences_3sets.pdf", wi = 7, he = 9)
scpp(2.5);
par(mar = c(3,5,3,5));
  for(grp in c('All', 'D1', 'D2')){
    if(grp %in% c('D1', 'D2')){
      normalized_df_plot = normalized_df %>% filter(Type == grp)
    }else{
      normalized_df_plot = normalized_df
    }
    
    mean_all = normalized_df_plot %>% 
      group_by(group, Striatal.Subregion) %>% 
      summarise(across(all_of(numeric_features), mean, na.rm = TRUE))
    
    mean_all_striatum = normalized_df_plot %>% 
      group_by(group) %>% 
      summarise(across(all_of(numeric_features), mean, na.rm = TRUE))
    
    # reshape data to create the desired line plot
    df_long = pivot_longer(mean_all, 
                           cols = numeric_features,  
                           names_to = "Feature", 
                           values_to = "Value")
    df_long_all_striatum = pivot_longer(mean_all_striatum, 
                                        cols = numeric_features,  
                                        names_to = "Feature", 
                                        values_to = "Value")
    df_long_all_striatum$Striatal.Subregion = 'All striatum'
    df_long_combined = bind_rows(df_long, df_long_all_striatum)
    
    # Reorder features
    df_long_combined$Feature = factor(df_long_combined$Feature, levels = desired_feature_order)
    df_long_combined$Striatal.Subregion = factor(df_long_combined$Striatal.Subregion, levels = c('All striatum', order_subregion))
    
    # facet by All/CPr/i/c
    p = ggplot(df_long_combined, aes(x = Feature, y = Value, group = group)) +
      geom_line(size = 0.8, aes(color = group)) +
      facet_wrap(~ Striatal.Subregion, ncol = 1) + 
      scale_color_manual(values = c("12m Q140" = "#FF33FF", "12m WT" = "#00cccc", "P56 WT" = "#ff9933")) + 
      labs(title = paste0("Line plot of normalized features, ", grp),
           x = "Features",
           y = "Averaged normalized values",
           color = "Group") +
      ylim(0.8, 1.25) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "black") + 
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        strip.text = element_text(face = "bold")
      )
    print(p)
  }
dev.off()

############################################################################################################
###################  Line plot: normalized feature regional pairwise differences ###########################
#################### (CPr vs. CPi, CPr vs. CPc, CPi vs. CPc) ################################################ 
#########                         by D1 and D2                           #####################
############################################################################################################
# show the regional pairwise differences instead of the raw normalized values
pdf(file = "line-plot_averaged_normalized_feature_subregion_pairwise_P56_12m.pdf", wi = 7, he = 9)
scpp(2.5);
par(mar = c(3,5,3,5));
  for(grp in c('All', 'D1', 'D2')){
    if(grp %in% c('D1', 'D2')){
      normalized_df_plot = normalized_df %>% filter(Type == grp)
    }else{
      normalized_df_plot = normalized_df
    }
    mean_all = normalized_df_plot %>% 
      group_by(group, Striatal.Subregion) %>% 
      summarise(across(all_of(numeric_features), mean, na.rm = TRUE))
    
    # reshape data to create the desired line plot
    df_long = pivot_longer(mean_all, 
                           cols = numeric_features,  
                           names_to = "Feature", 
                           values_to = "Value")
    
    differences = df_long %>%
      group_by(group, Feature) %>%
      summarise(
        CPr_CPi = Value[Striatal.Subregion == "CPr"] - Value[Striatal.Subregion == "CPi"],
        CPr_CPc = Value[Striatal.Subregion == "CPr"] - Value[Striatal.Subregion == "CPc"],
        CPi_CPc = Value[Striatal.Subregion == "CPi"] - Value[Striatal.Subregion == "CPc"],
        .groups = 'drop' 
      )
    
    differences_long = differences %>%
      pivot_longer(
        cols = c(CPr_CPi, CPr_CPc, CPi_CPc),  
        names_to = "subregion_pair",              
        values_to = "Value"                  
      )
    
    # rename CPr_CPi, CPr_CPc and CPi_CPc
    differences_long = differences_long %>%
      mutate(subregion_pair = case_when(
        subregion_pair == "CPr_CPi" ~ "CPr vs. CPi",
        subregion_pair == "CPr_CPc" ~ "CPr vs. CPc",
        subregion_pair == "CPi_CPc" ~ "CPi vs. CPc",
        TRUE ~ subregion_pair  
      ))
    
    differences_long$Feature = factor(differences_long$Feature, levels = desired_feature_order)
    differences_long$subregion_pair = factor(differences_long$subregion_pair, levels = c('CPr vs. CPi', 'CPr vs. CPc', 'CPi vs. CPc'))
    
    p = ggplot(differences_long, aes(x = Feature, y = Value, group = group)) +
      geom_line(size = 1.2, aes(color = group)) +
      facet_wrap(~ subregion_pair, ncol = 1) +  
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      ylim(-0.4, 0.4) + 
      scale_color_manual(values = c("12m Q140" = "#FF33FF", "12m WT" = "#00cccc", "P56 WT" = "#ff9933")) +  
      labs(title = paste0("Line plot of difference of averaged normalized features, ", grp),
           x = "Features",
           y = "Averaged normalized values",
           color = "subregion_pair") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        strip.text = element_text(face = "bold")
      )
    print(p)
  }
dev.off()

