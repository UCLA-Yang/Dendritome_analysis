##########################################################################################
#### create bar plot distribution of cell count for each brain ########################
#### also show the proportion of D1/D2 ################################################
##########################################################################################
# work directory
setwd("C:/Users/yanyanming77/Desktop/neuro_shape_analysis/Additional_analysis_for_paper_cleaned")

library(dplyr)
library(ggplot2)


data_2m = read.csv('C:/Users/yanyanming77/Desktop/Lmeasure_extraction/D1_D2/Lmeasure_extract_all_features/May-20-2024/morpho_2466CPneurons_with_community.csv')
nrow(data_2m) # 2466
data_12m = read.csv('C:/Users/yanyanming77/Desktop/Lmeasure_extraction/Q140_WT_withD1/Lmeasure_extract_all_features/all_htme_brains_with_registration_all.csv')
nrow(data_12m) # 1296

colnames(data_2m)

data_2m$Genotype = 'WT'
data_2m$Age = 'P56'

# create Brain_ind variable for brain-genotype-age
data_2m = data_2m %>%
  select(file_path, Brain, Type, Genotype, Age) %>%
  mutate(Brain_ind = paste0(Brain, ' (', Genotype, ' ', Age, ')'))

data_12m$Age = '12-m'
data_12m = data_12m %>%
  select(file_path, Brain, Type, Genotype, Age) %>%
  mutate(Brain_ind = paste0(Brain, ' (', Genotype, ' ', Age, ')'))

# combine all data
data_combined = rbind(data_2m, data_12m)
table(data_combined$Brain_ind)

# define order
order_Brain_ind = c('TME07-1 (WT P56)', 'TME08-1 (WT P56)', 'TME09-1 (WT P56)',  'TME10-1 (WT P56)', 'TME10-3 (WT P56)',
                    'hTME15-1 (WT 12-m)', 'hTME16-1 (WT 12-m)', 'hTME19-2 (WT 12-m)',
                    'hTME15-2 (Q140 12-m)',  'hTME18-1 (Q140 12-m)', 'hTME20-1 (Q140 12-m)', 'hTME24-2 (Q140 12-m)')
data_combined$Brain_ind = factor(data_combined$Brain_ind, levels = order_Brain_ind)

# calculate group size
df_summary = data_combined %>%
  group_by(Brain_ind, Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Brain_ind) %>%
  ungroup() 

total_cnt = df_summary %>% 
  group_by(Brain_ind)%>% 
  summarise(Total = sum(Count), .groups = 'drop')
df_summary = total_cnt %>% 
  left_join(df_summary, by = c("Brain_ind"))

# generate distribution bar plot
pdf(paste0("Type_Distribution_Brain_all_3762_MSNs.pdf"), wi = 9, he = 6)
p = ggplot(df_summary, aes(x = Brain_ind, y = Count, fill = Type, label = Count)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_text(position = position_stack(vjust = 0.5), size = 3, color = "white") + 
  geom_text(aes(y = Total, label=Total), vjust = -0.5, size = 3.5, color = "black") + 
  # facet_wrap(~Genotype) + 
  scale_fill_manual(values = c("D1" = "darkgreen", "D2" = "darkred")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() + 
  labs(title = "Distribution of Type by Brain", 
       x = "Brain", 
       y = "Count") +
  theme(
    plot.title = element_text(size = 16, face = 'bold', hjust=0.5), 
    axis.title = element_text(size = 14), 
    axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
    strip.text = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),  
    axis.line = element_line(color = "black") 
  )
print(p)
dev.off()