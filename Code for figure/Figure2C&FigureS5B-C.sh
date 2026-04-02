#######################################################
################Figure2C&FigureS5B-C###################
#######################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(forcats)
library(ggalluvial)
library(dplyr)
baseline_final=read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv")
setwd("/data1/gy/ATAC_for_review/Figure2C&FigureS5B-C/output")
##hp
sankey_hp <- baseline_final %>%  
  mutate(  
    HP_C14 = fct_explicit_na(HP_C14, na_level = "NA")  
  ) %>%  
  group_by(Cluster) %>%  
  count(HP_C14) %>%  
  mutate(percentage = n / sum(n) * 100) %>%  
  ungroup() 

##drinking
sankey_drink <- baseline_final %>%  
  mutate(  
    drink_status = fct_explicit_na(drink_status, na_level = "NA")  
  ) %>%  
  group_by(Cluster) %>%  
  count(drink_status) %>%  
  mutate(percentage = n / sum(n) * 100) %>%  
  ungroup() 

##Tea_consumption
sankey_Tea_consumption <- baseline_final %>%  
  mutate(  
    Tea_consumption = fct_explicit_na(Tea_consumption, na_level = "NA")  
  ) %>%  
  group_by(Cluster) %>%  
  count(Tea_consumption) %>%  
  mutate(percentage = n / sum(n) * 100) %>%  
  ungroup() 

##Gastric_lesion_state
sankey_Gastric_lesion <- baseline_final %>%  
  mutate(  
    Gastric_lesion_status = fct_explicit_na(Gastric_lesion_status, na_level = "NA")  
  ) %>%  
  group_by(Cluster) %>%  
  count(Gastric_lesion_status) %>%  
  mutate(percentage = n / sum(n) * 100) %>%  
  ungroup() 

###################color
drink_colors <- c(
  "Drinkers" = "#4575B4", 
  "Nondrinkers" = "#FFD700",
  "NA" = "#CCCCCC"
)
hp_colors <- c(
  "Positive" = "#FFA500", 
  "Negative" = "#B2B2B2" 
)
Gastric_lesion_colors <- c(
  "SG" = "#FFE4B5", 
  "AG" = "#FDAE6B",
  "IM" = "#F8766D"
)
Tea_consumption_colors <- c(
  "Tea-drinkers" = "#5B9BD5",
  "Nondrinkers" = "#C7E9B4", 
  "NA" = "#CCCCCC"
)

# 绘制 alluvial 图
##hp
sankey_hp <- sankey_hp %>%  
  mutate(HP_C14 = factor(HP_C14,  
                               levels = c("Negative", "Positive"),  
                               labels = c("Negative", "Positive")))
p_hp=ggplot(sankey_hp, aes(x = Cluster, y = percentage, fill = HP_C14, stratum = HP_C14, alluvium = HP_C14)) +
    geom_alluvium(width=0.5,alpha =0.3,knot.pos=0.2,color='white') +
    geom_col(width = 0.5,color='white') + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_discrete(limits = c("Cluster1","Cluster2","Cluster3")) + 
    scale_fill_manual(values = hp_colors, breaks = c("Positive","Negative")) + 
    xlab("") + ylab("Percentage(%)") + 
    theme_classic() + 
    theme(
      #
      axis.text.x = element_text(  
      size = 10,  
      color = "black",  
      angle = 45, 
      hjust = 1 
      ),  
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      plot.margin = margin(15, 15, 15, 15), 
      panel.grid = element_blank() 
)
ggsave(p_hp,filename="95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.Hp_cluster.Alluvialplot.pdf") ##Figure2C

##drinking
sankey_drink <- sankey_drink %>%  
  mutate(drink_status = factor(drink_status,  
                               levels = c("NA", "Nondrinkers", "Drinkers"),  
                               labels = c("NA", "Nondrinkers", "Drinkers")))
p_drink=ggplot(sankey_drink, aes(x = Cluster, y = percentage, fill = drink_status, stratum = drink_status, alluvium = drink_status)) +
    geom_alluvium(width=0.5,alpha =0.3,knot.pos=0.2,color='white') +
    geom_col(width = 0.5,color='white') + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_discrete(limits = c("Cluster1","Cluster2","Cluster3")) +
    scale_fill_manual(values = drink_colors, breaks = c("Drinkers","Nondrinkers","NA")) + 
    xlab("") + ylab("Percentage(%)") + 
    theme_classic() + 
    theme(
      #
      axis.text.x = element_text(  
      size = 10,  
      color = "black",  
      angle = 45, 
      hjust = 1 
      ),  
      axis.text.y = element_text(size = 10, color = "black"), 
      axis.title.y = element_text(size = 12, margin = margin(r = 10)), 
      legend.text = element_text(size = 10), 
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      plot.margin = margin(15, 15, 15, 15), 
      panel.grid = element_blank() 
)
ggsave(p_drink,filename="95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.drink_cluster.Alluvialplot.pdf") ##FigureS5B

##Gastric_lesion
sankey_Gastric_lesion <- sankey_Gastric_lesion %>%  
  mutate(Gastric_lesion_status = factor(Gastric_lesion_status,  
                               levels = rev(c("SG", "AG", "IM")), 
                               labels = rev(c("SG", "AG", "IM"))))
p_Gastric_lesion=ggplot(sankey_Gastric_lesion, aes(x = Cluster, y = percentage, fill = Gastric_lesion_status, stratum = Gastric_lesion_status, alluvium = Gastric_lesion_status)) + 
    geom_alluvium(width=0.5,alpha =0.3,knot.pos=0.2,color='white') + 
    geom_col(width = 0.5,color='white') + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_x_discrete(limits = c("Cluster1","Cluster2","Cluster3")) + 
    scale_fill_manual(values = Gastric_lesion_colors, breaks = c("SG", "AG", "IM")) +
    xlab("") + ylab("Percentage(%)") + 
    theme_classic() + 
    theme(
      #
      axis.text.x = element_text(  
      size = 10,  
      color = "black",  
      angle = 45,
      hjust = 1 
      ),  
      axis.text.y = element_text(size = 10, color = "black"),  
      axis.title.y = element_text(size = 12, margin = margin(r = 10)), 
      legend.text = element_text(size = 10), 
      legend.title = element_blank(), 
      legend.key.size = unit(1, "cm"),
      plot.margin = margin(15, 15, 15, 15), 
      panel.grid = element_blank() 
)
ggsave(p_Gastric_lesion,filename="95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.Gastric_lesion_cluster.Alluvialplot.pdf") ##FigureS5C

###################################################################################
############fisher tsest: Cluster3 compared with the other 2 clusters##############
###################################################################################
###################################################Hp
# 
fisher_vs <- function(df, ref) {  
  # 2x2 matrix 
  cluster3_data <- df %>% filter(Cluster == "Cluster3")  
  ref_data <- df %>% filter(Cluster == ref)  
  
  #
  mat <- matrix(  
    c(cluster3_data$n[cluster3_data$HP_C14 == "Positive"],  
      ref_data$n[ref_data$HP_C14 == "Positive"],  
      cluster3_data$n[cluster3_data$HP_C14 == "Negative"],  
      ref_data$n[ref_data$HP_C14 == "Negative"]),  
    nrow = 2,  
    byrow = FALSE,  
    dimnames = list(  
      c("Cluster3", ref),  
      c("Positive", "Negative")  
    )  
  )  
  
  # 
  ft <- fisher.test(mat)  
  
  #
  tibble::tibble(  
    comparison = paste("Cluster3 vs", ref),  
    OR = ft$estimate,  
    CI_low = ft$conf.int[1],  
    CI_high = ft$conf.int[2],  
    p = ft$p.value
  )  
}  

#
result <- bind_rows(  
  fisher_vs(sankey_hp, "Cluster1"),  
  fisher_vs(sankey_hp, "Cluster2")  
)  

#
print(as.data.frame(result))
#            comparison        OR    CI_low   CI_high            p
#1 Cluster3 vs Cluster1 66.444907 10.848776 798.99497 1.647600e-09
#2 Cluster3 vs Cluster2  7.435361  1.474063  74.03216 8.548188e-03

###################################################drink
#
sankey_drink_clean <- sankey_drink %>%
  filter(drink_status != "NA") %>%
  group_by(Cluster, drink_status) %>%
  summarise(n = sum(n), .groups = "drop")

# 
fisher_vs_cluster1 <- function(df, target_cluster) {
  #  2×2 matrix
  ref_data <- df %>% filter(Cluster == "Cluster1")
  target_data <- df %>% filter(Cluster == target_cluster)
  
  mat <- matrix(
    c(
      ref_data$n[ref_data$drink_status == "Drinkers"],
      ref_data$n[ref_data$drink_status == "Nondrinkers"],
      target_data$n[target_data$drink_status == "Drinkers"], 
      target_data$n[target_data$drink_status == "Nondrinkers"] 
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Cluster1", target_cluster),
      c("Drinkers", "Nondrinkers")
    )
  )
  
  #
  ft <- fisher.test(mat)
  
  #
  tibble::tibble(
    comparison = paste("Cluster1 vs", target_cluster),
    OR = ft$estimate,
    CI_low = ft$conf.int[1],
    CI_high = ft$conf.int[2],
    p = ft$p.value,
  )
}

# 
result <- bind_rows(
  fisher_vs_cluster1(sankey_drink_clean, "Cluster2"),
  fisher_vs_cluster1(sankey_drink_clean, "Cluster3")
)
#
print(as.data.frame(result))
#            comparison         OR      CI_low   CI_high           p
#1 Cluster1 vs Cluster2 0.07206617 0.001590206 0.5466114 0.002127460
#2 Cluster1 vs Cluster3 0.05657866 0.001195124 0.4710467 0.001230132

#################################################################
############prop.trend.test: Cluster1-3 SG/IM trend##############
#################################################################
############################IM
IM_counts <- c(1, 4, 5)
total_counts <- c(32, 39, 24)
trend_test <- prop.trend.test(IM_counts, total_counts)
trend_test
trend_test
#	Chi-squared Test for Trend in Proportions
#
#data:  IM_counts out of total_counts ,
# using scores: 1 2 3
#X-squared = 4.4995, df = 1, p-value = 0.0339

############################SG
SG_counts <- c(19, 21, 8)
total_counts <- c(32, 39, 24)
trend_test <- prop.trend.test(SG_counts, total_counts)
trend_test
trend_test
#	Chi-squared Test for Trend in Proportions
#
#data:  SG_counts out of total_counts ,
# using scores: 1 2 3
#X-squared = 3.5005, df = 1, p-value = 0.06135
