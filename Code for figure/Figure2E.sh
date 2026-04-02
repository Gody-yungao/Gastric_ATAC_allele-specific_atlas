###########################################
################Figure2E###################
###########################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(ggplot2)
library(ggpubr)  
library(tidyr)
library(dplyr)
##
baseline_final = read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv")
################################HALLMARK  
GSVA_Hallmark = read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/GSVA/Hallmark/95sample_TMM_log2.Hallmark.GSVA_score.csv",row.names=1)
GSVA_Hallmark_trans=as.data.frame(t(GSVA_Hallmark))
GSVA_Hallmark_trans$ID=rownames(GSVA_Hallmark_trans)
#
merged_data <- baseline_final %>%   
  left_join(GSVA_Hallmark_trans, by = "ID")
#
plot_data <- merged_data[, c("Cluster",
                             "HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT")]    
#
plot_data_long <- plot_data %>%  
  pivot_longer(  
    cols = -Cluster,  
    names_to = "pathway",  
    values_to = "score"  
  )  
#
comparisons <- list(c("Cluster1", "Cluster2"),   
                   c("Cluster1", "Cluster3"),  
                   c("Cluster2", "Cluster3"))  

#
setwd("/data1/gy/ATAC_for_review/Figure2E/output")
p <- ggplot(plot_data_long, aes(x = Cluster, y = score)) +  
  geom_boxplot(  
    aes(color = Cluster),   
    fill = "white",  
    width = 0.6,  
    outlier.shape = NA  
  ) +  
  geom_point(  
    aes(color = Cluster),  
    position = position_jitterdodge(  
      jitter.width = 0,  
      jitter.height = 0,  
      dodge.width = 0.8  
    ),  
    size = 1.3,  
    alpha = 0.6  
  ) +  
  stat_compare_means(  
    comparisons = comparisons,  
    method = "wilcox.test",  
    label = "p.format",
    label.sep = " ", 
    step.increase = 0.15, 
    size = 3.5,   
    vjust = 0.5   
  ) +  
  scale_color_manual(values = c("#045A8D", "#016C59", "#B30000")) +  
  labs(  
    x = "Cluster",  
    y = "Individual score"  
  ) +  
  facet_wrap(  
    ~ pathway,   
    nrow = 1,   
    ncol = 2,  
    scales = "free_x" 
  )
  
p <- p +
    scale_y_continuous(  
    breaks = seq(-1, 1, 0.5)  
    ) +  
    coord_cartesian(ylim = c(NA, 0.8)) +  
    theme(  
    #
    panel.border = element_blank(),  
    
    #
    axis.line = element_line(color = "black", linewidth = 0.3),  
    axis.ticks = element_line(color = "black", linewidth = 0.3),  
    
    #
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),  
    strip.background = element_rect(fill = "#F0F0F0"),  
    strip.text = element_text(  
      size = 8.5,   
      face = "bold",  
      margin = margin(2,0,2,0)  
    ),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  
    axis.text.y = element_text(size = 8),  
    legend.position = "none"  
  )  
ggsave(p,filename="95sample_TMM_log2.Hallmark.GSVA_score.2IM_GC_related_pathways.boxplot.pdf",width=6,height=3.8) ##Figure2E
