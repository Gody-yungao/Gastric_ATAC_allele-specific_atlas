##################################
############Figure2H##############
##################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(ggplot2)
library(ggpubr)  
library(tidyr)
library(dplyr)
##
gtf_v19 <- rtracklayer::import('/data1/gy/public/gtf/gencode.v19.annotation.gtf')
head(gtf_v19)
gtf_v19_df <- as.data.frame(gtf_v19)
geneid_df <- dplyr::select(gtf_v19_df,c(gene_name,gene_id))#,gene_biotype
geneid_df<-unique(geneid_df)
geneid_df$new <- geneid_df$gene_id
new_gtf_v19<-geneid_df
new_gtf_v19$gene_id<-sapply(stringr::str_split(new_gtf_v19$gene_id, "\\."), function(v)  return(v[1]))
new_gtf_v19<-unique(new_gtf_v19)
new_gtf_v19 <- new_gtf_v19[,c(1,2)]
##
TMM_norm <- as.data.frame(fread("/data1/gy/multistage_RNAseq/xuxianfeng/CombineCounts.FilterLowExpression-MergeMutiSample.TMM.tsv"))
TMM_norm=TMM_norm[,c(133,1:132)]
#
colnames(TMM_norm) <- gsub("_DGC|_IGC", "_GC", colnames(TMM_norm))
TMM_norm=merge(new_gtf_v19,TMM_norm,by="gene_id")
TMM_norm=TMM_norm[,-1]

##
condition <- sapply(strsplit(colnames(TMM_norm)[-1] , "_") , "[" , 2)
table(condition)
#condition
#    GC     IM Normal 
#    48     48     36
## 
TMM_norm_mat <- TMM_norm[,-1]
rownames(TMM_norm_mat) <- TMM_norm$gene_name

##########################################################################################
class_type <- c( "Normal" , "IM" , "GC")
TMM_norm_mat_trans=as.data.frame(t(TMM_norm_mat))
TMM_norm_mat_trans$ID=rownames(TMM_norm_mat_trans)
TMM_norm_mat_trans$lesion=condition
dim(TMM_norm_mat_trans)
#[1]   132 14752
TMM_norm_mat_trans=TMM_norm_mat_trans[,c(14751:14752,1:14750)]

# 
plot_data <- TMM_norm_mat_trans[, c("lesion",
                             "NRF1","SP2","IRF8",
                             "ZNF416","GMEB2","SPI1")]    
# 
plot_data_long <- plot_data %>%  
  pivot_longer(  
    cols = -lesion,  
    names_to = "TF",  
    values_to = "TMM"  
  )  
# 
comparisons <- list(c("Normal", "IM"),   
                   c("Normal", "GC"),  
                   c("IM", "GC"))  

plot_data_long$lesion=factor(plot_data_long$lesion, level=c( "Normal" , "IM" , "GC"))
plot_data_long$TF <- factor(plot_data_long$TF, levels = c("NRF1","SP2","IRF8","ZNF416","GMEB2","SPI1"))

#boxplot
setwd("/data1/gy/ATAC_for_review/Figure2H/output")
p <- ggplot(plot_data_long, aes(x = lesion, y = TMM)) +  
  geom_boxplot(  
    aes(color = lesion),   
    fill = "white",  
    width = 0.6,  
    outlier.shape = NA  
  ) +  
  geom_point(  
    aes(color = lesion),  
    position = position_jitterdodge(  
      jitter.width = 0.2,  
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
  scale_color_manual(values = c("#4E79A7", "#59A14F", "#E15759")) +  
  labs(  
    x = "Lesion",  
    y = "TF expression (TMM)"  
  ) +  
  facet_wrap(  
    ~ TF,   
    nrow = 2,   
    ncol = 3,  
    scales = "free_y" 
  )
  
p <- p +
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),  
    axis.text.y = element_text(size = 8),  
    legend.position = "none"  
  )  
ggsave(p,filename="Multistage_RNAseq.6Cluster3-specific_TFs.boxplot.pdf",width=7,height=6.2) ##Figure2H
