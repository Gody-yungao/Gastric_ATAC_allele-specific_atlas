#####################################################
################Figure2D&FigureS5D###################
#####################################################
##############################
######1.GSVA calculation######
##############################
conda activate GSVA
R
library(data.table)
###################
TMM=fread("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.log2TMM")
TMM=TMM[,-c(2:6)]
anno=fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.UCSC_hg19_known_gene.anno_by_chipseeker.with_symbol.txt")
anno=anno[,c(6,14,16)]
TMM=merge(anno,TMM,by.x="V4",by.y="Geneid")
colnames(TMM)[1]="PeakID"

##mean value for the same gene anno 
##aggregate
TMM <- aggregate(TMM[,-c(1,2,3)], list(TMM$geneId), FUN=mean)
dim(TMM)
#[1] 18108    96
rownames(TMM) = TMM$Group.1
TMM=as.matrix(TMM[,-1])

#####################cluster
baseline=read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv")

#####################GSVA
library(GSEABase)
library(GSVA)
library(clusterProfiler)
##################################################HALLMARK
Hallmark <- read.gmt("/data1/gy/public/GSEA_geneset_new/h.all.v2024.1.Hs.entrez.gmt")
Hallmark_list = split(Hallmark$gene, Hallmark$term)
GSVA_Hallmark <- gsva(TMM, Hallmark_list, kcdf="Gaussian",method = "gsva")
dir.create("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/GSVA/Hallmark")
fwrite(cbind("pathway"=rownames(GSVA_Hallmark),GSVA_Hallmark),"/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/GSVA/Hallmark/95sample_TMM_log2.Hallmark.GSVA_score.csv",col.names=T,row.names=F)
q()

###########################################################
#######2.GSVA cluster diff analysis（KW & wilcoxon）#######
###########################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(dplyr)
##
analyze_pathway <- function(pathway, data) {  
  df <- data[, c("Cluster", pathway)]  
  colnames(df)[2] <- "score"  
  df <- df[!is.na(df$score), ]  
  
  # Kruskal-Wallis
  kw_test <- kruskal.test(score ~ Cluster, data = df)  
  
  # 
  res <- list(  
    pathway = pathway,  
    kw_p = kw_test$p.value,  
    p_3v2 = NA, dir_3v2 = NA,  
    p_3v1 = NA, dir_3v1 = NA,  
    p_2v1 = NA, dir_2v1 = NA  
  )  
  
  # 
  comparisons <- list(  
    c("Cluster3", "Cluster2"),  
    c("Cluster3", "Cluster1"),  
    c("Cluster2", "Cluster1")  
  )  
  
  for (pair in comparisons) {  
    suffix <- paste0(substr(pair[1], 8, 8), "v", substr(pair[2], 8, 8))  
    p_col <- paste0("p_", suffix)  
    dir_col <- paste0("dir_", suffix)  
    
    if (all(pair %in% levels(df$Cluster))) {  
      sub_df <- df[df$Cluster %in% pair, ]  
      
      if (nrow(sub_df) >= 2) {  
        # Wilcoxon 
        p_value <- wilcox.test(score ~ Cluster, data = sub_df)$p.value  
        
        #
        medians <- aggregate(score ~ Cluster, data = sub_df,   
                           FUN = median, na.rm = TRUE)  
        medians <- setNames(medians$score, medians$Cluster)  
        
        #
        direction <- if (p_value >= 0.05) {  
          "nosig"  
        } else if (medians[pair[1]] > medians[pair[2]]) {  
          "Up"  
        } else {  
          "Down"  
        }  
        
        #
        res[[p_col]] <- p_value  
        res[[dir_col]] <- direction  
      }  
    }  
  }  
  return(res)  
}  

##
baseline_final = read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv")  

############################HALLMARK
GSVA_Hallmark = read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/GSVA/Hallmark/95sample_TMM_log2.Hallmark.GSVA_score.csv",row.names=1)
GSVA_Hallmark_trans=as.data.frame(t(GSVA_Hallmark))
GSVA_Hallmark_trans$ID=rownames(GSVA_Hallmark_trans)
#
merged_data <- baseline_final %>%   
  left_join(GSVA_Hallmark_trans, by = "ID") %>%   
  mutate(Cluster = factor(Cluster))
#
pathway_names <- setdiff(colnames(merged_data),   
                        c(names(baseline_final), "ID", "Cluster"))  
results <- lapply(pathway_names, function(p) {  
  analyze_pathway(p, merged_data)  
})   
#
result_df <- do.call(rbind, lapply(results, function(x) {  
  data.frame(  
    pathway = x$pathway,  
    Kruskal_Wallis_p = x$kw_p,  
    Cluster3_vs_Cluster2_p = x$p_3v2,  
    Cluster3_vs_Cluster2_sig = x$dir_3v2,  
    Cluster3_vs_Cluster1_p = x$p_3v1,  
    Cluster3_vs_Cluster1_sig = x$dir_3v1,  
    Cluster2_vs_Cluster1_p = x$p_2v1,  
    Cluster2_vs_Cluster1_sig = x$dir_2v1  
  )  
}))  

# BH
result_df$Kruskal_Wallis_fdr <- p.adjust(result_df$Kruskal_Wallis_p, method = "BH")
result_df$Cluster3_vs_Cluster2_fdr <- p.adjust(result_df$Cluster3_vs_Cluster2_p, method = "BH")
result_df$Cluster3_vs_Cluster1_fdr <- p.adjust(result_df$Cluster3_vs_Cluster1_p, method = "BH")
result_df$Cluster2_vs_Cluster1_fdr <- p.adjust(result_df$Cluster2_vs_Cluster1_p, method = "BH")

#
result_df <- result_df %>%
  mutate(
    Cluster3_vs_Cluster2_sig_fdr = case_when(
      is.na(Cluster3_vs_Cluster2_fdr) ~ NA_character_,
      Cluster3_vs_Cluster2_fdr >= 0.05 ~ "nosig",
      Cluster3_vs_Cluster2_sig == "Up" ~ "Up",
      Cluster3_vs_Cluster2_sig == "Down" ~ "Down"
    ),
    Cluster3_vs_Cluster1_sig_fdr = case_when(
      is.na(Cluster3_vs_Cluster1_fdr) ~ NA_character_,
      Cluster3_vs_Cluster1_fdr >= 0.05 ~ "nosig",
      Cluster3_vs_Cluster1_sig == "Up" ~ "Up",
      Cluster3_vs_Cluster1_sig == "Down" ~ "Down"
    ),
    Cluster2_vs_Cluster1_sig_fdr = case_when(
      is.na(Cluster2_vs_Cluster1_fdr) ~ NA_character_,
      Cluster2_vs_Cluster1_fdr >= 0.05 ~ "nosig",
      Cluster2_vs_Cluster1_sig == "Up" ~ "Up",
      Cluster2_vs_Cluster1_sig == "Down" ~ "Down"
    )
  )

#
result_df_final <- result_df %>%
  select(
    pathway,
    Kruskal_Wallis_p, Kruskal_Wallis_fdr,
    Cluster3_vs_Cluster2_p, Cluster3_vs_Cluster2_fdr, Cluster3_vs_Cluster2_sig_fdr,
    Cluster3_vs_Cluster1_p, Cluster3_vs_Cluster1_fdr, Cluster3_vs_Cluster1_sig_fdr,
    Cluster2_vs_Cluster1_p, Cluster2_vs_Cluster1_fdr, Cluster2_vs_Cluster1_sig_fdr
  )

# 
write.csv(result_df_final, "/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/GSVA/Hallmark/95sample_TMM_log2.Hallmark.GSVA_score.diff_result.csv",row.names=F,quote=F)  

############################
#######3.GSVA heatmap#######FigureS5D
############################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(ComplexHeatmap)
library(dplyr)
setwd("/data1/gy/ATAC_for_review/Figure2D&FigureS5D/output")
##
baseline_final = read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv", row.names = 1)  
baseline_final=baseline_final[,c(3,2,5:9,20)]
baseline_final$age_group <- factor(baseline_final$age_group, levels = c("40-49", "50-59", "60-69"))  
baseline_final$sex <- factor(baseline_final$sex, levels = c("Male", "Female"))  
baseline_final$Gastric_lesion_status <- factor(baseline_final$Gastric_lesion_status, levels = c("SG", "AG", "IM"))  
baseline_final$HP_C14 <- factor(baseline_final$HP_C14, levels = c("Positive", "Negative"))  
baseline_final$smoke_status <- factor(baseline_final$smoke_status, levels = c("Smokers", "Nonsmokers", "NA"))  
baseline_final$drink_status <- factor(baseline_final$drink_status, levels = c("Drinkers", "Nondrinkers", "NA"))  
baseline_final$Tea_consumption <- factor(baseline_final$Tea_consumption, levels = c("Tea-drinkers", "Nondrinkers", "NA"))
baseline_final$Cluster<- factor(baseline_final$Cluster, levels = c("Cluster1", "Cluster2", "Cluster3"))  
######
cluster_matrix=read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.matrix.csv",row.names=1)
# 
cluster_row_names <- rownames(cluster_matrix)  
# 
baseline_final_sorted <- baseline_final[match(cluster_row_names, rownames(baseline_final)), ]  

#
age_colors <- c("40-49" = "lightblue", "50-59" = "lightgreen", "60-69" = "lightyellow")  
sex_colors <- c("Male" = "skyblue", "Female" = "pink")  
smoke_colors <- c(
  "Smokers" = "#D73027",  
  "Nonsmokers" = "#1A9850",
  "NA" = "#CCCCCC" 
)
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
Cluster_colors <- c("Cluster1" = "#045A8D", "Cluster2" = "#016C59", "Cluster3" = "#B30000")  
#
top_anno <- HeatmapAnnotation(  
    df = baseline_final_sorted,  
    col = list(  
        age_group = age_colors,  
        sex = sex_colors,  
        HP_C14 = hp_colors,  
        Gastric_lesion_status = Gastric_lesion_colors,  
        smoke_status = smoke_colors,  
        drink_status = drink_colors,  
        Tea_consumption = Tea_consumption_colors,  
        Cluster = Cluster_colors  
    ),  
    which="column",
    annotation_label=c("Age", "Sex", "H.pylori Infection", "Gastric lesion", "Smoking", "Alcohol drinking", "Tea consumption", "Cluster"),
    annotation_legend_param = list(  
        age_group = list(title = "Age Group"),  
        sex = list(title = "Sex"),  
        HP_C14 = list(title = "H.pylori Infection"),  
        Gastric_lesion_status = list(title = "Gastric Lesion"),  
        smoke_status = list(title = "Smoking"),  
        drink_status = list(title = "Alcohol drinking"),  
        Tea_consumption = list(title = "Tea consumption"),  
        Cluster = list(title = "Cluster")  
    ),  
    annotation_name_side = "left",  
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
    border = TRUE,  
    show_annotation_name = TRUE  
)  

##############HALLMARK 
heatdata = read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/GSVA/Hallmark/95sample_TMM_log2.Hallmark.GSVA_score.csv",row.names=1)
heatdata = heatdata[,rownames(baseline_final_sorted)]

##
fh = function(x) hclust(dist(x), method="ward.D")
row_hclust <- fh(heatdata)
row_clusters <- cutree(row_hclust, k = 4)  
 
##################heatmap
#
library(circlize)
color_mapping <- colorRamp2(c(-0.5, 0, 0.5), c("#377EB8", "white", "#E41A1C")) 

p = Heatmap(as.matrix(heatdata),
            cluster_rows = fh , 
            cluster_columns = F ,
            show_column_dend = F,   
            show_row_dend = T, 
            show_column_names = F,
            show_row_names = T, 
            row_title = NULL,
            column_title = NULL,
            row_names_side = "right", 
            row_dend_side = "left",  
            row_split = row_clusters,  
            row_gap = unit(1, "mm"),   
            column_split = baseline_final_sorted$Cluster, 
            column_gap = unit(1, "mm"),
            top_annotation = top_anno,  
            heatmap_legend_param = list(title = 'GSVA score estimated by ATACseq'), 
            col = color_mapping,
            column_names_gp = gpar(fontsize = 8), 
            border = "black", 
            border_gp = gpar(col = "black", lwd = 2))
        
######pdf  
pdf("95sample_TMM_log2.Hallmark.GSVA_score.heatmap.pdf", width = 14, height = 13)  ##FigureS5D
draw(p, heatmap_legend_side = "bottom",annotation_legend_side = "bottom")  
dev.off()  
q()

####################################################
#######4.immune-response pathway GSVA boxplot#######Figure2D
####################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(ggplot2)
library(ggpubr)  
library(tidyr) 
library(dplyr)
setwd("/data1/gy/ATAC_for_review/Figure2D&FigureS5D/output")
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
plot_data <- merged_data[, c("Cluster","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
                              "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB")]    
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

plot_data_long$pathway <- factor(plot_data_long$pathway, levels = c("HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
                              "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB"))

# boxplot
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
    nrow = 2,   
    ncol = 3,  
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
ggsave(p,filename="95sample_TMM_log2.Hallmark.GSVA_score.6immune_response_pathways.boxplot.pdf") ##Figure2D






