###########################################
################Figure2B###################
###########################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(ConsensusClusterPlus)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(dplyr)  
#####baseline
baseline=openxlsx::read.xlsx("/data1/gy/ATACseq_RWAS/ATACseq/baseline/95sample_baseline.final.foranalysis.xlsx")
####Gastric_lesion
baseline$Gastric_lesion_status <- recode(
  baseline$Gastric_lesion_status,
  "Atrophic gastritis" = "AG",
  "Intestinal metaplasia" = "IM",
  "Superficial gastritis" = "SG"
)
###
rownames(baseline)=baseline$ID
baseline=baseline[,-c(1)]

#####TMM log2 matrix
TMM_log2=as.data.frame(fread("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.log2TMM"))
rownames(TMM_log2)=TMM_log2$Geneid
TMM_log2_mat=TMM_log2[,-c(1:6)]

#####ConsensusCluster kmeans cluster
maxK <- 5 
cluster_result_km <- ConsensusClusterPlus(as.matrix(TMM_log2_mat),
                                        maxK = maxK,
                                        reps = 1000, 
                                        pItem = 0.8, 
                                        pFeature = 1,
                                        clusterAlg = "km", 
                                        distance="euclidean", 
                                        title="/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap", 
                                        innerLinkage="complete",
                                        plot="pdf",
                                        seed = 123)

#cluster=3
i=3
num_clusters <- i
    
# annCol
annCol <- data.frame(results = paste0("Cluster", cluster_result_km[[i]][['consensusClass']]), row.names = colnames(TMM_log2_mat))  

# baseline添加cluster并导出  
baseline_final = cbind(baseline, annCol)  
colnames(baseline_final)[ncol(baseline_final)] = "Cluster"
    
#
baseline_final$Cluster <- factor(baseline_final$Cluster, levels = paste0("Cluster", 1:num_clusters))  
fwrite(cbind("ID" = rownames(baseline_final), baseline_final), paste0("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_", i, "cluster.baseline.csv")) 
    
#
heatdata <- cluster_result_km[[i]][["consensusMatrix"]]  
dimnames(heatdata) <- list(colnames(TMM_log2_mat), colnames(TMM_log2_mat))  
    
#
sample_order <- cluster_result_km[[i]]$consensusTree$order  
sample_ids <- colnames(TMM_log2_mat)[sample_order]
annCol=annCol[sample_ids,"results", drop = FALSE]
annCol=annCol[c(1:32,57:95,33:56),, drop = FALSE]

# 
heatdata <- heatdata[rownames(annCol), rownames(annCol)]
    
#
fwrite(heatdata, paste0("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_", i, "cluster.matrix.csv"), row.names = TRUE)  


######################################## diff test
library(dplyr)  
baseline_final=fread("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv")  
#
col_names <- colnames(baseline_final)[c(3:13,15:20)] ## rm "Fried foods intake"

#
ordered_vars <- c("age_group", "Gastric_lesion_status")   

# 
baseline_final$age_group <- factor(baseline_final$age_group,  
                                  ordered = TRUE,  
                                  levels = c("40-49", "50-59", "60-69"))

baseline_final$Gastric_lesion_status <- factor(baseline_final$Gastric_lesion_status,  
                                              ordered = TRUE,  
                                              levels = c("SG", "AG", "IM"))  

# 
results <- lapply(col_names, function(column_name) {  
  if (column_name %in% ordered_vars) {  
    # Kruskal-Wallis 
    kw_result <- kruskal.test(as.formula(paste(column_name, "~ Cluster")),  
                             data = baseline_final)  
    return(c(P_value = kw_result$p.value))  
    
  } else {  
    #2*2 table
    table_data <- table(baseline_final[[column_name]],   
                       baseline_final$Cluster)  
    # Fisher test  
    fisher_result <- fisher.test(table_data)  
    return(c(P_value = fisher_result$p.value))  
  }  
})  
#
results_df <- data.frame(Variable = col_names, P_value = unlist(results)) 
setwd("/data1/gy/ATAC_for_review/Figure2B/output")
fwrite(results_df,"95sample_TMM_log2.exp.ConsensusClusterPlus.Kruskal-Wallis.fishertest.among_3kmeans_cluster.csv")

######################################################## heatmap 
#
annCol=read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.matrix.csv",row.names=1)
baseline_final = read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv", row.names = 1)  
baseline_final = baseline_final[rownames(annCol),,drop=F]
baseline_final=baseline_final[,c(3,2,5:9,20)]
baseline_final$age_group <- factor(baseline_final$age_group, levels = c("40-49", "50-59", "60-69"))  
baseline_final$sex <- factor(baseline_final$sex, levels = c("Male", "Female"))  
baseline_final$Gastric_lesion_status <- factor(baseline_final$Gastric_lesion_status, levels = c("SG", "AG", "IM"))  
baseline_final$HP_C14 <- factor(baseline_final$HP_C14, levels = c("Positive", "Negative"))  
baseline_final$smoke_status <- factor(baseline_final$smoke_status, levels = c("Smokers", "Nonsmokers", "NA"))  
baseline_final$drink_status <- factor(baseline_final$drink_status, levels = c("Drinkers", "Nondrinkers", "NA"))  
baseline_final$Tea_consumption <- factor(baseline_final$Tea_consumption, levels = c("Tea-drinkers", "Nondrinkers", "NA"))
baseline_final$Cluster<- factor(baseline_final$Cluster, levels = c("Cluster1", "Cluster2", "Cluster3"))  
    
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
    df = baseline_final,  
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

left_anno <- HeatmapAnnotation(  
    df = baseline_final[, "Cluster", drop = FALSE],
    col = list(  
        Cluster = Cluster_colors  
    ),  
    which = "row",  
    annotation_label= " ",
    annotation_legend_param = list(   
        Cluster = list(title = "Cluster") 
    ),  
    annotation_name_side = "bottom",  
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
    border = TRUE,  
    show_annotation_name = TRUE  
)  

##  
heatdata=read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.matrix.csv",row.names=1)
p = Heatmap(as.matrix(heatdata),
            cluster_rows = F ,
            cluster_columns = F ,
            show_column_dend = F,
            show_row_dend = F,
            show_column_names = F,
            show_row_names = F,
            row_title = NULL, 
            column_title = NULL, 
            row_gap = unit(0, "mm"), 
            top_annotation = top_anno,  
            left_annotation = left_anno,
            heatmap_legend_param = list(title = ''), 
            col = colorRampPalette(c("white", "steelblue"))(100), 
            column_names_gp = gpar(fontsize = 8), 
            border = "black",
            border_gp = gpar(col = "black", lwd = 2))
        
###### pdf  
setwd("/data1/gy/ATAC_for_review/Figure2B/output")
pdf("95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.heatmap.pdf", width = 15, height = 13)  ##Figure2B
draw(p)  
dev.off()  
q()
