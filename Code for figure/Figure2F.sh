####################
######Figure2F######
####################
##########################
###### 1.TMM matrix ######
##########################
mkdir -p /data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_binary_file
###############
/Public/gaoyun/software/R-4.2.0/bin/R
#
library(data.table) 
library(limma)  
library(ComplexHeatmap) 
library(dplyr)
# 
TMM=read.delim("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.TMM",row.names=1,header=T)
TMM=as.matrix(TMM[,-c(1:5)])
TMM_trans=as.data.frame(t(TMM))
baseline=read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv")
TMM_trans$cluster <- baseline$Cluster 

################################################
###### 2.mean & se per OCR in each cluster ######
################################################
## row: sample, col: mean&sd per peak
cluster_stats <- TMM_trans %>%  
  group_by(cluster) %>%  
  summarise(across(everything(), list(mean = mean, sd = sd), .names = "{.col}_{.fn}")) 
fwrite(cluster_stats,"/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_binary_file/95sample_IterativeOverlapPeakSet.TMM_mean_sd.intra_cluster.txt",col.names=T,row.names=F,quote=F,sep="\t")
##计算peak数目
peak_num=(ncol(cluster_stats)-1)/2
peak_num
#[1] 109187

#############################
###### 3.binary matrix ######
#############################
binary_matrix <- list() 
##
for (i in 1:peak_num) {  
    #
    cluster_stats_sub = cluster_stats[, c(1, 2 * i, 2 * i + 1)]  
    
    #
    peak_id <- sub("_mean$", "", names(cluster_stats_sub)[2])  
    
    #
    sorted_cluster_stats_sub <- cluster_stats_sub %>%  
        arrange(.[[2]])

    #
    break_point <- NA  
  
    #
    for (k in 1:(nrow(sorted_cluster_stats_sub) - 1)) {  
      #
      max_mean_sd <- max(sorted_cluster_stats_sub[1:k, 2] + sorted_cluster_stats_sub[1:k, 3])  
      
      #
      if (sorted_cluster_stats_sub[k + 1, 2] > max_mean_sd) {  
        break_point <- k + 1  
        break  
      }  
    }  

    #
    if (!is.na(break_point)) {   
        sorted_cluster_stats_sub[, peak_id] = c(rep("0", break_point - 1), rep("1", nrow(sorted_cluster_stats_sub) - break_point + 1))  
        binary_matrix[[i]] = sorted_cluster_stats_sub[, c(1, which(names(sorted_cluster_stats_sub) == peak_id))]  
    }  
}  

######################################
library(purrr)  
#
filtered_binary_matrix <- binary_matrix[!sapply(binary_matrix, is.null)]  
#
sorted_filtered_binary_matrix <- lapply(filtered_binary_matrix, function(tibble) {  
    tibble %>%  
        arrange(cluster) 
})  
# 
cluster_column <- sorted_filtered_binary_matrix[[1]]$cluster  
#
second_columns <- lapply(sorted_filtered_binary_matrix, function(tibble) {  
    tibble[, 2, drop = FALSE]
})  
#
combined_second_columns <- do.call(cbind, second_columns)  
#
final_binary_matrix <- data.frame(cluster = cluster_column, combined_second_columns)
#
colnames(final_binary_matrix)[-1] <- sub("\\.", ":", colnames(final_binary_matrix)[-1]) 
colnames(final_binary_matrix)[-1] <- sub("\\.", "-", colnames(final_binary_matrix)[-1])
##
dim(final_binary_matrix)
#[1]     3 44422
fwrite(final_binary_matrix,"/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_binary_file/95sample_IterativeOverlapPeakSet.TMM.binary_matrix.txt",col.names=T,row.names=F,quote=F,sep="\t")

#######################################
# 
final_binary_matrix[-1] <- apply(final_binary_matrix[-1], 2, as.numeric)  
# 
column_sums <- colSums(final_binary_matrix[, -1], na.rm = TRUE)  
#
columns_to_keep <- names(column_sums[column_sums == 1])  
# 
final_filtered_binary_matrix <- final_binary_matrix[, c("cluster", columns_to_keep)] 
dim(final_filtered_binary_matrix)
#[1]     3 12800
fwrite(final_filtered_binary_matrix,"/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/95sample_IterativeOverlapPeakSet.TMM.binary_matrix.single_cluster_specific.txt",col.names=T,row.names=F,quote=F,sep="\t")

######################
#######4.limma########
######################
####################################final_filtered_binary_matrix
final_filtered_binary_matrix=fread("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_binary_file/95sample_IterativeOverlapPeakSet.TMM.binary_matrix.single_cluster_specific.txt")

#Cluster-specific row
cluster1_row <- final_filtered_binary_matrix[final_filtered_binary_matrix$cluster == "Cluster1", ]   
cluster2_row <- final_filtered_binary_matrix[final_filtered_binary_matrix$cluster == "Cluster2", ]  
cluster3_row <- final_filtered_binary_matrix[final_filtered_binary_matrix$cluster == "Cluster3", ] 
#
cluster1_peaks <- colnames(cluster1_row)[cluster1_row == 1]  
cluster2_peaks <- colnames(cluster2_row)[cluster2_row == 1]
cluster3_peaks <- colnames(cluster3_row)[cluster3_row == 1]

###################################
TMM_log2=read.delim("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.log2TMM",row.names=1,header=T)
TMM_log2_mat=TMM_log2[,-c(1:5)]
##
TMM_log2_mat_norm=t(normalize.quantiles(t(TMM_log2_mat))) 
rownames(TMM_log2_mat_norm)=rownames(TMM_log2_mat)
colnames(TMM_log2_mat_norm)=colnames(TMM_log2_mat)
#保留cluster-specific
cluster1_TMM_log2_mat_norm=TMM_log2_mat[rownames(TMM_log2_mat_norm)%in%cluster1_peaks,]
cluster2_TMM_log2_mat_norm=TMM_log2_mat[rownames(TMM_log2_mat_norm)%in%cluster2_peaks,]
cluster3_TMM_log2_mat_norm=TMM_log2_mat[rownames(TMM_log2_mat_norm)%in%cluster3_peaks,]

##################################
baseline=read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.baseline.csv",row.names=1)
##
identical(rownames(baseline), colnames(TMM_log2_mat_norm)) 
#[1] TRUE
# Cluster1_binary  
baseline$Cluster1_binary <- ifelse(baseline$Cluster == "Cluster1", "Yes", "No")  
# Cluster2_binary  
baseline$Cluster2_binary <- ifelse(baseline$Cluster == "Cluster2", "Yes", "No")  
# Cluster3_binary  
baseline$Cluster3_binary <- ifelse(baseline$Cluster == "Cluster3", "Yes", "No")
#
cluster1_design <- model.matrix(~ 0 + factor(baseline$Cluster1_binary))
cluster2_design <- model.matrix(~ 0 + factor(baseline$Cluster2_binary))
cluster3_design <- model.matrix(~ 0 + factor(baseline$Cluster3_binary))
colnames(cluster1_design) = levels(factor(baseline$Cluster1_binary))
colnames(cluster2_design) = levels(factor(baseline$Cluster2_binary))
colnames(cluster3_design) = levels(factor(baseline$Cluster3_binary))
rownames(cluster1_design) = colnames(cluster1_TMM_log2_mat_norm)
rownames(cluster2_design) = colnames(cluster2_TMM_log2_mat_norm)
rownames(cluster3_design) = colnames(cluster3_TMM_log2_mat_norm)

##################################
dir.create("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_binary_limma")
setwd("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_binary_limma")
###############cluster1
cluster1_fit <- lmFit(cluster1_TMM_log2_mat_norm, cluster1_design)
cluster1_contrast <- makeContrasts('Yes-No', levels = cluster1_design)
cluster1_contrast
cluster1_fit2 <- contrasts.fit(cluster1_fit, cluster1_contrast)
cluster1_fit2 <- eBayes(cluster1_fit2)
qqt(cluster1_fit2$t, df=cluster1_fit2$df.prior+cluster1_fit2$df.residual, pch=16, cex=0.2)
cluster1_results <- topTable(cluster1_fit2, coef = 1, number = Inf, adjust.method = "BH")
fwrite(cbind(peakID=rownames(cluster1_results),cluster1_results),"95sample_IterativeOverlapPeakSet.TMM.cluster1_specific.limma_result.txt",sep="\t",row.names=F,col.names=T,quote=F)
cluster1_results_FDR0.001=subset(cluster1_results,adj.P.Val<0.001)
dim(cluster1_results_FDR0.001)
#[1] 5849    6
fwrite(cbind(peakID=rownames(cluster1_results_FDR0.001),cluster1_results_FDR0.001),"95sample_IterativeOverlapPeakSet.TMM.cluster1_specific.limma_result,FDR0.001.txt",sep="\t",row.names=F,col.names=T,quote=F)
###############cluster2
cluster2_fit <- lmFit(cluster2_TMM_log2_mat_norm, cluster2_design)
cluster2_contrast <- makeContrasts('Yes-No', levels = cluster2_design)
cluster2_contrast
cluster2_fit2 <- contrasts.fit(cluster2_fit, cluster2_contrast)
cluster2_fit2 <- eBayes(cluster2_fit2)
qqt(cluster2_fit2$t, df=cluster2_fit2$df.prior+cluster2_fit2$df.residual, pch=16, cex=0.2)
cluster2_results <- topTable(cluster2_fit2, coef = 1, number = Inf, adjust.method = "BH")
fwrite(cbind(peakID=rownames(cluster2_results),cluster2_results),"95sample_IterativeOverlapPeakSet.TMM.cluster2_specific.limma_result.txt",sep="\t",row.names=F,col.names=T,quote=F)
cluster2_results_FDR0.001=subset(cluster2_results,adj.P.Val<0.001)
dim(cluster2_results_FDR0.001)
#[1] 47  6
fwrite(cbind(peakID=rownames(cluster2_results_FDR0.001),cluster2_results_FDR0.001),"95sample_IterativeOverlapPeakSet.TMM.cluster2_specific.limma_result,FDR0.001.txt",sep="\t",row.names=F,col.names=T,quote=F)
###############cluster3
cluster3_fit <- lmFit(cluster3_TMM_log2_mat_norm, cluster3_design)
cluster3_contrast <- makeContrasts('Yes-No', levels = cluster3_design)
cluster3_contrast
cluster3_fit2 <- contrasts.fit(cluster3_fit, cluster3_contrast)
cluster3_fit2 <- eBayes(cluster3_fit2)
qqt(cluster3_fit2$t, df=cluster3_fit2$df.prior+cluster3_fit2$df.residual, pch=16, cex=0.2)
cluster3_results <- topTable(cluster3_fit2, coef = 1, number = Inf, adjust.method = "BH")
fwrite(cbind(peakID=rownames(cluster3_results),cluster3_results),"95sample_IterativeOverlapPeakSet.TMM.cluster3_specific.limma_result.txt",sep="\t",row.names=F,col.names=T,quote=F)
cluster3_results_FDR0.001=subset(cluster3_results,adj.P.Val<0.001)
dim(cluster3_results_FDR0.001)
#[1] 6878    6
fwrite(cbind(peakID=rownames(cluster3_results_FDR0.001),cluster3_results_FDR0.001),"95sample_IterativeOverlapPeakSet.TMM.cluster3_specific.limma_result,FDR0.001.txt",sep="\t",row.names=F,col.names=T,quote=F)

##################################################################
#######5.distal binarization cluster-specific peak heatmap########
##################################################################
#########################################################baseline
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
######
cluster_matrix=read.csv("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_consensus_heatmap/95sample_TMM_log2.exp.ConsensusClusterPlus.kmeans_3cluster.matrix.csv",row.names=1)
# 
cluster_row_names <- rownames(cluster_matrix)  
#
baseline_final_sorted <- baseline_final[match(cluster_row_names, rownames(baseline_final)), ]  

#########################################################specific_peak
cluster1_specific=fread("95sample_IterativeOverlapPeakSet.TMM.cluster1_specific.limma_result,FDR0.001.txt")
cluster1_specific$cluster="Cluster1"
cluster2_specific=fread("95sample_IterativeOverlapPeakSet.TMM.cluster2_specific.limma_result,FDR0.001.txt")
cluster2_specific$cluster="Cluster2"
cluster3_specific=fread("95sample_IterativeOverlapPeakSet.TMM.cluster3_specific.limma_result,FDR0.001.txt")
cluster3_specific$cluster="Cluster3"
sorted_final_filtered_binary_df =rbind(cluster1_specific,cluster2_specific,cluster3_specific)
#获取排序后的peakid  
sorted_peakid <- sorted_final_filtered_binary_df$peakID  

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
Cluster_specific_peak_colors <- c("Cluster1" = "#0570B0", "Cluster2" = "#02818A", "Cluster3" = "#D7301F") 

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
    df = sorted_final_filtered_binary_df[, "cluster", drop = FALSE], 
    col = list(  
        cluster = Cluster_specific_peak_colors
    ),  
    which = "row",  
    annotation_label= "Cluster-specific OCRs",
    annotation_legend_param = list(   
        cluster = list(title = "Cluster-specific OCRs")  
    ),  
    annotation_name_side = "bottom",  
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
    border = TRUE,  
    show_annotation_name = TRUE  
)  

## heatmap
heatdata=TMM_log2_mat
zscore_heatdata <- t(apply(heatdata, 1, scale))
colnames(zscore_heatdata)=colnames(heatdata)
#
sorted_zscore_heatdata <- zscore_heatdata[,cluster_row_names]
#
sorted_zscore_heatdata_final <- sorted_zscore_heatdata [sorted_peakid,]

#
library(circlize)
color_mapping <- colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C")) 
##
p = Heatmap(as.matrix(sorted_zscore_heatdata_final),
            cluster_rows = F ,
            cluster_columns = F ,
            show_column_dend = F, 
            show_row_dend = F, 
            show_column_names = F,
            show_row_names = F,
            row_title = NULL, 
            column_title = NULL,
            column_split = baseline_final_sorted$Cluster, 
            row_split = sorted_final_filtered_binary_df$cluster,
            column_gap = unit(1, "mm"),
            row_gap = unit(1, "mm"), 
            top_annotation = top_anno,
            left_annotation = left_anno,  
            heatmap_legend_param = list(title = 'ATAC Zscore'),
            col = color_mapping,
            column_names_gp = gpar(fontsize = 8), 
            border = "black",
            border_gp = gpar(col = "black", lwd = 2)) 
        
######pdf  
setwd("/data1/gy/ATAC_for_review/Figure2F/output")
pdf("95sample_TMM_log2.exp.distal_binarization.heatmap.pdf", width = 15, height = 14)  ##Figure2F
draw(p)  
dev.off()  
