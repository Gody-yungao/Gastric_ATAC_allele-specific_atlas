####################
######Figure1E######
####################
##########################
###### 1.TMM matrix ######
##########################
mkdir -p /data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file
###############
/Public/gaoyun/software/R-4.2.0/bin/R
#
library(data.table) 
library(limma)      
library(ComplexHeatmap) 
library(dplyr)
#
TMM=read.delim("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/final_TMM/pantissue.pantissue_OCR.IterativeOverlapPeakSet.final.TMM",row.names=1,header=T)
TMM_trans=as.data.frame(t(TMM))
library(openxlsx)
baseline=read.xlsx("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/QCsample_metadata/pantissue_allsample_afterQC_metadata.xlsx")
#
order_check <- all(baseline$sample_new == rownames(TMM_trans))
order_check
#[1] TRUE
#
TMM_trans$cluster <- baseline$tissue

##################################################
###### 2.mean & se per OCR in each cluster  ######
##################################################
##
cluster_stats <- TMM_trans %>%  
  group_by(cluster) %>%  
  summarise(across(everything(), list(mean = mean, sd = sd), .names = "{.col}_{.fn}")) 
TMM_trans_dt <- as.data.table(TMM_trans)
cluster_stats <- TMM_trans_dt[, lapply(.SD, function(x) .(mean = mean(x), sd = sd(x))), by = cluster]
#
cluster_stats_mean <- cluster_stats[seq(1, nrow(cluster_stats), by = 2), ]
#
cluster_stats_sd <- cluster_stats[seq(2, nrow(cluster_stats), by = 2), ]
# 
colnames(cluster_stats_mean)[-1] <- paste0(colnames(cluster_stats_mean)[-1], "_mean")
# 
colnames(cluster_stats_sd)[-1] <- paste0(colnames(cluster_stats_sd)[-1], "_sd")
# 
merged_stats <- cbind(cluster_stats_mean, cluster_stats_sd[ , -1])
#
#
new_colnames <- colnames(merged_stats)[-1]
#
sorted_colnames <- sort(new_colnames)
#
setcolorder(merged_stats, c("cluster", sorted_colnames))
fwrite(merged_stats,"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM_mean_sd.intra_cluster.txt",col.names=T,row.names=F,quote=F,sep="\t")
##
peak_num=(ncol(merged_stats)-1)/2
peak_num
#[1] 375455

#############################
###### 3.binary matrix ######
#############################
setDT(merged_stats)

binary_matrix <- list() 
##
for (i in 1:peak_num) {  
    #
    cluster_stats_sub <- merged_stats[, .(cluster, 
                                          mean_value = unlist(merged_stats[[2 * i]]), 
                                          sd_value = unlist(merged_stats[[2 * i + 1]]))]
    
    
    #
    peak_id <- sub("_mean$", "", names(merged_stats)[2 * i])  
    
    #
    sorted_cluster_stats_sub <- cluster_stats_sub[order(cluster_stats_sub[[2]])]

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
        sorted_cluster_stats_sub[, (peak_id) := c(rep("0", break_point - 1), 
                                                   rep("1", nrow(sorted_cluster_stats_sub) - break_point + 1))]  
        binary_matrix[[i]] <- sorted_cluster_stats_sub[, .(cluster, get(peak_id))] 
        colnames(binary_matrix[[i]])[2] = peak_id
    }  
}

######################################binary matrix
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
#[1]     19 303174
fwrite(final_binary_matrix,"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.binary_matrix.txt",col.names=T,row.names=F,quote=F,sep="\t")

#######################################
#
final_binary_matrix[-1] <- apply(final_binary_matrix[-1], 2, as.numeric)  
#
column_sums <- colSums(final_binary_matrix[, -1], na.rm = TRUE)  
#
columns_to_keep <- names(column_sums[column_sums <= 4])  
#
final_filtered_binary_matrix <- final_binary_matrix[, c("cluster", columns_to_keep)] 
dim(final_filtered_binary_matrix)
#[1]     19 189863
fwrite(final_filtered_binary_matrix,"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.binary_matrix.below4_tissue_specific.txt",col.names=T,row.names=F,quote=F,sep="\t")

######################
#######4.limma########
######################
####################################
final_filtered_binary_matrix=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.binary_matrix.below4_tissue_specific.txt")
final_filtered_binary_matrix[1, cluster := gsub("Retina", "retina", cluster)]
#
final_filtered_binary_matrix <- final_filtered_binary_matrix[order(cluster)]

###################################
TMM_log2=read.delim("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/final_TMM/pantissue.pantissue_OCR.IterativeOverlapPeakSet.final.log2TMM",row.names=1,header=T)
TMM_log2_mat=TMM_log2
##
library(preprocessCore)
TMM_log2_mat_norm=t(normalize.quantiles(t(TMM_log2_mat))) 
rownames(TMM_log2_mat_norm)=rownames(TMM_log2_mat)
colnames(TMM_log2_mat_norm)=colnames(TMM_log2_mat)
#
specific_TMM_log2_mat_norm=TMM_log2_mat[rownames(TMM_log2_mat_norm) %in% colnames(final_filtered_binary_matrix),]
##################################
library(openxlsx)
baseline=read.xlsx("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/QCsample_metadata/pantissue_allsample_afterQC_metadata.xlsx")
baseline$tissue <- gsub("Retina", "retina", baseline$tissue)
rownames(baseline)=baseline$sample_new
##
identical(rownames(baseline), colnames(TMM_log2_mat_norm)) 
#[1] TRUE

#
baseline_binary=merge(baseline,final_filtered_binary_matrix,by.x="tissue",by.y="cluster")

##############################
library(limma)
#
TMM_log2_mat_norm <- TMM_log2_mat_norm[, baseline_binary$sample_new]
#
peak_list <- colnames(baseline_binary)[5:ncol(baseline_binary)]
length(peak_list)
#[1] 189862

#
results_list <- list()
for (peak_id in peak_list) {
    #
    group <- as.factor(baseline_binary[[peak_id]])  
    names(group) <- baseline_binary$sample_new 

    #
    expr_vec <- as.numeric(TMM_log2_mat_norm[peak_id,])
    
    #
    design <- model.matrix(~ group)
    fit <- lmFit(expr_vec, design)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=2, adjust.method="fdr", number=1) # coef=2为group1 vs group0
    
    #
    results_list[[peak_id]] <- res
}

#
all_results <- do.call(rbind, results_list)
dim(all_results)
#[1] 189862      6
all_results$adj.P.Val=p.adjust(all_results$P.Value,method="BH")
fwrite(cbind(peakID=rownames(all_results),all_results),"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.below4_tissue_specific.limma_result.txt",sep="\t",row.names=F,col.names=T,quote=F)
all_results_0.001=subset(all_results,adj.P.Val<0.001)
dim(all_results_0.001)
#[1] 189279      6
fwrite(cbind(peakID=rownames(all_results_0.001),all_results_0.001),"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.below4_tissue_specific.limma_result_FDR0.001.txt",sep="\t",row.names=F,col.names=T,quote=F)

#
selected_peaks <- intersect(colnames(final_filtered_binary_matrix)[-1], rownames(all_results_0.001))
#
final_filtered_binary_matrix_limma <- final_filtered_binary_matrix[, c("cluster", selected_peaks), with = FALSE]
dim(final_filtered_binary_matrix_limma)
#[1]     19 189280
fwrite(final_filtered_binary_matrix_limma,"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.below4_tissue_specific.limma_result_FDR0.001.binary_matrix.txt",sep="\t",row.names=F,col.names=T,quote=F)

#
peak_count <- rowSums(final_filtered_binary_matrix_limma[,-1] == 1)
result <- data.frame(cluster = final_filtered_binary_matrix_limma$cluster, peak_count = peak_count)
#          cluster peak_count
#1   adrenal-gland      23001
#2           brain      57009
#3           colon      21731
#4       esophagus       7290
#5  fallopian-tube       9507
#6             fat       3078
#7           heart      13197
#8          kidney      17624
#9           liver      25091
#10           lung       8164
#11         muscle      23411
#12          nerve       8463
#13          ovary       8491
#14       pancreas      20376
#15         retina      38578
#16           skin       6189
#17         spleen       9159
#18        stomach       9598
#19  thyroid-gland      19342

###############stomach
row_stomach <- which(final_filtered_binary_matrix_limma$cluster == "stomach")
#
cols_stomach_1 <- which(as.numeric(final_filtered_binary_matrix_limma[row_stomach, -1, with=FALSE]) == 1) + 1  # +1补上cluster列
#
subset_dt <- final_filtered_binary_matrix_limma[, c(1, cols_stomach_1), with=FALSE]
fwrite(subset_dt,"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.below4_tissue_specific.stomach.limma_result_FDR0.001.binary_matrix.txt",sep="\t",row.names=F,col.names=T,quote=F)


##################################################################
#######5.distal binarization cluster-specific peak heatmap########
##################################################################
####
baseline=read.xlsx("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/QCsample_metadata/pantissue_allsample_afterQC_metadata.xlsx")
baseline$tissue <- gsub("Retina", "retina", baseline$tissue)
rownames(baseline)=baseline$sample_new

####
TMM_log2_mat_df <- as.data.frame(t(TMM_log2_mat))
##
identical(rownames(baseline), rownames(TMM_log2_mat_df)) 
#[1] TRUE
#
library(tidyr) 
TMM_log2_mat_df <- TMM_log2_mat_df %>%
  mutate(tissue = baseline$tissue) %>%
  gather(key = "peak", value = "expression", -tissue)
#
mean_expression <- TMM_log2_mat_df %>%
  #
  group_by(peak, tissue) %>%
  #
  summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = 'drop') %>%
  #
  pivot_wider(names_from = tissue, values_from = mean_expression)
mean_expression=as.data.table(mean_expression)
dim(mean_expression)
#[1] 375455     20

####specific_peak
tissue_specific_peak=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.below4_tissue_specific.limma_result_FDR0.001.txt")
mean_expression_tissue_specific=subset(mean_expression,peak %in% tissue_specific_peak$peakID)
dim(mean_expression_tissue_specific)
#[1] 189279     20

####
df <- as.data.frame(mean_expression_tissue_specific)
#
max_tissue <- colnames(df)[-1][apply(df[,-1], 1, which.max)]
#
df$max_tissue <- max_tissue
specific_peak_max_tissue=df[, c("peak", "max_tissue")]
specific_peak_max_tissue <- specific_peak_max_tissue[order(specific_peak_max_tissue$max_tissue), ]

#
tissue_order <- levels(factor(specific_peak_max_tissue$max_tissue, levels=unique(specific_peak_max_tissue$max_tissue)))
#
tissue_order_new <- c(tissue_order[1:(length(tissue_order) - 2)],
                      tissue_order[length(tissue_order)],
                      tissue_order[length(tissue_order) - 1])
#
specific_peak_max_tissue$max_tissue <- factor(
  specific_peak_max_tissue$max_tissue,
  levels = tissue_order_new
)
#
specific_peak_max_tissue <- specific_peak_max_tissue[order(specific_peak_max_tissue$max_tissue), ]
#
sorted_peakid <- specific_peak_max_tissue$peak

#
specific_peak_stomach=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.below4_tissue_specific.stomach.limma_result_FDR0.001.binary_matrix.txt")
##
df <- as.data.frame(specific_peak_stomach)
df$total_count <- rowSums(df[ , -1, drop = FALSE])

##
specific_peakid_stomach=colnames(specific_peak_stomach)[2:length(specific_peak_stomach)]
sorted_peakid_new=c(
  sorted_peakid[!sorted_peakid %in% specific_peakid_stomach],
  sorted_peakid[ sorted_peakid %in% specific_peakid_stomach]
)

########################################################heatmap
heatdata=as.data.frame(mean_expression_tissue_specific)
rownames(heatdata)=heatdata$peak
heatdata=heatdata[,-1]
zscore_heatdata <- t(apply(heatdata, 1, scale))
colnames(zscore_heatdata)=colnames(heatdata)
#
sorted_tissue=as.character(unique(specific_peak_max_tissue$max_tissue))
sorted_tissue
# [1] "adrenal-gland"  "brain"          "colon"          "esophagus"     
# [5] "fallopian-tube" "fat"            "heart"          "kidney"        
# [9] "liver"          "lung"           "muscle"         "nerve"         
#[13] "ovary"          "pancreas"       "retina"         "skin"          
#[17] "spleen"         "thyroid-gland"  "stomach"  

##
sorted_zscore_heatdata <- zscore_heatdata[,sorted_tissue]
sorted_zscore_heatdata <- sorted_zscore_heatdata[sorted_peakid_new,]

#
library(circlize)
color_mapping <- colorRamp2(c(-3, 0, 3), c("#377EB8", "white", "#E41A1C")) 
#
mat <- t(as.matrix(sorted_zscore_heatdata))

##
setwd("/data1/gy/ATAC_for_review/Figure1E/output")
p = Heatmap(mat,
            cluster_rows = F ,
            cluster_columns = F ,
            show_column_dend = F,
            show_row_dend = F,
            show_column_names = F,
            show_row_names = T, 
            row_title = NULL, 
            column_title = NULL, 
            #column_split = baseline_final_sorted$Cluster,
            #row_split = sorted_final_filtered_binary_df$cluster,
            column_gap = unit(1, "mm"),
            row_gap = unit(1, "mm"),
            heatmap_legend_param = list(
              title = 'ATAC Zscore',
              at = c(-3,-1.5,0,1.5,3),
              labels = c("-3","-1.5","0","1.5","3")
            ),
            col = color_mapping,
            column_names_gp = gpar(fontsize = 8),
            border = "black",
            use_raster=T,
            border_gp = gpar(col = "black", lwd = 2))
        
######
pdf("pantissue_TMM_log2.exp.distal_binarization.heatmap.pdf", width = 10, height = 8)  ##Figure1E
draw(p)  
dev.off()  
