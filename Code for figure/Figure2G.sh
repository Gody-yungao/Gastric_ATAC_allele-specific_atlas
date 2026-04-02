###########################################
################Figure2G###################
###########################################
##########################################
#######1.cluster1-specific TF(homer)######
##########################################
##bed file
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
data=fread("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_binary_limma/95sample_IterativeOverlapPeakSet.TMM.cluster1_specific.limma_result,FDR0.001.txt")
#
bed_data <- data.frame(  
  chrom = sub(":.*", "", data$peakID),
  start = as.numeric(sub(".*:(\\d+)-.*", "\\1", data$peakID))-1,
  end = as.numeric(sub(".*-(\\d+)", "\\1", data$peakID)), 
  name = data$peakID 
  #score = data$logFC, 
  #strand = "."
)  
#
bed_data_sorted <- bed_data[order(bed_data$chrom, bed_data$start), ]  
#
write.table(bed_data_sorted, "/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_specific_TF/cluster1/homer_input/95sample_IterativeOverlapPeakSet.TMM.cluster1_specific.limma_result,FDR0.001.bed", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
all_OCR=fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed")
all_OCR$V4=paste0(all_OCR$V1,":",all_OCR$V2+1,"-",all_OCR$V3)
bg_OCR=subset(all_OCR,!(V4 %in% bed_data_sorted$name))
write.table(bg_OCR[,1:4], "/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_specific_TF/cluster1/homer_input/95sample_IterativeOverlapPeakSet.TMM.cluster1_specific.limma_result,FDR0.001.bg.bed", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
q()

###########homer
cd /Public/gaoyun/software/homer 
#perl configureHomer.pl -install hg19
export PATH=$PATH:/Public/gaoyun/software/homer/bin/  
cd /data1/gy/ATAC_for_review/Figure2G/output/cluster1/homer_output
findMotifsGenome.pl \
/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_specific_TF/cluster1/homer_input/95sample_IterativeOverlapPeakSet.TMM.cluster1_specific.limma_result,FDR0.001.bed \
hg19 \
./ \
-bg /data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_specific_TF/cluster1/homer_input/95sample_IterativeOverlapPeakSet.TMM.cluster1_specific.limma_result,FDR0.001.bg.bed \
-size given -p 5 -mask \
-preparsedDir ./ 

############################################
#######2.cluster3-specific TF(homer)########Figure2G
############################################
##bed file
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
data=fread("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_binary_limma/95sample_IterativeOverlapPeakSet.TMM.cluster3_specific.limma_result,FDR0.001.txt")
#
bed_data <- data.frame(  
  chrom = sub(":.*", "", data$peakID), 
  start = as.numeric(sub(".*:(\\d+)-.*", "\\1", data$peakID))-1, 
  end = as.numeric(sub(".*-(\\d+)", "\\1", data$peakID)),
  name = data$peakID 
  #score = data$logFC, 
  #strand = "." 
)  
#
bed_data_sorted <- bed_data[order(bed_data$chrom, bed_data$start), ]  
# 
write.table(bed_data_sorted, "/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_specific_TF/cluster3/homer_input/95sample_IterativeOverlapPeakSet.TMM.cluster3_specific.limma_result,FDR0.001.bed", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)  
all_OCR=fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed")
all_OCR$V4=paste0(all_OCR$V1,":",all_OCR$V2+1,"-",all_OCR$V3)
bg_OCR=subset(all_OCR,!(V4 %in% bed_data_sorted$name))
write.table(bg_OCR[,1:4], "/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_specific_TF/cluster3/homer_input/95sample_IterativeOverlapPeakSet.TMM.cluster3_specific.limma_result,FDR0.001.bg.bed", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
q()

###########homer
cd /Public/gaoyun/software/homer 
#perl configureHomer.pl -install hg19
export PATH=$PATH:/Public/gaoyun/software/homer/bin/  
cd /data1/gy/ATAC_for_review/Figure2G/output/cluster3/homer_output
findMotifsGenome.pl \
/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_specific_TF/cluster3/homer_input/95sample_IterativeOverlapPeakSet.TMM.cluster3_specific.limma_result,FDR0.001.bed \
hg19 \
./ \
-bg /data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM_distal_heatmap/cluster_specific_TF/cluster3/homer_input/95sample_IterativeOverlapPeakSet.TMM.cluster3_specific.limma_result,FDR0.001.bg.bed \
-size given -p 10 -mask \
-preparsedDir ./ \
-len 8,10,12


