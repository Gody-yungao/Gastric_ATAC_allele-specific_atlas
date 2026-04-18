###################################################
##########OCR accessibility calculation############
###################################################
####saf file
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
peakset=fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed",header=F)
peakset$V2=peakset$V2+1
peakset$ID=paste0(peakset$V1,":",peakset$V2,"-",peakset$V3)
peakset$Strand="+"
peakset=peakset[,c(6,1,2,3,7)]
colnames(peakset)=c("#ID","Chr","Start","End","Strand")
peakset$'#ID'=paste0(peakset$Chr,":",peakset$Start,"-",peakset$End)
write.table(peakset,"/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/readcount/95sample_IterativeOverlapPeakSet.saf",col.names=T,row.names=F,quote=F,sep="\t")
q()

#######count Tn5 insertion
#####https://www.biostars.org/p/320394/#320408
conda activate /Public/gaoyun/miniconda3/envs/total_RNAseq
saf=/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/readcount/95sample_IterativeOverlapPeakSet.saf
featureCounts \
-a $saf \
-F SAF \
--read2pos 5 \
-p \
-T 15 \
-o /data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/readcount/95sample_IterativeOverlapPeakSet.readcount \
/data1/gy/ATACseq_RWAS/ATACseq/bam/*.last.shift.sort.bam

###
awk 'BEGIN { OFS = "\t" }
     NR == 1 { if ($0 ~ /^#/) { next } }
     NR == 2 { for (i = 1; i <= NF; ++i) { sub(".*/", "", $i); sub("\\.last\\.shift\\.sort\\.bam", "", $i) } }
     { print }' /data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/readcount/95sample_IterativeOverlapPeakSet.readcount \
     > /data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/readcount/95sample_IterativeOverlapPeakSet.final.readcount

#########TMM normalization (for downstream analysis)
mkdir /data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
###
counts <- read.delim("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/readcount/95sample_IterativeOverlapPeakSet.final.readcount",header = T)
library(edgeR)
y = DGEList(counts=counts[,7:101],genes=counts$Geneid)
y = calcNormFactors(y, method = "TMM")
y$samples$lib.size = colSums(y$counts)
results1 = (t( t(y$counts) / y$samples$lib.size / y$samples$norm.factors)) * 1000000
results2 = edgeR::cpm(y,log = F) ##results1=results2
write.table(cbind(counts[,1:6],results1),"/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.TMM",row.names = F,col.names = T,quote=F,sep = "\t")
write.table(cbind(counts[,1:6],log2(results1+1)),"/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.log2TMM",row.names = F,col.names = T,quote=F,sep = "\t")

#####TMM inverse normal transformation (for caQTL mapping)
##
inverse_qnorm <- function(row) {
  qnorm_value <- qnorm((rank(row, ties.method = "first") - 0.5) / length(row))
  return(qnorm_value)
}
##
results1_new <- t(apply(results1, 1, inverse_qnorm))
write.table(cbind(counts[,1:6],results1_new),"/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.TMM_qnorm",row.names = F,col.names = T,quote=F,sep = "\t")
q()
