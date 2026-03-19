#########################################
################Figure1B#################
#########################################
conda activate /Public/gaoyun/miniconda3/envs/chipseeker
R
setwd("/data1/gy/ATAC_for_review/Figure1B/output")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
peak=read.delim("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed",header=F)
dim(peak)
#[1] 109187      5
peak <- readPeakFile("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed")
peak_Anno <- annotatePeak(peak,
                          level = "gene",
                          tssRegion=c(-1000, 1000),
                          TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
peak_Anno_df = as.data.frame(peak_Anno)
peak_Anno_df$V4=paste0(peak_Anno_df$seqnames,":",peak_Anno_df$start,"-",peak_Anno_df$end)
write.table(peak_Anno_df,"/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.UCSC_hg19_known_gene.anno_by_chipseeker.txt",col.names=T,row.names=F,quote=F,sep="\t")

###genomic anno
library(dplyr)
library(stringr)  
anno_df <- as.data.frame(peak_Anno@anno) %>%  
  #
  mutate(annotation_clean = str_remove(annotation, "\\s*\\(.*$")) %>%  
  #
  mutate(region = case_when(  
    annotation_clean == "Promoter"            ~ "Promoter",  
    annotation_clean == "5' UTR"              ~ "5'UTR",  
    annotation_clean == "3' UTR"              ~ "3'UTR",  
    annotation_clean == "Downstream"          ~ "Distal",
    annotation_clean == "Distal Intergenic"   ~ "Distal",  
    str_detect(annotation_clean, "Exon")      ~ "Exon",  
    str_detect(annotation_clean, "Intron")    ~ "Intron",  
    TRUE                                      ~ "Other"  
  ))  
table(anno_df$region)
#   3'UTR    5'UTR   Distal     Exon   Intron Promoter 
#    3840     4226    29778     7553    35341    28449
prop.table(table(anno_df$region))
#     3'UTR      5'UTR     Distal       Exon     Intron   Promoter 
#0.03516902 0.03870424 0.27272477 0.06917490 0.32367406 0.26055300
data=as.data.frame(table(anno_df$region))
class_colors <- c(
 "3'UTR" = "#EF3B2C",
 "5'UTR" = "#FA9FB5",
 "Distal" = "#41AB5D",
 "Exon" = "#67A9CF",
 "Intron" = "#6A51A3",
 "Promoter" = "#FEB24C"
)

####################################################pieplot
total <- sum(data$Freq)  

# percentage
data$percent <- data$Freq / total * 100

##order
data$Var1=factor(data$Var1,levels=c("Exon","Promoter","5'UTR","Intron","3'UTR","Distal"))

library(ggplot2)

p <- ggplot(data, aes(x = 1, y = Freq, fill = Var1)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = class_colors,
                    guide = guide_legend(reverse = FALSE)) +
  ylab("") + xlab("") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.border     = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    axis.ticks       = element_blank()
  )
ggsave(p,
       filename="95sample_IterativeOverlapPeakSet.UCSC_hg19_known_gene.anno_by_chipseeker.genomic_anno_distribution.pieplot.pdf",
       width=6,height=5) ##Figure1B
q() 
