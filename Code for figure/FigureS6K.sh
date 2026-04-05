#############################################
##################FigureS6K##################
#############################################
###########################
#####1.ChIPseeker anno#####
###########################
conda activate /Public/gaoyun/miniconda3/envs/chipseeker
R
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

##asOCR
library(data.table)
asOCR=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1.bed")
anno_df$V4=paste0(anno_df$seqnames,":",anno_df$start,"-",anno_df$end)
anno_df_asOCR=subset(anno_df,anno_df$V4 %in% asOCR$V4)
dir.create("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/chipseeker_anno")
fwrite(anno_df_asOCR,
       "/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/chipseeker_anno/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1.UCSC_hg19_known_gene.anno_by_chipseeker.txt",
       sep="\t",quote=F,col.names=T,row.names=F)
##non-asOCR
anno_df_nonasOCR=subset(anno_df,!(anno_df$V4 %in% asOCR$V4))

##
table(anno_df_asOCR$region)
#   3'UTR    5'UTR   Distal     Exon   Intron Promoter 
#     121      193     1891      286     1858     1309
table(anno_df_nonasOCR$region)
#   3'UTR    5'UTR   Distal     Exon   Intron Promoter 
#    3719     4033    27887     7267    33483    27140

##
data_asOCR=as.data.frame(table(anno_df_asOCR$region))
data_asOCR$class="asOCR"
data_nonasOCR=as.data.frame(table(anno_df_nonasOCR$region))
data_nonasOCR$class="non-asOCR"
##
data_asOCR <- data_asOCR %>%
  group_by(class) %>%
  mutate(pct = Freq / sum(Freq) * 100)
data_nonasOCR <- data_nonasOCR %>%
  group_by(class) %>%
  mutate(pct = Freq / sum(Freq) * 100)
data=rbind(data_asOCR,data_nonasOCR)

##############################
#####2.Enrichment barplot#####
##############################
library(dplyr)  
library(tidyr)  
library(purrr)
#
wide <- data %>%  
  select(Var1, Freq, class) %>%  
  pivot_wider(  
    names_from  = class,  
    values_from = Freq  
  ) %>%  
  rename(  
    asOCR    = asOCR,  
    nonasOCR = `non-asOCR`  
  )  

#
total_asOCR    <- sum(wide$asOCR)  
total_nonasOCR <- sum(wide$nonasOCR)  

#
res <- wide %>%  
  rowwise() %>%  
  mutate(  
    #
    tab = list(matrix(  
      c(asOCR,  
        total_asOCR - asOCR,  
        nonasOCR,  
        total_nonasOCR - nonasOCR),  
      nrow = 2,  
      byrow = TRUE,  
      dimnames = list(  
        sample = c("asOCR", "non-asOCR"),  
        status = c("inRegion", "outRegion")  
      )  
    )),  
    ft = list(fisher.test(tab)),  
    OR      = ft$estimate,  
    p.value = ft$p.value,
    CI_Lower = ft$conf.int[1],  
    CI_Upper = ft$conf.int[2]
  ) %>%  
  ungroup() %>%  
  select(Var1, asOCR, nonasOCR, OR, p.value,CI_Lower,CI_Upper)  

##
library(dplyr)
#
res <- res %>%     
  mutate(  
    sig = case_when(  
      p.value < 0.0001 ~ "****",
      p.value < 0.001 ~ "***",  
      p.value < 0.01  ~ "**",  
      p.value < 0.05  ~ "*",  
      TRUE         ~ "ns"  
    )  
  )  

#
res <- res %>%   
  mutate(
    Var1 = factor(Var1, levels = Var1[order(OR, decreasing = T)]),
    bar_start      = ifelse(OR >= 1, 1, OR),  
    bar_end        = ifelse(OR >= 1, OR, 1),  
  ) 
fwrite(res,
       "/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/chipseeker_anno/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1.UCSC_hg19_known_gene.anno_by_chipseeker.genomic_anno_enrichment.result",
       quote=F,sep="\t",col.names=T,row.names=F)

##
library(ggplot2)
setwd("/data1/gy/ATAC_for_review/FigureS6K/output")
class_colors <- c(
 "3'UTR" = "#EF3B2C",
 "5'UTR" = "#FA9FB5",
 "Distal" = "#41AB5D",
 "Exon" = "#67A9CF",
 "Intron" = "#6A51A3",
 "Promoter" = "#FEB24C"
)
p <- ggplot(res, aes(x = Var1)) +  
  geom_rect(aes(  
    xmin = as.numeric(Var1) - 0.35,  
    xmax = as.numeric(Var1) + 0.35,  
    ymin = bar_start,  
    ymax = bar_end,  
    fill = Var1  
  )) +  
  geom_errorbar(aes(  
    x    = as.numeric(Var1),  
    ymin = CI_Lower,  
    ymax = CI_Upper  
  ), width = 0.1, color = "black") +  
  geom_hline(yintercept = 1, color = "black") +  
  #
  geom_text(aes(  
    x     = as.numeric(Var1),  
    #
    y     = ifelse(OR >= 1, CI_Upper + 0.02, CI_Lower - 0.02),  
    label = sig  
  ), size = 3) +  
  scale_x_discrete() +  
  scale_fill_manual(values = class_colors) +  
  scale_y_continuous(  
    limits = c(min(res$CI_Lower)-0.02, max(res$CI_Upper)+0.02)  
  ) +  
  theme_minimal() +  
  theme(  
    panel.grid    = element_blank(),  
    axis.line     = element_line(color = "black"),  
    axis.ticks    = element_line(color = "black"),  
    axis.text.x   = element_text(size = 11,angle=30,hjust=1),  
    axis.text.y   = element_text(size = 11),  
    axis.title.y  = element_text(size = 12),  
    axis.title.x  = element_text(size = 12),  
    legend.position = "none"  
  ) +  
  labs(  
    x = "Genomic annotation",  
    y = "OR(compared with non-asOCRs)"  
  )
ggsave(p,
       filename="95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1.UCSC_hg19_known_gene.anno_by_chipseeker.genomic_anno_enrichment.pdf",
       width = 6, height = 5) ##FigureS6K
