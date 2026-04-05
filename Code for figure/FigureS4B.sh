#######################
#######FigureS4B#######
#######################
conda activate /Public/gaoyun/miniconda3/envs/chipseeker
R
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
##peakfile
setwd("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4")
tissue_files <- list(
  pantissue  = "../10.final_IterativeOverlapPeakSet_pantissue/pantissue_IterativeOverlapPeakSet.sort.bed",
  adrenal.gland      = "./adrenal-gland/adrenal-gland_IterativeOverlapPeakSet.sort.bed",
  brain      = "./brain/brain_IterativeOverlapPeakSet.sort.bed",
  colon      = "./colon/colon_IterativeOverlapPeakSet.sort.bed",
  esophagus  = "./esophagus/esophagus_IterativeOverlapPeakSet.sort.bed",
  fallopian.tube  = "./fallopian-tube/fallopian-tube_IterativeOverlapPeakSet.sort.bed",
  fat  = "./fat/fat_IterativeOverlapPeakSet.sort.bed",
  heart  = "./heart/heart_IterativeOverlapPeakSet.sort.bed",
  kidney     = "./kidney/kidney_IterativeOverlapPeakSet.sort.bed",
  liver      = "./liver/liver_IterativeOverlapPeakSet.sort.bed",
  lung      = "./lung/lung_IterativeOverlapPeakSet.sort.bed",
  muscle      = "./muscle/muscle_IterativeOverlapPeakSet.sort.bed",
  nerve      = "./nerve/nerve_IterativeOverlapPeakSet.sort.bed",
  ovary      = "./ovary/ovary_IterativeOverlapPeakSet.sort.bed",
  pancreas   = "./pancreas/pancreas_IterativeOverlapPeakSet.sort.bed",
  Retina     = "./Retina/Retina_IterativeOverlapPeakSet.sort.bed",
  skin       = "./skin/skin_IterativeOverlapPeakSet.sort.bed",
  spleen       = "./spleen/spleen_IterativeOverlapPeakSet.sort.bed",
  thyroid.gland = "./thyroid-gland/thyroid-gland_IterativeOverlapPeakSet.sort.bed",
  stomach    = "./stomach/stomach_IterativeOverlapPeakSet.sort.bed"
)

#
OCRs <- lapply(tissue_files, readPeakFile)

#
OCRs_anno <- lapply(OCRs, function(peaks) {
  annotatePeak(peaks,
               level = "gene",
               tssRegion = c(-1000, 1000),
               TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
})
#
pantissue_OCR_anno_df = as.data.frame(OCRs_anno["pantissue"])
library(dplyr)
pantissue_OCR_anno_df <- pantissue_OCR_anno_df %>%
  rename_with(~gsub("^pantissue\\.", "", .))
pantissue_OCR_anno_df$V4=paste0(pantissue_OCR_anno_df$seqnames,":",pantissue_OCR_anno_df$start,"-",pantissue_OCR_anno_df$end)
write.table(pantissue_OCR_anno_df,"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/iterative_peak_filtering/10.final_IterativeOverlapPeakSet_pantissue/pantissue_IterativeOverlapPeakSet.UCSC_hg19_known_gene.anno_by_chipseeker.txt",col.names=T,row.names=F,quote=F,sep="\t")

##
library(dplyr)
library(stringr)

#
process_annotation <- function(anno) {
  anno_df <- as.data.frame(anno@anno) %>%
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
  
  return(anno_df)
}

#
anno_results <- lapply(OCRs_anno, process_annotation)

#
region_summary <- lapply(anno_results, function(df) {
  table(df$region)
})

#
region_summary_df <- as.data.frame(do.call(rbind, region_summary))
region_summary_df$tissue=rownames(region_summary_df)

#
library(tidyr)
data_long <- region_summary_df %>%
  pivot_longer(
    cols = -tissue,         
    names_to = "region",   
    values_to = "Freq" 
  )
#
data_long$tissue <- factor(data_long$tissue,
                           levels = c("pantissue","adrenal.gland","brain","colon","esophagus","fallopian.tube","fat","heart","kidney","liver","lung","muscle",
                                      "nerve","ovary","pancreas","Retina","skin","spleen","thyroid.gland","stomach"))
##
class_colors <- c(
 "3'UTR" = "#EF3B2C",
 "5'UTR" = "#FA9FB5",
 "Distal" = "#41AB5D",
 "Exon" = "#67A9CF",
 "Intron" = "#6A51A3",
 "Promoter" = "#FEB24C"
)
#
library(scales)
library(ggplot2)
#
total_counts <- data_long %>%
  group_by(tissue) %>%
  summarize(Total = sum(Freq), .groups = 'drop')  
total_counts
## A tibble: 20 × 2
#   tissue          Total
#   <fct>           <int>
# 1 pantissue      375455
# 2 adrenal.gland   67428
# 3 brain          181131
# 4 colon           75613
# 5 esophagus       37812
# 6 fallopian.tube  40134
# 7 fat             32782
# 8 heart           87690
# 9 kidney          65068
#10 liver           88905
#11 lung            58058
#12 muscle          59438
#13 nerve           38213
#14 ovary           49449
#15 pancreas        75147
#16 Retina          77317
#17 skin            22875
#18 spleen          35002
#19 thyroid.gland   44348
#20 stomach        109187

#
total_counts <- data_long %>%
  group_by(tissue) %>%
  summarise(Total = sum(Freq))

#
data_long$region=factor(data_long$region,levels=c("Exon","Promoter","5'UTR","Intron","3'UTR","Distal"))

######plot
setwd("/data1/gy/ATAC_for_review/FigureS4B/output")
p <- ggplot(data_long, aes(x = tissue, y = Freq, fill = region)) +  
  geom_col(width = 0.7) +  
  scale_fill_manual(values = class_colors, guide = guide_legend(reverse = FALSE)) +  
  ylab("Number of curated OCRs") +  
  xlab("") +  
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(total_counts$Total) * 1.1)) +  
  theme_minimal(base_size = 14) +  
  theme(  
    panel.grid       = element_blank(),  
    panel.background = element_blank(),  
    plot.background  = element_blank(),  
    panel.border     = element_blank(),  
    axis.line        = element_line(size = 0.5),  
    axis.ticks.x     = element_line(size = 0.5),  
    axis.ticks.y     = element_line(size = 0.5),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1)  
  ) +
  #
  geom_text(
    data = total_counts,
    aes(x = tissue, y = Total, label = Total),
    vjust = -0.3,
    size = 2,
    inherit.aes = FALSE 
  )
ggsave(p,
       filename="pantissue_IterativeOverlapPeakSet.UCSC_hg19_known_gene.anno_by_chipseeker.genomic_anno_distribution.barplot.pdf",
       width=8,height=6)
