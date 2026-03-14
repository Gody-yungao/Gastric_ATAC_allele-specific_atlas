#########################################
############Figure1C&FigureS2C###########
#########################################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATAC_for_review/Figure1C&FigureS2C/output")
library(data.table)
#####background region generation
library(BSgenome.Hsapiens.UCSC.hg19)  
library(GenomicRanges)  
library(regioneR)  

#  
OCR <- rtracklayer::import("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed")  
#   
## 1. valid_chrs
valid_chrs <- paste0("chr", 1:22)  
## 2. chrom_len
chrom_len  <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[valid_chrs]  
## 3. mini_genome
mini_genome <- GRanges(seqnames = valid_chrs,  
                       ranges   = IRanges(start = 1, end = chrom_len),  
                       seqinfo  = Seqinfo(seqnames   = valid_chrs,  
                                          seqlengths = chrom_len,  
                                          genome     = "hg19"))  
## 4. random_bg
set.seed(123)  
random_bg <- randomizeRegions(OCR,  
                              genome  = mini_genome,
                              allow.overlaps = FALSE)
random_bg_df <- as.data.frame(random_bg)[,c("seqnames","start","end")]
random_bg_df$name=paste0(random_bg_df$seqnames,":",random_bg_df$start,"-",random_bg_df$end)
#sort
chr.order <- paste0("chr", 1:22)   
random_bg_df_sorted <- random_bg_df[  
  order( match(random_bg_df$seqnames, chr.order),   
         random_bg_df$start                         
  ),  
]  
#write.table(random_bg_df_sorted, file="/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.bg_random.sort.bed", 
#            quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
#
state_mapping1 <- c(  
  '2' = 'Active TSS',  
  '1' = 'Flanking active TSS',  
  '7' = 'Poised promoter',  
  '4' = 'Active enhancer',  
  '6' = 'Weak enhancer',  
  '3' = 'Poised enhancer',  
  '8' = 'Repressed polycomb',  
  '9' = 'Repressed polycomb',  
  '11' = 'Heterochromatin',  
  '5' = 'Quiescent',  
  '10' = 'Quiescent'  
)  
state_mapping2 <- c(  
  '2' = 'Promoter',  
  '1' = 'Promoter',  
  '7' = 'Promoter',  
  '4' = 'Enhancer',  
  '6' = 'Enhancer',  
  '3' = 'Enhancer',  
  '8' = 'Repressed',  
  '9' = 'Repressed',  
  '11' = 'Heterochromatin',  
  '5' = 'Quiescent',  
  '10' = 'Quiescent'  
)  
#stomach_ENCODE
state_ENCODE=fread("/data1/gy/chromHMM_stomach/segmentation/state11/stomach_37Y-ENCODE_11_stomach_dense.bed")
state_ENCODE=state_ENCODE[,1:4]
colnames(state_ENCODE)=c("chr","start","end","state")
state_ENCODE[, state_name1 := state_mapping1[as.character(state)]]  
state_ENCODE[, state_name2 := state_mapping2[as.character(state)]]  
#stomach_NJ-GaEpi
state_NJMU=fread("/data1/gy/chromHMM_stomach/segmentation/state11/1443364-N_11_stomach_dense.bed")
state_NJMU=state_NJMU[,1:4]
colnames(state_NJMU)=c("chr","start","end","state")
state_NJMU[, state_name1 := state_mapping1[as.character(state)]]  
state_NJMU[, state_name2 := state_mapping2[as.character(state)]]  
#OCR
OCR=fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed")
OCR$V4=paste0(OCR$V1,":",OCR$V2+1,"-",OCR$V3)
OCR$V3=OCR$V3-250
OCR$V2=OCR$V3-1
#bg
bg=fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.bg_random.sort.bed")
bg$V3=bg$V3-250
bg$V2=bg$V3-1

#######################anno
library(GenomicRanges) 
# OCR
gr_OCR <- GRanges(  
  seqnames = OCR$V1,  
  ranges = IRanges(start = OCR$V2+1, end = OCR$V3),
  name = OCR$V4,
  SPM = OCR$V5
)  
# bg
gr_bg <- GRanges(  
  seqnames = bg$V1,  
  ranges = IRanges(start = bg$V2+1, end = bg$V3),
  name = bg$V4
)  
# state_ENCODE
gr_state_ENCODE <- GRanges(  
  seqnames = state_ENCODE$chr,  
  ranges = IRanges(start = state_ENCODE$start+1, end = state_ENCODE$end),  
  state_name1 = state_ENCODE$state_name1,  
  state_name2 = state_ENCODE$state_name2  
)  
# state_NJ-GaEpi
gr_state_NJMU <- GRanges(  
  seqnames = state_NJMU$chr,  
  ranges = IRanges(start = state_NJMU$start+1, end = state_NJMU$end),  
  state_name1 = state_NJMU$state_name1,  
  state_name2 = state_NJMU$state_name2  
)  

# overlap  
hits_gr_OCR_with_state_ENCODE <- findOverlaps(gr_OCR, gr_state_ENCODE)
hits_gr_OCR_with_state_NJMU <- findOverlaps(gr_OCR, gr_state_NJMU)  
hits_gr_bg_with_state_ENCODE <- findOverlaps(gr_bg, gr_state_ENCODE)
hits_gr_bg_with_state_NJMU <- findOverlaps(gr_bg, gr_state_NJMU)  


###############################ENCODE
##OCR
df_annot_OCR_with_state_ENCODE <- data.frame(  
    OCR_index = queryHits(hits_gr_OCR_with_state_ENCODE),  
    state_index = subjectHits(hits_gr_OCR_with_state_ENCODE)  
)  
OCR_anno_with_state_ENCODE <- cbind(  
    OCR[df_annot_OCR_with_state_ENCODE$OCR_index, ],  
    state_ENCODE[df_annot_OCR_with_state_ENCODE$state_index, c("state_name1", "state_name2")]  
)  
table(OCR_anno_with_state_ENCODE$state_name2)
#       Enhancer Heterochromatin        Promoter       Quiescent       Repressed 
#          48213             497           40816           14601            5060
table(OCR_anno_with_state_ENCODE$state_name1)
#    Active enhancer          Active TSS Flanking active TSS     Heterochromatin 
#              29727               18433               14515                 497 
#    Poised enhancer     Poised promoter           Quiescent  Repressed polycomb 
#               2417                7868               14601                5060 
#      Weak enhancer 
#              16069
#fwrite(OCR_anno_with_state_ENCODE,
#       "/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/chromHMM_overlap/95sample_IterativeOverlapPeakSet.ENCODE_chromHMM_state.txt",
#       sep="\t",quote=F,col.names=T,row.names=F)
##bg
df_annot_bg_with_state_ENCODE <- data.frame(  
    bg_index = queryHits(hits_gr_bg_with_state_ENCODE),  
    state_index = subjectHits(hits_gr_bg_with_state_ENCODE)  
)  
bg_anno_with_state_ENCODE <- cbind(  
    OCR[df_annot_bg_with_state_ENCODE$OCR_index, ],  
    state_ENCODE[df_annot_bg_with_state_ENCODE$state_index, c("state_name1", "state_name2")]  
)  
table(bg_anno_with_state_ENCODE$state_name2)
#     Enhancer Heterochromatin        Promoter       Quiescent       Repressed 
#         7769            4591            1886           83046           11895 
table(bg_anno_with_state_ENCODE$state_name1)
#    Active enhancer          Active TSS Flanking active TSS     Heterochromatin 
#               2666                 543                 753                4591 
#    Poised enhancer     Poised promoter           Quiescent  Repressed polycomb 
#                621                 590               83046               11895 
#      Weak enhancer 
#               4482 

################################NJ-GaEpi
##OCR
df_annot_OCR_with_state_NJMU <- data.frame(  
    OCR_index = queryHits(hits_gr_OCR_with_state_NJMU),  
    state_index = subjectHits(hits_gr_OCR_with_state_NJMU)  
)  
OCR_anno_with_state_NJMU <- cbind(  
    OCR[df_annot_OCR_with_state_NJMU$OCR_index, ],  
    state_NJMU[df_annot_OCR_with_state_NJMU$state_index, c("state_name1", "state_name2")]  
)  
table(OCR_anno_with_state_NJMU$state_name2)
#       Enhancer Heterochromatin        Promoter       Quiescent       Repressed 
#          36115             335           58803            9880            4054 
table(OCR_anno_with_state_NJMU$state_name1)
#    Active enhancer          Active TSS Flanking active TSS     Heterochromatin 
#              10722               29698               22165                 335 
#    Poised enhancer     Poised promoter           Quiescent  Repressed polycomb 
#              13895                6940                9880                4054 
#      Weak enhancer 
#              11498 
#fwrite(OCR_anno_with_state_NJMU,
#       "/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/chromHMM_overlap/95sample_IterativeOverlapPeakSet.NJMU_chromHMM_state.txt",
#       sep="\t",quote=F,col.names=T,row.names=F)

##bg
df_annot_bg_with_state_NJMU <- data.frame(  
    bg_index = queryHits(hits_gr_bg_with_state_NJMU),  
    state_index = subjectHits(hits_gr_bg_with_state_NJMU)  
)  
bg_anno_with_state_NJMU <- cbind(  
    bg[df_annot_bg_with_state_NJMU$bg_index, ],  
    state_NJMU[df_annot_bg_with_state_NJMU$state_index, c("state_name1", "state_name2")]  
)  
table(bg_anno_with_state_NJMU$state_name2)
#       Enhancer Heterochromatin        Promoter       Quiescent       Repressed 
#           4224             798            1840           92297           10028
table(bg_anno_with_state_NJMU$state_name1)
#   Active enhancer          Active TSS Flanking active TSS     Heterochromatin 
#                509                1039                 500                 798 
#    Poised enhancer     Poised promoter           Quiescent  Repressed polycomb 
#               2462                 301               92297               10028 
#      Weak enhancer 
#               1253


##################################################
########OCR summit enrichment in chromHMM#########
##################################################
#######################################ENCODE
# 1. count  
ocr_counts <- table(OCR_anno_with_state_ENCODE$state_name1)  
bg_counts  <- table(bg_anno_with_state_ENCODE$state_name1)  

# 2. total num
total_ocr <- nrow(OCR_anno_with_state_ENCODE)  
total_bg  <- nrow(bg_anno_with_state_ENCODE)  

# 3. Fisher test
states <- union(names(ocr_counts), names(bg_counts))  

res_list <- lapply(states, function(st) {  
  a <- ifelse(!is.na(ocr_counts[st]), ocr_counts[st], 0)  # OCR  
  b <- ifelse(!is.na(bg_counts[st]),  bg_counts[st],  0)  # bg  
  c <- total_ocr - a  
  d <- total_bg  - b    

  # 2x2 matirx: row OCR/bg，column in_state/not_in_state  
  mat <- matrix(c(a, c, b, d), nrow = 2,  
                dimnames = list(  
                  Group = c("OCR", "bg"),  
                  State = c("in", "out")  
                ))  
  ft <- fisher.test(mat)  
  
  data.frame(  
    state      = st,  
    ocr_in     = a,  
    bg_in      = b,  
    odds_ratio = unname(ft$estimate),  
    p.value    = ft$p.value,  
    conf_low   = ft$conf.int[1],  
    conf_high  = ft$conf.int[2],  
    stringsAsFactors = FALSE  
  )  
})  

# 4. res rbind  
res_df <- do.call(rbind, res_list)   

#  
library(dplyr)  
res_df <- res_df %>%     
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
res_df <- res_df %>%   
  mutate(  
    state = factor(state, levels = rev(c(  
      "Active TSS","Flanking active TSS","Poised promoter",  
      "Active enhancer","Weak enhancer","Poised enhancer",  
      "Repressed polycomb","Heterochromatin","Quiescent"  
    ))),  
    OR_log2        = log2(odds_ratio),  
    bar_start      = ifelse(OR_log2 >= 0, 0, OR_log2),  
    bar_end        = ifelse(OR_log2 >= 0, OR_log2, 0),  
    conf_high_log2 = log2(conf_high),  
    conf_low_log2  = log2(conf_low)
  ) %>%   
  arrange(state)  
#
res_df <- res_df %>%   
  mutate(  
    class          = rev(c(  
                       rep("Promoter",3),  
                       rep("Enhancer",3),  
                       "Repressed","Heterochromatin","Quiescent"  
                     ))
  )
#
class_colors <- c(  
  "Promoter"      = "#D73027",  
  "Enhancer"      = "#FDAE61",  
  "Repressed"     = "#3288BD",  
  "Heterochromatin" = "#5E4FA2",  
  "Quiescent"     = "#878787"  
)  

library(ggplot2)
p1 <- ggplot(res_df, aes(x = state)) +  
  geom_rect(aes(  
    xmin = as.numeric(state) - 0.35,  
    xmax = as.numeric(state) + 0.35,  
    ymin = bar_start,  
    ymax = bar_end,  
    fill = class  
  )) +  
  geom_errorbar(aes(  
    x    = as.numeric(state),  
    ymin = conf_low_log2,  
    ymax = conf_high_log2  
  ), width = 0.1, color = "black") +  
  geom_hline(yintercept = 0, color = "black") +  
  # 
  geom_text(aes(  
    x     = as.numeric(state),  
    # 
    y     = ifelse(OR_log2 >= 0, bar_end + 0.2, bar_start - 0.2),  
    label = sig  
  ), size = 5) +  
  scale_x_discrete() +  
  scale_fill_manual(values = class_colors) +  
  scale_y_continuous(  
    limits = c(min(res_df$conf_low_log2)-0.5, max(res_df$conf_high_log2)+0.5)  
  ) +  
  theme_minimal() +  
  theme(  
    panel.grid    = element_blank(),  
    axis.line     = element_line(color = "black"),  
    axis.ticks    = element_line(color = "black"),  
    axis.text.x   = element_text(size = 11),  
    axis.text.y   = element_text(size = 11),  
    axis.title.y  = element_text(size = 12),  
    axis.title.x  = element_text(size = 12),  
    legend.position = "none"  
  ) +  
  labs(  
    x = "Chromatin state (from ENCODE stomach)",  
    y = "log2(OR)"  
  ) +  
  coord_flip()  

# 
ggsave(p1, file="95sample_IterativeOverlapPeakSet.ENCODE_chromHMM_state.enrichment.barplot.pdf", 
       width = 6, height = 7)  ##FigureS2C

#######################################NJ-GaEpi
# 1. count
ocr_counts <- table(OCR_anno_with_state_NJMU$state_name1)  
bg_counts  <- table(bg_anno_with_state_NJMU$state_name1)  

# 2. total num
total_ocr <- nrow(OCR_anno_with_state_NJMU)  
total_bg  <- nrow(bg_anno_with_state_NJMU)  

# 3. Fisher test  
states <- union(names(ocr_counts), names(bg_counts))  

res_list <- lapply(states, function(st) {  
  a <- ifelse(!is.na(ocr_counts[st]), ocr_counts[st], 0)  # OCR
  b <- ifelse(!is.na(bg_counts[st]),  bg_counts[st],  0)  # bg
  c <- total_ocr - a 
  d <- total_bg  - b 

  # 2x2 matirx: row OCR/bg，column in_state/not_in_state  
  mat <- matrix(c(a, c, b, d), nrow = 2,  
                dimnames = list(  
                  Group = c("OCR", "bg"),  
                  State = c("in", "out")  
                ))  
  ft <- fisher.test(mat)  
  
  data.frame(  
    state      = st,  
    ocr_in     = a,  
    bg_in      = b,  
    odds_ratio = unname(ft$estimate),  
    p.value    = ft$p.value,  
    conf_low   = ft$conf.int[1],  
    conf_high  = ft$conf.int[2],  
    stringsAsFactors = FALSE  
  )  
})  

# 4. res rbind
res_df <- do.call(rbind, res_list)   

#
library(dplyr)
res_df <- res_df %>%     
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
res_df <- res_df %>%   
  mutate(  
    state = factor(state, levels = rev(c(  
      "Active TSS","Flanking active TSS","Poised promoter",  
      "Active enhancer","Weak enhancer","Poised enhancer",  
      "Repressed polycomb","Heterochromatin","Quiescent"  
    ))),  
    OR_log2        = log2(odds_ratio),  
    bar_start      = ifelse(OR_log2 >= 0, 0, OR_log2),  
    bar_end        = ifelse(OR_log2 >= 0, OR_log2, 0),  
    conf_high_log2 = log2(conf_high),  
    conf_low_log2  = log2(conf_low)
  ) %>%   
  arrange(state)  
#
res_df <- res_df %>%   
  mutate(  
    class          = rev(c(  
                       rep("Promoter",3),  
                       rep("Enhancer",3),  
                       "Repressed","Heterochromatin","Quiescent"  
                     ))
  )
#
class_colors <- c(  
  "Promoter"      = "#D73027",  
  "Enhancer"      = "#FDAE61",  
  "Repressed"     = "#3288BD",  
  "Heterochromatin" = "#5E4FA2",  
  "Quiescent"     = "#878787"  
)  

##
p2 <- ggplot(res_df, aes(x = state)) +  
  geom_rect(aes(  
    xmin = as.numeric(state) - 0.35,  
    xmax = as.numeric(state) + 0.35,  
    ymin = bar_start,  
    ymax = bar_end,  
    fill = class  
  )) +  
  geom_errorbar(aes(  
    x    = as.numeric(state),  
    ymin = conf_low_log2,  
    ymax = conf_high_log2  
  ), width = 0.1, color = "black") +  
  geom_hline(yintercept = 0, color = "black") +  
  # 
  geom_text(aes(  
    x     = as.numeric(state),  
    #
    y     = ifelse(OR_log2 >= 0, bar_end + 0.2, bar_start - 0.2),  
    label = sig  
  ), size = 5) +  
  scale_x_discrete() +  
  scale_fill_manual(values = class_colors) +  
  scale_y_continuous(  
    limits = c(min(res_df$conf_low_log2)-0.5, max(res_df$conf_high_log2)+0.5)  
  ) +  
  theme_minimal() +  
  theme(  
    panel.grid    = element_blank(),  
    axis.line     = element_line(color = "black"),  
    axis.ticks    = element_line(color = "black"),  
    axis.text.x   = element_text(size = 11),  
    axis.text.y   = element_text(size = 11),  
    axis.title.y  = element_text(size = 12),  
    axis.title.x  = element_text(size = 12),  
    legend.position = "none"  
  ) +  
  labs(  
    x = "Chromatin state (from NJMU stomach)",  
    y = "log2(OR)"  
  ) +  
  coord_flip()  

#
ggsave(p2, file="95sample_IterativeOverlapPeakSet.NJ-GaEpi_chromHMM_state.enrichment.barplot.pdf", 
       width = 6, height = 7)  ##Figure1C

 
