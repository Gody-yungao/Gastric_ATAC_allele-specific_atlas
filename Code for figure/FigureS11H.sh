##########################
########FigureS11H########
##########################
##############################################
########1.DiffBind: diff atac analysis########
##############################################
conda activate /Public/gaoyun/miniconda3/envs/diffbind
mkdir /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/siGATA_ATAC/20251115/diff_peak
R
library(DiffBind)
rm(list = ls())
setwd("/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/siGATA_ATAC/20251115/diff_peak")

##############################
SampleID <- c(paste("NC",1:2,sep = "-"),paste("siGATA4_3",1:2,sep = "-"),paste("siGATA6_2",1:2,sep = "-"))
Tissue <- rep("AGS",times = 6)
Factor <- rep(NA,times = 6)
Condition <- rep(NA,times = 6)
Treatment <- rep(c("NC","siGATA4","siGATA6"),c(2,2,2))
Replicate <- rep(1:2,length = 6)
PeakCaller <- rep("narrow",times = 6)

bam_file_path <- "/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/siGATA_ATAC/20251115/mapping"
bamReads <- c(paste(bam_file_path,"AGS-GATA-NC-1.last.shift.sort.bam",sep = "/"),
              paste(bam_file_path,"AGS-GATA-NC-2.last.shift.sort.bam",sep = "/"),
              paste(bam_file_path,"AGS-GATA-4-3-1.last.shift.sort.bam",sep = "/"),
              paste(bam_file_path,"AGS-GATA-4-3-2.last.shift.sort.bam",sep = "/"),
              paste(bam_file_path,"AGS-GATA-6-2-1.last.shift.sort.bam",sep = "/"),
              paste(bam_file_path,"AGS-GATA-6-2-2.last.shift.sort.bam",sep = "/"))

peak_file_path <- "/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/siGATA_ATAC/20251115/macs2_peak"
Peaks <- c(paste(peak_file_path,"AGS-GATA-NC-1_peaks.narrowPeak",sep = "/"),
           paste(peak_file_path,"AGS-GATA-NC-2_peaks.narrowPeak",sep = "/"),
           paste(peak_file_path,"AGS-GATA-4-3-1_peaks.narrowPeak",sep = "/"),
           paste(peak_file_path,"AGS-GATA-4-3-2_peaks.narrowPeak",sep = "/"),
           paste(peak_file_path,"AGS-GATA-6-2-1_peaks.narrowPeak",sep = "/"),
           paste(peak_file_path,"AGS-GATA-6-2-2_peaks.narrowPeak",sep = "/"))

samples <- data.frame(SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,Peaks,PeakCaller)

##############################
atac <- dba(sampleSheet = samples,minOverlap = 2)

##############################
atac_count <- dba.count(atac,minOverlap = 2,summits=250)
info = dba.show(atac_count)
libsizes = cbind(LibReads=info$Reads, FRiP=info$FRiP,
                  PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) = info$ID

##############################
atac_count_norm <- dba.normalize(atac_count,normalize=DBA_NORM_LIB)
norm <- dba.normalize(atac_count_norm, bRetrieve=TRUE)
normlibs <- cbind(FullLibSize=norm$lib.sizes, 
                  NormFacs=norm$norm.factors, 
                  NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID

##############################
atac_diff_contrast <- dba.contrast(atac_count_norm,
                      minMembers=2,  
                      reorderMeta = list(Treatment="NC"))
atac_diff_contrast
#Design: [~Treatment] | 3 Contrasts:
#     Factor   Group Samples  Group2 Samples2
#1 Treatment      NC       2 siGATA4        2
#2 Treatment      NC       2 siGATA6        2
#3 Treatment siGATA6       2 siGATA4        2

##############################
atac_diff <- dba.analyze(atac_diff_contrast, bRetrieveAnalysis=F)
dba.show(atac_diff, bContrasts=TRUE)
#Design: [~Treatment] | 3 Contrasts:
#     Factor   Group Samples  Group2 Samples2 DB.DESeq2
#1 Treatment      NC       2 siGATA4        2     57555
#2 Treatment      NC       2 siGATA6        2     35245
#3 Treatment siGATA6       2 siGATA4        2     32328
  
#
atac_diff.DB.siGATA4 <- dba.report(atac_diff,contrast=1)
atac_diff.DB.siGATA4.DF <- as.data.frame(atac_diff.DB.siGATA4)
atac_diff.DB.siGATA6 <- dba.report(atac_diff,contrast=2)
atac_diff.DB.siGATA6.DF <- as.data.frame(atac_diff.DB.siGATA6)

##########################################siGATA4 vs. siNC
atac_diff.DB_all.siGATA4 <- dba.report(atac_diff,th=1,method=DBA_DESEQ2, contrast = 1,bCounts=TRUE,bNormalized=TRUE)
atac_diff.DB_all.siGATA4.DF <- as.data.frame(atac_diff.DB_all.siGATA4)
atac_diff.DB_all.siGATA4.DF$peak=paste0(atac_diff.DB_all.siGATA4.DF$seqnames,":",atac_diff.DB_all.siGATA4.DF$start,"-",atac_diff.DB_all.siGATA4.DF$end)
atac_diff.DB_all.siGATA4.DF$Fold=-atac_diff.DB_all.siGATA4.DF$Fold 
write.table(atac_diff.DB_all.siGATA4.DF,file = "AGS_siGATA4_diff_peaks_ATACseq_all_with_normreadcounts.txt",
            quote = F,sep = "\t",
            row.names = F,col.names = T)
subset(atac_diff.DB_all.siGATA4.DF,peak=="chr4:48073178-48074157") ##rs875179 local OCR
#       seqnames    start      end width strand     Conc  Conc_NC Conc_siGATA4
#146426     chr4 48073178 48074157   980      * 7.316372 8.099047     5.478235
#            Fold     p.value          FDR   NC.1  NC.2 siGATA4_3.1 siGATA4_3.2
#146426 -2.573274 2.36803e-30 6.996361e-28 279.29 269.1       39.74       49.41
#                         peak
#146426 chr4:48073178-48074157

##########################################siGATA6 vs. siNC
atac_diff.DB_all.siGATA6 <- dba.report(atac_diff,th=1,method=DBA_DESEQ2, contrast = 2,bCounts=TRUE,bNormalized=TRUE)
atac_diff.DB_all.siGATA6.DF <- as.data.frame(atac_diff.DB_all.siGATA6)
atac_diff.DB_all.siGATA6.DF$peak=paste0(atac_diff.DB_all.siGATA6.DF$seqnames,":",atac_diff.DB_all.siGATA6.DF$start,"-",atac_diff.DB_all.siGATA6.DF$end)
atac_diff.DB_all.siGATA6.DF$Fold=-atac_diff.DB_all.siGATA6.DF$Fold
write.table(atac_diff.DB_all.siGATA6.DF,file = "AGS_siGATA6_diff_peaks_ATACseq_all_with_normreadcounts.txt",
            quote = F,sep = "\t",
            row.names = F,col.names = T)
subset(atac_diff.DB_all.siGATA6.DF,peak=="chr4:48073178-48074157") ##rs875179 local OCR
#       seqnames    start      end width strand     Conc  Conc_NC Conc_siGATA6
#146426     chr4 48073178 48074157   980      * 7.791615 8.099047      7.40042
#             Fold      p.value         FDR   NC.1  NC.2 siGATA6_2.1 siGATA6_2.2
#146426 -0.6235045 0.0001347932 0.001975932 279.29 269.1      181.68      156.22
#                         peak
#146426 chr4:48073178-48074157
q()

###################################################
########2. boxplot for siNC siGATA4 siGATA6########
###################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(ggplot2)
library(ggpubr)
#
dt <- data.frame(
  Group = c("NC", "NC", "siGATA4", "siGATA4", "siGATA6", "siGATA6"),
  Value = c(279.29, 269.1, 39.74, 49.41, 181.68, 156.22)
)

#
dt$Group <- factor(dt$Group, levels = c("NC", "siGATA4", "siGATA6"))

#
comparisons <- list( c("NC", "siGATA4"), c("NC", "siGATA6") )

#######
setwd("/data1/gy/ATAC_for_review/FigureS11H/output")
p <- ggplot(dt, aes(x = Group, y = Value)) +  
  geom_boxplot(  
    aes(color = Group),   
    fill = "white",  
    width = 0.6,  
    outlier.shape = NA  
  ) +  
  stat_compare_means(  
    comparisons = comparisons,  
    method = "t.test",  
    label = "p.format",
    label.sep = " ",
    step.increase = 0.15,
    size = 3.5,
    vjust = 0.5
  ) +  
  #
  scale_color_manual(values = c("NC" = "#bcbcbc", 
                                "siGATA4" = "#61b096", 
                                "siGATA6" = "#e0816c")) +  
  #
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(  
    x = "",  
    y = "Normalized ATAC read count"  
  ) 
  
p <- p +
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
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8),  
    axis.text.y = element_text(size = 8),  
    legend.position = "none"  
  )  
ggsave(p,file="AGS_siGATA4_siGATA6_diff_peaks_ATACseq.boxplot.pdf",width=3,height=4) ##FigureS11H
q()
