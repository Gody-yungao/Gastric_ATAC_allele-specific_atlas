############################
#########FigureS4A##########
############################
conda activate R_base
R
setwd("/data1/gy/ATAC_for_review/FigureS4A/output")
library(data.table)
################################
#########pan-tissue QC##########
################################
#############################################
###############i.TSS enrichment##############
#############################################
GEO_brain=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/TSS_enrichment_score/brain/GEO_brain.tss_enrich.qc")
GEO_colon=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/TSS_enrichment_score/colon/GEO_colon.tss_enrich.qc")
GEO_kidney=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/TSS_enrichment_score/kidney/GEO_kidney.tss_enrich.qc")
GEO_liver=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/TSS_enrichment_score/liver/GEO_liver.tss_enrich.qc")
GEO_pancreas=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/TSS_enrichment_score/pancreas/GEO_pancreas.tss_enrich.qc")
GEO_retina=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/TSS_enrichment_score/retina/GEO_retina.tss_enrich.qc")
GEO_skin=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/TSS_enrichment_score/skin/GEO_skin.tss_enrich.qc")
ERR_esophagus=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/ERR/TSS_enrichment_score/esophagus/ERR_esophagus.tss_enrich.qc")
ERR_esophagus$origin="ERR"
ENCODE=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/ENCODE/TSS_enrichment_score/ENCODE.tss_enrich.qc")
ENCODE$origin="ENCODE"
pantissue_TSS=rbind(GEO_brain,GEO_colon,GEO_kidney,GEO_liver,GEO_pancreas,GEO_retina,GEO_skin)
pantissue_TSS$origin="GEO"
pantissue_TSS=rbind(pantissue_TSS,ERR_esophagus,ENCODE)
colnames(pantissue_TSS)=c("sample","TSS_enrichment_score","origin")
pantissue_TSS$summitfile=paste0(pantissue_TSS$origin,"_",pantissue_TSS$sample,"_summits.bed")

##################################################
###############ii.library complexity###############
##################################################
GEO_brain=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/library_complexity/brain/ATACseq_pbc.qc")
GEO_colon=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/library_complexity/colon/ATACseq_pbc.qc")
GEO_kidney=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/library_complexity/kidney/ATACseq_pbc.qc")
GEO_liver=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/library_complexity/liver/ATACseq_pbc.qc")
GEO_pancreas=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/library_complexity/pancreas/ATACseq_pbc.qc")
GEO_retina=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/library_complexity/retina/ATACseq_pbc.qc")
GEO_skin=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/library_complexity/skin/ATACseq_pbc.qc")
ERR_esophagus=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/ERR/library_complexity/esophagus/ATACseq_pbc.qc")
ENCODE=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/ENCODE/library_complexity/ATACseq_pbc.qc")
pantissue_pbc=rbind(GEO_brain,GEO_colon,GEO_kidney,GEO_liver,GEO_pancreas,GEO_retina,GEO_skin,ERR_esophagus,ENCODE)
pantissue_pbc=pantissue_pbc[,c(1,6,7,8)]
colnames(pantissue_pbc)=c("sample","NRF","PBC1","PBC2")

####################################
###############iii.FRiP###############
####################################
GEO_brain=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/FRiP/brain/FRiP.txt")
GEO_colon=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/FRiP/colon/FRiP.txt")
GEO_kidney=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/FRiP/kidney/FRiP.txt")
GEO_liver=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/FRiP/liver/FRiP.txt")
GEO_pancreas=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/FRiP/pancreas/FRiP.txt")
GEO_retina=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/FRiP/retina/FRiP.txt")
GEO_skin=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/GEO/FRiP/skin/FRiP.txt")
ERR_esophagus=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/ERR/FRiP/esophagus/FRiP.txt")
ENCODE=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/ENCODE/FRiP/FRiP.txt")
pantissue_FRiP=rbind(GEO_brain,GEO_colon,GEO_kidney,GEO_liver,GEO_pancreas,GEO_retina,GEO_skin,ERR_esophagus,ENCODE)
pantissue_FRiP=pantissue_FRiP[,c(1,4)]
colnames(pantissue_FRiP)=c("sample","FRiP")

###
pantissue=merge(pantissue_TSS,pantissue_pbc,by="sample")
pantissue=merge(pantissue,pantissue_FRiP,by="sample")

#
library(dplyr)
library(stringr)
pantissue <- pantissue %>%
  mutate(
    tissue = str_extract(sample, "^[^_]+")  
  ) %>%
  mutate(
    tissue = if_else(startsWith(sample, "RPE"), "Retina", tissue)
  )
pantissue=subset(pantissue,tissue != "stomach")
table(pantissue$tissue)
# adrenal-gland         artery      bile-duct          brain         breast 
#             8              4              3            140              3 
#         colon      esophagus fallopian-tube            fat          heart 
#            18              6              3              3             44 
#        kidney          liver           lung         muscle          nerve 
#             7             28              8              9              5 
#         ovary       pancreas         Retina           skin         spleen 
#             7             23             19             19              5 
# thyroid-gland           vein 
#             3              3
library(openxlsx)
write.xlsx(pantissue,"GEO_ERR_ENCODE_sample_bfQC.xlsx")
fwrite(pantissue,"GEO_ERR_ENCODE_sample_bfQC.txt",sep="\t",col.names=T,row.names=F,quote=F)

##QC
pantissueQC=subset(pantissue,PBC1 >= 0.7 & PBC2 >= 1 & NRF >= 0.7 & TSS_enrichment_score >=6 & FRiP >=0.2)
table(pantissueQC$tissue)
# adrenal-gland         artery      bile-duct          brain         breast 
#             8              1              1            137              2 
#         colon      esophagus fallopian-tube            fat          heart 
#            15              5              3              3             37 
#        kidney          liver           lung         muscle          nerve 
#             7             27              8              8              3 
#         ovary       pancreas         Retina           skin         spleen 
#             6             13             19              3              5 
# thyroid-gland           vein 
#             3              2
write.xlsx(pantissueQC,"GEO_ERR_ENCODE_sample_afterQC.xlsx")
fwrite(pantissueQC,"GEO_ERR_ENCODE_sample_afterQC.txt",sep="\t",col.names=T,row.names=F,quote=F)

##
setDT(pantissueQC)
pantissueQC <- pantissueQC[, if (.N >= 3) .SD, by = tissue] ##at least 3 samples for each tissue-type
table(pantissueQC$tissue)
# adrenal-gland          brain          colon      esophagus fallopian-tube 
#             8            137             15              5              3 
#           fat          heart         kidney          liver           lung 
#             3             37              7             27              8 
#        muscle          nerve          ovary       pancreas         Retina 
#             8              3              6             13             19 
#          skin         spleen  thyroid-gland 
#             3              5              3
write.xlsx(pantissueQC,"GEO_ERR_ENCODE_sample_afterQC.3more.final.xlsx")
fwrite(pantissueQC,"GEO_ERR_ENCODE_sample_afterQC.3more.final.txt",sep="\t",col.names=T,row.names=F,quote=F)


##############################
#########2.bfQC plot##########
##############################
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(patchwork)
library(colorspace)
#
my22colors <- c(
  '#f3715c', '#8e7437', '#74905d', '#494e8f', '#bd6758',
  '#fab27b', '#b2d235', '#DC050C', '#E8601C', '#F1932D',
  '#F6C141', '#F7EE55', '#CAE0AB', '#7BAFDE', '#5289C7',
  '#1965B0', '#D6C1DE', '#B178A6',
  '#882E72', '#8A2BE2',
  '#FF69B4', '#7CFC00'
)
#
col <- c(
  "adrenal-gland" = '#f3715c',
  "brain" = '#494e8f',
  "colon" = '#fab27b',
  "esophagus" = '#b2d235',
  "fallopian-tube" = '#DC050C',
  "fat" = '#E8601C',
  "heart" = '#F1932D',
  "kidney" = '#F6C141',
  "liver" = '#F7EE55',
  "lung" = '#CAE0AB',
  "muscle" = '#7BAFDE',
  "nerve" = '#5289C7',
  "ovary" = '#1965B0',
  "pancreas" = '#D6C1DE',
  "Retina" = '#B178A6',
  "skin" = '#882E72',
  "spleen" = '#8A2BE2',
  "thyroid-gland" = '#FF69B4',
  "stomach" = '#8e7437'
)
#
#darkened_my22colors <- darken(my22colors, amount = 0.2) 

#####
pantissue <- pantissue %>%
  mutate(
    QC = if_else(sample %in% pantissueQC$sample,
                 "high_quality", "low_quality")
  )
pantissue=subset(pantissue,pantissue$tissue %in% pantissueQC$tissue)

#####################################1.TSS enrichment
p1 <- ggplot(pantissue, aes(x = tissue, y = TSS_enrichment_score, fill = NULL, color=tissue)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(
    aes(shape = QC),  
    width = 0.2, size = 1.5, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  scale_y_continuous(limits = c(0, max(pantissue$TSS_enrichment_score) + 0.5)) + 
  geom_hline(yintercept = 6, linetype = "dashed", color = "grey", size = 0.5) +
  scale_color_manual(values = col, guide = "none") + 
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(y = "TSS enrichment") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), 
        legend.position = "right"  )

#####################################2.FRiP
p2 <- ggplot(pantissue, aes(x = tissue, y = FRiP, fill = NULL, color=tissue)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(
    aes(shape = QC),  
    width = 0.2, size = 1.5, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey", size = 0.5) + 
  scale_color_manual(values = col, guide = "none") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(y = "FRiP") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), 
        legend.position = "right"  ) 

#####################################3.NRF
p3 <- ggplot(pantissue, aes(x = tissue, y = NRF, fill = NULL, color=tissue)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(
    aes(shape = QC),
    width = 0.2, size = 1.5, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  scale_y_continuous(limits = c(0,max(pantissue$NRF)+0.05)) + 
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "grey", size = 0.5) +
  scale_color_manual(values = col, guide = "none") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(y = "NRF") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "right"  ) 
		
#####################################4.PBC1
p4 <- ggplot(pantissue, aes(x = tissue, y = PBC1, fill = NULL, color=tissue)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(
    aes(shape = QC),
    width = 0.2, size = 1.5, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  scale_y_continuous(limits = c(0,max(pantissue$PBC1)+0.05)) + 
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "grey", size = 0.5) +
  scale_color_manual(values = col, guide = "none") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 12)) +
  labs(y = "PBC1") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "right"  )

#####################################5.PBC2
p5 <- ggplot(pantissue, aes(x = tissue, y = PBC2, fill = NULL, color=tissue)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(
    aes(shape = QC),
    width = 0.2, size = 1.5, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  scale_y_continuous(limits = c(0,max(pantissue$PBC2)+0.05)) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey", size = 0.5) + 
  scale_color_manual(values = col, guide = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(y = "PBC2") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "right"  ) 
		
#
combined_plot <- p1 + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
                 p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
                 p4 + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
                 p5 + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
                 p3 + 
                 plot_layout(nrow = 5, heights = c(0.7, 0.7, 0.7, 0.7, 0.7))


#
ggsave("pantissue.GEO_ERR_ENCODE_sample_QC.boxplot.pdf",combined_plot,width=9.5,height=10) ##Figure4A
