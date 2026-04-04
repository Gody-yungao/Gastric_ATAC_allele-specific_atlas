##########################
#########Figure4E#########
##########################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
########################################################################
########1.epithelium-specific OCR caQTL-epi-sceQTL coloc result#########
########################################################################
OCR_annotated=fread("/data1/gy/ATACseq_RWAS/ATACseq/OCR_celltype_anno/macs2_type_peak_anno/95sample_IterativeOverlapPeakSet.type_anno.txt")
#
OCR_single_type <- OCR_annotated[  
  type_specific != "" &   
  type_specific != FALSE
]  
Type <- "Epithelium"
OCR_single_type_sub <- subset(OCR_single_type, type == Type)  
caOCR=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.bed")
OCR_single_type_sub_caOCR=subset(OCR_single_type_sub,OCR %in% caOCR$V4)
##epithelium-specific caOCR num
length(OCR_single_type_sub_caOCR$OCR)
#[1] 2825

####epithelium eGene num
sceQTL_epi=fread("/data1/gy/sceQTL/eQTL_CHN100K/epi/epi/epi_sceQTL_result.allpairs.withsymbol.withFDR.FDR0.05.txt")
length(unique(sceQTL_epi$gene_id))
#[1] 2959

##
Type1="epi"
coloc_path <- paste0(  
    "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/",  
    Type1, "_caQTL_sceQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv"  
  )  
coloc_epi=fread(coloc_path)
coloc_epi=coloc_epi %>%  
  filter(!(is.na(Distance_sceQTL) & is.na(Distance_eQTL)))

###################statistics
both <- coloc_epi %>%  
  filter(!is.na(Distance_sceQTL) & !is.na(Distance_eQTL)) %>%  
  nrow()  
#
sceQTL_coloc_only <- coloc_epi %>%  
  filter(!is.na(Distance_sceQTL) & is.na(Distance_eQTL)) %>%  
  nrow()  
#
stats_matrix <- matrix(  
  c(both, sceQTL_coloc_only),  
  ncol = 2,  
  dimnames = list("Count", c("Both", "sceQTL_coloc_only"))  
)  
stats_matrix=as.data.frame(stats_matrix)
stats_matrix
#      Both sceQTL_coloc_only
#Count  102                47


###################epithelium-specific OCR caQTL-epi-sceQTL coloc result caOCR & eGene num
coloc_epi_sceQTL=na.omit(coloc_epi[,1:19])
length(unique(coloc_epi_sceQTL$Peak))
#[1] 130
length(unique(coloc_epi_sceQTL$Geneid))
#[1] 113

##################Adobe Illustrator plot vennplot Figure4E (top)


#################epithelium-specific coloc with/without bulk coloc anno
sceQTL_coloc_only <- coloc_epi %>%  
  filter(!is.na(Distance_sceQTL) & is.na(Distance_eQTL))
both <- coloc_epi %>%  
  filter(!is.na(Distance_sceQTL) & !is.na(Distance_eQTL))
fwrite(sceQTL_coloc_only, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/epi_caQTL_sceQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.without_bulkcoloc.chr1_22.csv", 
       col.names = TRUE, row.names = FALSE, quote = FALSE)
fwrite(both, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/epi_caQTL_sceQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.with_bulkcoloc.chr1_22.csv", 
       col.names = TRUE, row.names = FALSE, quote = FALSE)

##################################################################################################
########2.pieplot of proportion of epithelium-specific coloc with/without bulk coloc anno#########
##################################################################################################
setwd("/data1/gy/ATAC_for_review/Figure4E/output")
library(ggplot2)  
##
plot_data <- data.frame(  
  Category = factor(c("non-colocalized", "colocalized"),  
                   levels = c("non-colocalized", "colocalized")),  
  Count = c(stats_matrix$sceQTL_coloc_only[1], stats_matrix$Both[1])  
) %>%   
  mutate(  
    Percent = Count/sum(Count)*100,  
    Label = sprintf("%d\n(%.2f%%)", Count, Percent)
  )  
##
p=ggplot(plot_data, aes(x = "", y = Count, fill = Category)) +  
  geom_bar(width = 1, stat = "identity", color = "white") +  
  coord_polar(theta = "y") +  
  geom_text(aes(label = Label),  
            position = position_stack(vjust = 0.5),  
            color = "black",   
            size = 5) +  
  scale_fill_manual(  
    values = c("#D6604D", "#F4A582"),  
    labels = c("non-colocalized bulk caOCR-eGene", "colocalized bulk caOCR-eGene")
  ) +  
  theme_void() +  
  theme(  
    legend.position = "bottom",  
    legend.title = element_blank(),  
    legend.text = element_text(size = 11, margin = margin(r = 15)),  
    legend.key.size = unit(0.8, "cm")  
  ) +  
  guides(fill = guide_legend(nrow = 1))
ggsave(p,filename="epi.caQTL_sceQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.with_without_bulk_coloc.pieplot.pdf") ##Figure4E (bottom)
