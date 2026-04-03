#########################
#########Figure3A########
#########################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
caQTL_sig=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt")
caQTL_sig_bed=unique(caQTL_sig[,c(2,3,4,1)])
caQTL_sig_bed$pheno_start=caQTL_sig_bed$pheno_start-1
write.table(caQTL_sig_bed,"/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.bed",col.names=F,row.names=F,quote=F,sep="\t")
dim(caQTL_sig_bed)
#[1] 13665     4
caQTL_sig_inpeak=subset(caQTL_sig,caQTL_sig$pheno_var_dist==0)
caQTL_sig_inpeak_bed=unique(caQTL_sig_inpeak[,c(2,3,4,1)])
caQTL_sig_inpeak_bed$pheno_start=caQTL_sig_inpeak_bed$pheno_start-1
write.table(caQTL_sig_inpeak_bed,"/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.inpeak.FDR0.1.bed",col.names=F,row.names=F,quote=F,sep="\t")
dim(caQTL_sig_inpeak_bed)
#[1] 3229    4
ASCA_sig=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1")
ASCA_sig_bed=unique(ASCA_sig[,c(1,4,5,6)])
write.table(ASCA_sig_bed,"/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1.bed",col.names=F,row.names=F,quote=F,sep="\t")
dim(ASCA_sig_bed)
#[1] 5658    4

########################barplot
library(ggplot2)  

#
data <- data.frame(  
  Peak = c("caOCR", "caOCR", "asOCR"),  
  Value = c(3229, 13665 - 3229, 5658),  
  Section = c("caOCR with caSNP in corrosponding OCR", "caOCR without caSNP in corrosponding OCR", "asOCR")  
)  

##
data$Peak <- factor(data$Peak, levels = c("caOCR", "asOCR"))  
data$Section <- factor(data$Section, levels = c("caOCR without caSNP in corrosponding OCR", "caOCR with caSNP in corrosponding OCR", "asOCR")) 

#
setwd("/data1/gy/ATAC_for_review/Figure3A/output")
p=ggplot(data, aes(x = Peak, y = Value, fill = Section)) +  
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c(
  "caOCR with caSNP in corrosponding OCR" = "#08519C",
  "caOCR without caSNP in corrosponding OCR" = "#6BAED6",
  "asOCR" = "#238B45"
  )) + 
  labs(title = "", x = "", y = "Peak") + 
  theme_minimal() +
  theme_minimal() +
  theme(  
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black")
   )+  
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))
ggsave(p,file="significant_caQTL_ASCA_FDR0.1.OCR_barplot.pdf",height=4,width=5.5) ##Figure3A (left)


# venn plot
library(ggvenn)
y= list('caOCR with caSNP in corrosponding OCR'= unique(caQTL_sig_inpeak_bed$pheno_id),'asOCR'= unique(ASCA_sig_bed$NAME))
venn=ggvenn(
  y, 
  show_percentage = F,
  fill_color = c("#08519C", "#238B45"),
  stroke_size = 0.5, set_name_size = 4,
  auto_scale = T
  )
ggsave(venn,file="significant_caQTL_ASCA_FDR0.1.inOCR.OCR_vennplot.pdf",height=4,width=7) ##Figure3A (right)
