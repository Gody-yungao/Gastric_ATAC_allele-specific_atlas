###############################################
##################FigureS6G-H##################
###############################################
########################################################
#######1.vSampler sample con-caSNPs as control set######
########################################################
##http://www.mulinlab.org/vsampler/
#Query Format VCF-like 
#Population 1kgEAS
#Minor Allele Frequency Deviation ±0.05
#Number of Variants in LD Deviation ±50 using r2 0.1
#Match Coding/Noncoding Region
#Exclude Input SNPs
#Match Variant Type
#Sampling Number 1
#Annotation Number 1
#Random Seed 1

################################
#######2.input SNP bed文件######
################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(tidyr)
##
options(scipen=999) 
##
caSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_caQTL_FDR0.1.SNP.vcf")
caSNP$start=caSNP$POS-1
caSNP=caSNP[,c(1,6,2,3)]
##
caSNP_vSampler=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/vSampler_control/95sample_caQTL_FDR0.1.SNP.without_chr.vSampler_control.vcf")
caSNP_vSampler$start=caSNP_vSampler$POS-1
caSNP_vSampler=caSNP_vSampler[,c(1,6,2,3)]
caSNP_vSampler$'#CHROM'=paste0("chr",caSNP_vSampler$'#CHROM')
##
caSNP = caSNP[order(caSNP$`#CHROM`, as.numeric(caSNP$POS)), ]  
caSNP_vSampler = caSNP_vSampler[order(caSNP_vSampler$`#CHROM`, as.numeric(caSNP_vSampler$POS)), ]  
##
fwrite(caSNP,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/chromatin_feature_anno/input/95sample_caQTL_FDR0.1.SNP.bed",row.names=F,col.names=F,quote=F,sep="\t")
fwrite(caSNP_vSampler,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/chromatin_feature_anno/input/95sample_caQTL_FDR0.1.SNP.vSampler_control.bed",row.names=F,col.names=F,quote=F,sep="\t")

###################################################################
##################3 enrichment of caSNP in NJ-GaEpi##################
###################################################################
#############################################################################
#######i.bedtools counting overlaps between SNPs and different features######
#############################################################################
conda activate /Public/gaoyun/miniconda3/envs/ATACseq
cd /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/chromatin_feature_anno/input/
####################################caSNP
bedtools intersect \
-a 95sample_caQTL_FDR0.1.SNP.bed \
-b \
/data1/gy/fresh_tissue_epigeno/H3K4me1/normal_sample_merged_macs2/stomach_H3K4me1_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/H3K4me3/normal_sample_merged_macs2/stomach_H3K4me3_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/H3K27ac/normal_sample_merged_macs2/stomach_H3K27ac_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/H3K27me3/normal_sample_merged_macs2/stomach_H3K27me3_merged_with_peakname.broadPeak \
/data1/gy/fresh_tissue_epigeno/H3K9me3/normal_sample_merged_macs2/stomach_H3K9me3_merged_with_peakname.broadPeak \
/data1/gy/fresh_tissue_epigeno/RNAPOLII/normal_sample_merged_macs2/stomach_RNAPOLII_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/BRD4/normal_sample_merged_macs2/stomach_BRD4_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/CTCF/macs2_peak/stomach_CTCF_merged_with_peakname.narrowPeak \
-C -filenames \
> ../NJMU_stomach_output/95sample_caQTL_FDR0.1.SNP.intersected_with_NJMU_stomach_feature
#
awk -F'\t' '{  
    split($5, a, "/"); 
    #
    feature = a[5];  
    #
    print $1, $2, $3, $4, feature, $6;  
}' OFS="\t" ../NJMU_stomach_output/95sample_caQTL_FDR0.1.SNP.intersected_with_NJMU_stomach_feature \
> ../NJMU_stomach_output/95sample_caQTL_FDR0.1.SNP.intersected_with_NJMU_stomach_feature.final
####################################caSNP_vSampler
bedtools intersect \
-a 95sample_caQTL_FDR0.1.SNP.vSampler_control.bed \
-b \
/data1/gy/fresh_tissue_epigeno/H3K4me1/normal_sample_merged_macs2/stomach_H3K4me1_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/H3K4me3/normal_sample_merged_macs2/stomach_H3K4me3_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/H3K27ac/normal_sample_merged_macs2/stomach_H3K27ac_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/H3K27me3/normal_sample_merged_macs2/stomach_H3K27me3_merged_with_peakname.broadPeak \
/data1/gy/fresh_tissue_epigeno/H3K9me3/normal_sample_merged_macs2/stomach_H3K9me3_merged_with_peakname.broadPeak \
/data1/gy/fresh_tissue_epigeno/RNAPOLII/normal_sample_merged_macs2/stomach_RNAPOLII_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/BRD4/normal_sample_merged_macs2/stomach_BRD4_merged_with_peakname.narrowPeak \
/data1/gy/fresh_tissue_epigeno/CTCF/macs2_peak/stomach_CTCF_merged_with_peakname.narrowPeak \
-C -filenames \
> ../NJMU_stomach_output/95sample_caQTL_FDR0.1.SNP.vSampler_control.intersected_with_NJMU_stomach_feature
#
awk -F'\t' '{  
    split($5, a, "/");
    feature = a[5];
    print $1, $2, $3, $4, feature, $6;  
}' OFS="\t" ../NJMU_stomach_output/95sample_caQTL_FDR0.1.SNP.vSampler_control.intersected_with_NJMU_stomach_feature \
> ../NJMU_stomach_output/95sample_caQTL_FDR0.1.SNP.vSampler_control.intersected_with_NJMU_stomach_feature.final

######################################
#######ii.long type to wide type######
######################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)  
library(dplyr)  
library(tidyr)  
############################caSNP
caSNP <- fread("../NJMU_stomach_output/95sample_caQTL_FDR0.1.SNP.intersected_with_NJMU_stomach_feature.final", header = F)  
colnames(caSNP) <- c("CHROM", "START", "END", "ID", "FEATURE", "VALUE")  
#
caSNP_wide <- caSNP %>%  
  pivot_wider(names_from = FEATURE, values_from = VALUE, values_fill = list(VALUE = 0))  
############################caSNP_vSampler
caSNP_vSampler <- fread("../NJMU_stomach_output/95sample_caQTL_FDR0.1.SNP.vSampler_control.intersected_with_NJMU_stomach_feature.final", header = F)  
colnames(caSNP_vSampler) <- c("CHROM", "START", "END", "ID", "FEATURE", "VALUE")  
#
caSNP_vSampler_wide <- caSNP_vSampler %>%  
  pivot_wider(names_from = FEATURE, values_from = VALUE, values_fill = list(VALUE = 0))  

#######################################
#######iii.enrichment fisher test######
#######################################
#
results <- data.frame(FEATURE = character(),  
                      Odds_Ratio = numeric(),  
                      P_Value = numeric(),  
                      CI_Lower = numeric(),  
                      CI_Upper = numeric(),  
                      stringsAsFactors = FALSE)  

# 
features <- colnames(caSNP_wide)[5:ncol(caSNP_wide)]

#
for (feature in features) {  
  #
  a <- sum(caSNP_wide[[feature]] == 1)          
  b <- sum(caSNP_wide[[feature]] == 0)          
  c <- sum(caSNP_vSampler_wide[[feature]] == 1)     
  d <- sum(caSNP_vSampler_wide[[feature]] == 0)    

  # 2x2 
  contingency_table <- matrix(c(a, b, c, d), nrow = 2)  

  # fisher test
  fisher_test <- fisher.test(contingency_table, alternative = "two.sided", conf.int = TRUE)  

  #
  results <- rbind(results, data.frame(FEATURE = feature,  
                                        Odds_Ratio = fisher_test$estimate,  
                                        P_Value = fisher_test$p.value,  
                                        CI_Lower = fisher_test$conf.int[1],  
                                        CI_Upper = fisher_test$conf.int[2]))  
}  
#
results <- results %>% 
  mutate(  
      FEATURE = factor(FEATURE, levels = FEATURE[order(Odds_Ratio, decreasing = TRUE)]), 
      bar_start = ifelse(Odds_Ratio >= 1, 1, Odds_Ratio),  
      bar_end = ifelse(Odds_Ratio>= 1,Odds_Ratio, 1),  
      Direction = ifelse(Odds_Ratio >= 1, "Above", "Below"),
      p_label = dplyr::case_when(
      P_Value < 1e-3 ~ "***",
      P_Value < 1e-2 ~ "**",
      P_Value < 5e-2 ~ "*",
      TRUE           ~ "ns")
  )  

#############barplot
library(ggplot2)
pad <- diff(range(results$CI_Lower, results$CI_Upper)) * 0.03

results <- results %>%
  mutate(
    #
    lbl_y = if_else(Odds_Ratio >= 1, bar_end + pad, bar_start - pad),   
    lbl_vjust = if_else(Odds_Ratio >= 1, 0, 1)                        
  )

p1 <- ggplot(results, aes(x = FEATURE)) +
  geom_rect(aes(xmin = as.numeric(FEATURE) - 0.35,
                xmax = as.numeric(FEATURE) + 0.35,
                ymin = bar_start, ymax = bar_end, fill = Direction)) +
  geom_errorbar(aes(x = as.numeric(FEATURE), ymin = CI_Lower, ymax = CI_Upper),
                width = 0.1, color = "black") +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  geom_text(aes(x = as.numeric(FEATURE), y = lbl_y, label = p_label, vjust = lbl_vjust),
            size = 4, color = "black") +
  scale_x_discrete() +
  scale_fill_manual(values = c("Above" = "#E31A1C", "Below" = "#225EA8")) +
  scale_y_continuous(
  breaks = seq(0.5, max(results$CI_Upper)+1, by = 0.5),
  limits = c(min(results$CI_Lower)*0.88, max(results$CI_Upper)+1)) + 
  theme_minimal() +  
  theme(  
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 11),  
    axis.text.y = element_text(size = 11),  
    axis.title.y = element_text(size = 12),  
    axis.title.x = element_text(size = 12),  
    plot.title = element_text(size = 14, hjust = 0.5),  
    legend.position = "none"
  ) +  
  labs(title = "", x = "Chromatin feature (from NJ-GaEpi stomach)", y = "OR (compared with non-caSNPs)")
  #coord_flip() 

#
setwd("/data1/gy/ATAC_for_review/FigureS6G-H/output") 
ggsave(p1, file="95sample_caQTL_FDR0.1.SNP.NJ-GaEpi_feature.enrichment.barplot.pdf", width = 7, height = 6)  ##FigureS6G


###################################################################
##################4 enrichment of caSNP in ENCODE##################
###################################################################
#############################################################################
#######i.bedtools counting overlaps between SNPs and different features######
#############################################################################
conda activate /Public/gaoyun/miniconda3/envs/ATACseq
cd /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/chromatin_feature_anno/input
####################################caSNP
bedtools intersect \
  -a 95sample_caQTL_FDR0.1.SNP.bed \
  -b /data1/gy/public/stomach_ENCODE_region/merged/*.bed \
  -C -filenames \
> ../ENCODE_stomach_output/95sample_caQTL_FDR0.1.SNP.intersected_with_ENCODE_stomach_feature 
#
awk 'BEGIN{OFS="\t"} \
{split($5, a, "/"); split(a[length(a)], b, "."); $5 = b[1]; print}' \
../ENCODE_stomach_output/95sample_caQTL_FDR0.1.SNP.intersected_with_ENCODE_stomach_feature \
> ../ENCODE_stomach_output/95sample_caQTL_FDR0.1.SNP.intersected_with_ENCODE_stomach_feature.final
####################################caSNP_vSampler
bedtools intersect \
  -a 95sample_caQTL_FDR0.1.SNP.vSampler_control.bed \
  -b /data1/gy/public/stomach_ENCODE_region/merged/*.bed \
  -C -filenames \
> ../ENCODE_stomach_output/95sample_caQTL_FDR0.1.SNP.vSampler_control.intersected_with_ENCODE_stomach_feature
#
awk 'BEGIN{OFS="\t"} \
{split($5, a, "/"); split(a[length(a)], b, "."); $5 = b[1]; print}' \
../ENCODE_stomach_output/95sample_caQTL_FDR0.1.SNP.vSampler_control.intersected_with_ENCODE_stomach_feature \
> ../ENCODE_stomach_output/95sample_caQTL_FDR0.1.SNP.vSampler_control.intersected_with_ENCODE_stomach_feature.final

######################################
#######ii.long type to wide type######
######################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)  
library(dplyr)  
library(tidyr)  
############################caSNP
caSNP <- fread("../ENCODE_stomach_output/95sample_caQTL_FDR0.1.SNP.intersected_with_ENCODE_stomach_feature.final", header = F)  
colnames(caSNP) <- c("CHROM", "START", "END", "ID", "FEATURE", "VALUE")  
#
caSNP_wide <- caSNP %>%  
  pivot_wider(names_from = FEATURE, values_from = VALUE, values_fill = list(VALUE = 0))  
############################caSNP_vSampler
caSNP_vSampler <- fread("../ENCODE_stomach_output/95sample_caQTL_FDR0.1.SNP.vSampler_control.intersected_with_ENCODE_stomach_feature.final", header = F)  
colnames(caSNP_vSampler) <- c("CHROM", "START", "END", "ID", "FEATURE", "VALUE")  
# 将第五列拆分并转为宽格式  
caSNP_vSampler_wide <- caSNP_vSampler %>%  
  pivot_wider(names_from = FEATURE, values_from = VALUE, values_fill = list(VALUE = 0))  

#######################################
#######iii.enrichment fisher test######
#######################################
#
results <- data.frame(FEATURE = character(),  
                      Odds_Ratio = numeric(),  
                      P_Value = numeric(),  
                      CI_Lower = numeric(),  
                      CI_Upper = numeric(),  
                      stringsAsFactors = FALSE)  

#
features <- colnames(caSNP_wide)[5:ncol(caSNP_wide)] 

#
for (feature in features) {  
  #
  a <- sum(caSNP_wide[[feature]] == 1)      
  b <- sum(caSNP_wide[[feature]] == 0)  
  c <- sum(caSNP_vSampler_wide[[feature]] == 1)
  d <- sum(caSNP_vSampler_wide[[feature]] == 0)

  # 2x2 
  contingency_table <- matrix(c(a, b, c, d), nrow = 2)  

  # fisher test
  fisher_test <- fisher.test(contingency_table, alternative = "two.sided", conf.int = TRUE)  

  #
  results <- rbind(results, data.frame(FEATURE = feature,  
                                        Odds_Ratio = fisher_test$estimate,  
                                        P_Value = fisher_test$p.value,  
                                        CI_Lower = fisher_test$conf.int[1],  
                                        CI_Upper = fisher_test$conf.int[2]))  
}  

#
results <- results %>% 
  mutate(  
      FEATURE = factor(FEATURE, levels = FEATURE[order(Odds_Ratio, decreasing = TRUE)]), 
      bar_start = ifelse(Odds_Ratio >= 1, 1, Odds_Ratio),  
      bar_end = ifelse(Odds_Ratio>= 1,Odds_Ratio, 1),  
      Direction = ifelse(Odds_Ratio >= 1, "Above", "Below"),
      p_label = dplyr::case_when(
      P_Value < 1e-3 ~ "***",
      P_Value < 1e-2 ~ "**",
      P_Value < 5e-2 ~ "*",
      TRUE           ~ "ns") 
  )  

#####barplot
library(ggplot2)
#
pad <- diff(range(results$CI_Lower, results$CI_Upper)) * 0.03

results <- results %>%
  mutate(
    #
    lbl_y = if_else(Odds_Ratio >= 1, bar_end + pad, bar_start - pad),
    lbl_vjust = if_else(Odds_Ratio >= 1, 0, 1) 
  )

p2 <- ggplot(results, aes(x = FEATURE)) +
  geom_rect(aes(xmin = as.numeric(FEATURE) - 0.35,
                xmax = as.numeric(FEATURE) + 0.35,
                ymin = bar_start, ymax = bar_end, fill = Direction)) +
  geom_errorbar(aes(x = as.numeric(FEATURE), ymin = CI_Lower, ymax = CI_Upper),
                width = 0.1, color = "black") +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  geom_text(aes(x = as.numeric(FEATURE), y = lbl_y, label = p_label, vjust = lbl_vjust),
            size = 4, color = "black") +
  scale_x_discrete() + 
  scale_fill_manual(values = c("Above" = "#E31A1C", "Below" = "#225EA8")) +
  scale_y_continuous(
  breaks = seq(0.5, max(results$CI_Upper)+1, by = 0.5),
  limits = c(min(results$CI_Lower)*0.9, max(results$CI_Upper)+1)) +
  theme_minimal() +  
  theme(  
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"), 
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 11),  
    axis.text.y = element_text(size = 11),  
    axis.title.y = element_text(size = 12),  
    axis.title.x = element_text(size = 12),  
    plot.title = element_text(size = 14, hjust = 0.5),  
    legend.position = "none" 
  ) +  
  labs(title = "", x = "Chromatin feature (from ENCODE stomach)", y = "OR (compared with non-caSNPs)")

#
setwd("/data1/gy/ATAC_for_review/FigureS6G-H/output") 
ggsave(p2, file="95sample_caQTL_FDR0.1.SNP.ENCODE_feature.enrichment.barplot.pdf", width = 10, height = 6)  ##FigureS6H
