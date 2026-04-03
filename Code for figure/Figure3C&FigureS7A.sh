##########################################
############Figure3C&FigureS7A############
##########################################
###################################
##########1.Sei prediction#########[GPU]
###################################
cd /home/data/t080529/software/sei-framework
################################caSNPs in_corresponding_peak
bash 1_variant_effect_prediction.sh \
/home/data/t080529/RWAS/Sei/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.vcf \
hg19 \
/home/data/t080529/RWAS/Sei/output/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.variant_effect_prediction \
--cuda
bash 2_varianteffect_sc_score.sh \
/home/data/t080529/RWAS/Sei/output/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.variant_effect_prediction/chromatin-profiles-hdf5/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.ref_predictions.h5 \
/home/data/t080529/RWAS/Sei/output/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.variant_effect_prediction/chromatin-profiles-hdf5/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.alt_predictions.h5 \
/home/data/t080529/RWAS/Sei/output/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.variant_effect_prediction
#################################asSNPs
bash 1_variant_effect_prediction.sh \
/home/data/t080529/RWAS/Sei/input/95sample_ASCA_FDR0.1.SNP.vcf \
hg19 \
/home/data/t080529/RWAS/Sei/output/95sample_ASCA_FDR0.1.SNP.variant_effect_prediction \
--cuda
bash 2_varianteffect_sc_score.sh \
/home/data/t080529/RWAS/Sei/output/95sample_ASCA_FDR0.1.SNP.variant_effect_prediction/chromatin-profiles-hdf5/95sample_ASCA_FDR0.1.SNP.ref_predictions.h5 \
/home/data/t080529/RWAS/Sei/output/95sample_ASCA_FDR0.1.SNP.variant_effect_prediction/chromatin-profiles-hdf5/95sample_ASCA_FDR0.1.SNP.alt_predictions.h5 \
/home/data/t080529/RWAS/Sei/output/95sample_ASCA_FDR0.1.SNP.variant_effect_prediction 

####################################
#########2.concordance score########
####################################
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/DeepSEA/Sei_result_enrichment
######
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
#############################################################################
###################i.caSNP in_corresponding_peak prediction####################
#############################################################################
caQTL_inpeak_predict=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/DeepSEA/output/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.variant_effect_prediction/sorted.95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.chromatin_profile_diffs.tsv")
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt")
caQTL_inpeak=subset(caQTL,pheno_var_dist==0)
caQTL_inpeak=caQTL_inpeak[,c("var_id","beta")]
caQTL_inpeak_predict=caQTL_inpeak_predict %>%
               left_join(caQTL_inpeak,by=c("name"="var_id"))


########
# 
cols_to_remove <- grep("DNase\\.fdr|DNase\\.all|DNase\\.hot", colnames(caQTL_inpeak_predict))
# 
caQTL_inpeak_predict <- caQTL_inpeak_predict[, -cols_to_remove, with = FALSE]


###
# Define the range of columns to check  
start_col <- 11  
end_col <- 21761  
# Iterate over each row and modify the specified columns  
caQTL_inpeak_predict[, (start_col:end_col) := lapply(.SD, function(x) {  
  ifelse((beta < 0 & x < 0) | (beta >= 0 & x >= 0), abs(x), -abs(x))  
}), .SDcols = start_col:end_col]  

################################################
# 
mean_result <- data.table(column = character(), mean_value = numeric())  
#
for (col in start_col:end_col) {  
  mean_value <- mean(caQTL_inpeak_predict[[col]], na.rm = TRUE)  
  mean_result <- rbind(mean_result, data.table(column = colnames(caQTL_inpeak_predict)[col], mean_value = mean_value))  
} 
#
overall_mean <- mean(mean_result$mean_value, na.rm = TRUE)  
overall_sd <- sd(mean_result$mean_value, na.rm = TRUE)  
# z-score & p
mean_result[, `:=`(  
  z_score = (mean_value - overall_mean) / overall_sd, 
  #p_value = 2 * (1 - pnorm(abs((mean_value - overall_mean) / overall_sd))) 
  log_p_value = pnorm(abs((mean_value - overall_mean) / overall_sd), lower.tail = FALSE, log.p = TRUE) + log(2) 
)]  
mean_result$p_value = exp(mean_result$log_p_value)  
fwrite(mean_result,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/DeepSEA/Sei_result_enrichment/95sample_caQTL.in_corresponding_peak.6039caSNP.Sei_prediction_consist_result.csv")

#########################################################################
############################ii.asSNP prediction##########################
#########################################################################
ASCA_predict=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/DeepSEA/output/95sample_ASCA_FDR0.1.SNP.variant_effect_prediction/sorted.95sample_ASCA_FDR0.1.SNP.chromatin_profile_diffs.tsv")
ASCA=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1")
ASCA=ASCA[,c("RSID","C0.AF")]
ASCA_predict=ASCA_predict %>%
               left_join(ASCA,by=c("name"="RSID"))

########
# 
cols_to_remove <- grep("DNase\\.fdr|DNase\\.all|DNase\\.hot", colnames(ASCA_predict))
#
ASCA_predict <- ASCA_predict[, -cols_to_remove, with = FALSE]

###
ASCA_predict$C0.AF.direction=ifelse(ASCA_predict$C0.AF<0.5,-1,1)
# Define the range of columns to check  
start_col <- 11  
end_col <- 21761  
# Iterate over each row and modify the specified columns  
ASCA_predict[, (start_col:end_col) := lapply(.SD, function(x) {  
  ifelse((C0.AF.direction < 0 & x < 0) | (C0.AF.direction >= 0 & x >= 0), abs(x), -abs(x))  
}), .SDcols = start_col:end_col]  

##################################################
# 
mean_result <- data.table(column = character(), mean_value = numeric())  
# 
for (col in start_col:end_col) {  
  mean_value <- mean(ASCA_predict[[col]], na.rm = TRUE)  
  mean_result <- rbind(mean_result, data.table(column = colnames(ASCA_predict)[col], mean_value = mean_value))  
} 
#
overall_mean <- mean(mean_result$mean_value, na.rm = TRUE)  
overall_sd <- sd(mean_result$mean_value, na.rm = TRUE)  
# z-score & p
mean_result[, `:=`(  
  z_score = (mean_value - overall_mean) / overall_sd, 
  #p_value = 2 * (1 - pnorm(abs((mean_value - overall_mean) / overall_sd))) 
  log_p_value = pnorm(abs((mean_value - overall_mean) / overall_sd), lower.tail = FALSE, log.p = TRUE) + log(2) 
)]  
mean_result$p_value = exp(mean_result$log_p_value) 
fwrite(mean_result,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/DeepSEA/Sei_result_enrichment/95sample_ASCA.10470asSNP.Sei_prediction_consist_result.csv")
q()

#########################
#########4.QQplot########
##########################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(qqman)
library(dplyr)
#################################################caSNP in OCR
setwd("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/DeepSEA/Sei_result_enrichment_QQplot")
mean_result=fread("../Sei_result_enrichment/95sample_caQTL.in_corresponding_peak.6039caSNP.Sei_prediction_consist_result.csv")
mean_result = mean_result[order(mean_result$p_value,decreasing = F),]
mean_result$obs = -log10(mean_result$p_value)
mean_result$exp = -log10(ppoints(length(mean_result$p_value)))
#
mean_result <- mean_result %>%  
  mutate(  
    color = case_when(  
      grepl("Gastric|Stomach", column) ~ "#E31A1C",  
      TRUE ~ "black"  
    ),  
    legend_label = case_when(  
      grepl("Gastric|Stomach", column) ~ column,  
      TRUE ~ NA_character_  
    )  
  )  

#
gastric_colors <- c("#E31A1C", "#FC4E2A", "#EC7014", "#FE9929", "#FEC44F")

# 
non_na_indices <- which(!is.na(mean_result$legend_label))
# 
num_colors <- length(gastric_colors)
indices_to_color <- head(non_na_indices, num_colors)

#
mean_result$color[indices_to_color] <- gastric_colors[1:length(indices_to_color)]

#
top10 <- mean_result %>%  
  arrange(desc(obs)) %>%  
  slice(1:10)  
##
library(ggplot2)

#
legend_data <- top10[!is.na(top10$legend_label), c("color", "legend_label")]
legend_data <- legend_data[!duplicated(legend_data$legend_label), ]

##
setwd("/data1/gy/ATAC_for_review/Figure3C&FigureS7A/output")
p1 <- ggplot(mean_result, aes(x = exp, y = obs)) +  
  geom_point(data = mean_result, aes(color = "black"), size = 3.5, shape = 16) +  
  geom_point(data = top10, aes(color = color), size = 3.5, shape = 16) +  
  geom_abline(intercept = 0, slope = 1) +  
  scale_x_continuous(expand = c(0.01, 0.01)) +  
  scale_y_continuous(breaks = seq(0, 15, by = 2.5)) +  
  labs(x = expression(paste("Expected -log"[10], "(", italic('P'), ")")),  
       y = expression(paste("Observed -log"[10], "(", italic('P'), ")"))) +  
  theme_bw() +  
  theme(panel.grid = element_blank(),  
        panel.border = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 12, color = 'black', face = 'plain'),  
        axis.title = element_text(size = 12, color = 'black', face = 'plain'),  
        legend.text = element_text(size = 12),  
        legend.title = element_blank(),  
        legend.position = c(0.8, 0.13)) +  
  scale_color_identity(
    guide = "legend",
    breaks = legend_data$color,       
    labels = legend_data$legend_label
  )
ggsave(p1,file="95sample_caQTL.in_corresponding_peak.6039caSNP.Sei_prediction_consist_result.QQplot.top10anno.pdf") ##Figure3C

#################################################asSNP
setwd("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/DeepSEA/Sei_result_enrichment_QQplot")
mean_result=fread("../Sei_result_enrichment/95sample_ASCA.10470asSNP.Sei_prediction_consist_result.csv")
mean_result = mean_result[order(mean_result$p_value,decreasing = F),]
mean_result$obs = -log10(mean_result$p_value)
mean_result$exp = -log10(ppoints(length(mean_result$p_value)))
#
mean_result <- mean_result %>%  
  mutate(  
    color = case_when(  
      grepl("Gastric|Stomach", column) ~ "#E31A1C",  
      TRUE ~ "black"  
    ),  
    legend_label = case_when(  
      grepl("Gastric|Stomach", column) ~ column,  
      TRUE ~ NA_character_  
    )  
  )  

#
gastric_colors <- c("#E31A1C", "#FC4E2A", "#EC7014", "#FE9929", "#FEC44F")

#
non_na_indices <- which(!is.na(mean_result$legend_label))
#
num_colors <- length(gastric_colors)
indices_to_color <- head(non_na_indices, num_colors)

#
mean_result$color[indices_to_color] <- gastric_colors[1:length(indices_to_color)]

#
top10 <- mean_result %>%  
  arrange(desc(obs)) %>%  
  slice(1:10)  
##
library(ggplot2)

#
legend_data <- top10[!is.na(top10$legend_label), c("color", "legend_label")]
legend_data <- legend_data[!duplicated(legend_data$legend_label), ]

##
setwd("/data1/gy/ATAC_for_review/Figure3C&FigureS7A/output")
p2 <- ggplot(mean_result, aes(x = exp, y = obs)) +  
  geom_point(data = mean_result, aes(color = "black"), size = 3.5, shape = 16) +  
  geom_point(data = top10, aes(color = color), size = 3.5, shape = 16) +  
  geom_abline(intercept = 0, slope = 1) +  
  scale_x_continuous(expand = c(0.01, 0.01)) +  
  scale_y_continuous(breaks = seq(0, 15, by = 2.5)) +  
  labs(x = expression(paste("Expected -log"[10], "(", italic('P'), ")")),  
       y = expression(paste("Observed -log"[10], "(", italic('P'), ")"))) +  
  theme_bw() +  
  theme(panel.grid = element_blank(),  
        panel.border = element_blank(),  
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 12, color = 'black', face = 'plain'),  
        axis.title = element_text(size = 12, color = 'black', face = 'plain'),  
        legend.text = element_text(size = 12),  
        legend.title = element_blank(),  
        legend.position = c(0.8, 0.13)) +  
  scale_color_identity(
    guide = "legend",
    breaks = legend_data$color,       
    labels = legend_data$legend_label 
  )
ggsave(p2,file="95sample_ASCA.10470asSNP.Sei_prediction_consist_result.QQplot.top10anno.pdf") ##FigureS7A

