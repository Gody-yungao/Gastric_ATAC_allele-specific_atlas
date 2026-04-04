########################################
###########Figure3E&FigureS7B###########
########################################
#################
###1.SuRE data###
#################
##https://osf.io/pjxm4/files/osfstorage 
##https://www.nature.com/articles/nbt.3754
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
#########################SuRE
df=fread("/data1/gy/public/SURE/SuRE_SNP_table_LP190708.txt")
df$chrbp=paste0(df$chr,":",df$SNPabspos)
df=df[,c(2,6:8,10:12,14:16)]
df$hepg2.SuRE_activity=df$hepg2.alt.mean-df$hepg2.ref.mean

#######################################################
###2.1 caSNP(in corresponding peak) in SuRE database###
#######################################################
##########################caSNP (in corresponding peak)
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt")
caQTL=caQTL %>% filter(pheno_var_dist == 0)
caQTL=caQTL[,c(1,8,12,13)]
caQTL <- caQTL %>%  
  separate(var_id, into = c("chr_caSNP", "pos_caSNP", "ref_caSNP", "alt_caSNP"), sep = ":", remove = F)
caQTL$chrbp=paste0(caQTL$chr_caSNP,":",caQTL$pos_caSNP) 

##
caQTL=inner_join(caQTL,df,by=c("chrbp"="chrbp"))
dim(caQTL)
#[1] 4711   20
caQTL <- caQTL %>%  
  mutate(hepg2.SuRE_activity_new = case_when(  
    ref_caSNP == ref & alt_caSNP == alt ~ hepg2.SuRE_activity,   
    ref_caSNP == alt & alt_caSNP == ref ~ -hepg2.SuRE_activity, 
    TRUE ~ NA_real_     
  ))
##
any(is.na(caQTL)) 
#[1] FALSE
length(unique(caQTL$var_id))
#[1] 4711
write.csv(caQTL, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/caSNP_in_corresponing_peak_SuRE_cor_intersected.csv", row.names = FALSE)

###########################correlation
#######################HepG2 
#
thresholds <- c(1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7)
#
midpoints <- sqrt(thresholds[-length(thresholds)] * thresholds[-1])
#
new_thresholds <- sort(c(thresholds, midpoints), decreasing = TRUE)
#
results <- data.frame(threshold = integer(),  
                      pair_count = integer(),  
                      correlation_coefficient = numeric(),  
                      p_value = numeric(),  
                      ci_lower = numeric(),  
                      ci_upper = numeric(),  
                      stringsAsFactors = FALSE)  

# 
for (threshold in new_thresholds) {  
  # 
  filtered_data <- caQTL %>%  
    filter(hepg2.wilcox.p.value <= threshold)  
  
  #
  pair_count <- nrow(filtered_data)  
  
  #
  if (pair_count > 1) { 
    cor_test_result <- cor.test(filtered_data$beta,
                                filtered_data$hepg2.SuRE_activity_new,
                                use = "complete.obs")
  
    #
    correlation_coefficient <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
  
    #
    n <- sum(complete.cases(filtered_data$beta,
                            filtered_data$hepg2.SuRE_activity_new))
  
    #
    sem <- (1 - correlation_coefficient^2) / sqrt(n - 2)
  
    #
    ci_lower <- correlation_coefficient - sem
    ci_upper <- correlation_coefficient + sem
  
  } else {
    correlation_coefficient <- NA
    p_value <- NA
    ci_lower <- NA
    ci_upper <- NA
  }
  
  #
  results <- rbind(results, data.frame(threshold = threshold,  
                                        pair_count = pair_count,  
                                        correlation_coefficient = correlation_coefficient,  
                                        p_value = p_value,  
                                        ci_lower = ci_lower,  
                                        ci_upper = ci_upper))  
}  
# 
write.csv(results, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/caSNP_in_corresponing_peak_HepG2_SuRE_activity_correlation_results_by_threshold.csv", row.names = FALSE)

###########################################################
###2.2 non-caSNP(in corresponding peak) in SuRE database###
###########################################################
##########################non-caSNP (in corresponding peak)
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
noncaQTL=caQTL %>% filter(pheno_var_dist == 0)%>% filter(FDR>=0.1)
noncaQTL=noncaQTL[,c(1,8,12,13)]
noncaQTL <- noncaQTL %>%  
  separate(var_id, into = c("chr_noncaSNP", "pos_noncaSNP", "ref_noncaSNP", "alt_noncaSNP"), sep = ":", remove = F)
noncaQTL$chrbp=paste0(noncaQTL$chr_noncaSNP,":",noncaQTL$pos_noncaSNP) 

##
noncaQTL=inner_join(noncaQTL,df,by=c("chrbp"="chrbp"))
dim(noncaQTL)
#[1] 67818    20
noncaQTL <- noncaQTL %>%  
  mutate(hepg2.SuRE_activity_new = case_when(  
    ref_noncaSNP == ref & alt_noncaSNP == alt ~ hepg2.SuRE_activity,  
    ref_noncaSNP == alt & alt_noncaSNP == ref ~ -hepg2.SuRE_activity,
    TRUE ~ NA_real_     
  ))
##
any(is.na(noncaQTL)) 
#[1] FALSE
length(unique(noncaQTL$var_id))
#[1] 67818
write.csv(noncaQTL, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/noncaSNP_in_corresponing_peak_SuRE_cor_intersected.csv", row.names = FALSE)

###########################correlation
##########################HepG2 
#
thresholds <- c(1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7)
# 
midpoints <- sqrt(thresholds[-length(thresholds)] * thresholds[-1])
# 
new_thresholds <- sort(c(thresholds, midpoints), decreasing = TRUE)
# 
results <- data.frame(threshold = integer(),  
                      pair_count = integer(),  
                      correlation_coefficient = numeric(),  
                      p_value = numeric(),  
                      ci_lower = numeric(),  
                      ci_upper = numeric(),  
                      stringsAsFactors = FALSE)  

#
for (threshold in new_thresholds) {  
  #
  filtered_data <- noncaQTL %>%  
    filter(hepg2.wilcox.p.value <= threshold)  
  
  #
  pair_count <- nrow(filtered_data)  
  
  #
  if (pair_count > 1) { 
    cor_test_result <- cor.test(filtered_data$beta,
                                filtered_data$hepg2.SuRE_activity_new,
                                use = "complete.obs")
  
    #
    correlation_coefficient <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
  
    #
    n <- sum(complete.cases(filtered_data$beta,
                            filtered_data$hepg2.SuRE_activity_new))
  
    #
    sem <- (1 - correlation_coefficient^2) / sqrt(n - 2)
  
    #
    ci_lower <- correlation_coefficient - sem
    ci_upper <- correlation_coefficient + sem
  
  } else {
    correlation_coefficient <- NA
    p_value <- NA
    ci_lower <- NA
    ci_upper <- NA
  }
  
  #
  results <- rbind(results, data.frame(threshold = threshold,  
                                        pair_count = pair_count,  
                                        correlation_coefficient = correlation_coefficient,  
                                        p_value = p_value,  
                                        ci_lower = ci_lower,  
                                        ci_upper = ci_upper))  
}  
#
write.csv(results, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/noncaSNP_in_corresponing_peak_HepG2_SuRE_activity_correlation_results_by_threshold.csv", row.names = FALSE)

#######################################################
###3.1.asSNP(in corresponding peak) in SuRE database###
#######################################################
library(purrr)
##########################caSNP
ASCA=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1")
ASCA=ASCA[,c(3,6,10,11)]
ASCA$ALT_REF_ratio_log2 <- log2(ASCA$ALL.AF / (1 - ASCA$ALL.AF)) 
ASCA <- ASCA %>%  
  separate(RSID, into = c("chr_asSNP", "pos_asSNP", "ref_asSNP", "alt_asSNP"), sep = ":", remove = F)
ASCA$chrbp=paste0(ASCA$chr_asSNP,":",ASCA$pos_asSNP) 
##
ASCA=inner_join(ASCA,df,by=c("chrbp"="chrbp"))
ASCA <- ASCA %>%  
  mutate(hepg2.SuRE_activity_new = case_when(  
    ref_asSNP == ref & alt_asSNP == alt ~ hepg2.SuRE_activity,  
    ref_asSNP == alt & alt_asSNP == ref ~ -hepg2.SuRE_activity, 
    TRUE ~ NA_real_  
  ))
##
any(is.na(ASCA)) 
#[1] FALSE
length(unique(ASCA$RSID))
#[1] 8061
write.csv(ASCA, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/asSNP_SuRE_intersected.csv", row.names = FALSE)

###########################correlation
#######################HepG2
#
thresholds <- c(1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7)
#
midpoints <- sqrt(thresholds[-length(thresholds)] * thresholds[-1])
# 
new_thresholds <- sort(c(thresholds, midpoints), decreasing = TRUE)
#
results <- data.frame(threshold = integer(),  
                      pair_count = integer(),  
                      correlation_coefficient = numeric(),  
                      p_value = numeric(),  
                      ci_lower = numeric(),  
                      ci_upper = numeric(),  
                      stringsAsFactors = FALSE)  

# 
for (threshold in new_thresholds) {  
  # 
  filtered_data <- ASCA %>%  
    filter(hepg2.wilcox.p.value <= threshold)  
  
  #
  pair_count <- nrow(filtered_data)  
  
  # 
  if (pair_count > 1) { 
    cor_test_result <- cor.test(filtered_data$ALT_REF_ratio_log2,
                                filtered_data$hepg2.SuRE_activity_new,
                                use = "complete.obs")
  
    #
    correlation_coefficient <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
  
    #
    n <- sum(complete.cases(filtered_data$ALT_REF_ratio_log2,
                            filtered_data$hepg2.SuRE_activity_new))
  
    #
    sem <- (1 - correlation_coefficient^2) / sqrt(n - 2)
  
    #
    ci_lower <- correlation_coefficient - sem
    ci_upper <- correlation_coefficient + sem
  
  } else {
    correlation_coefficient <- NA
    p_value <- NA
    ci_lower <- NA
    ci_upper <- NA
  }
  
  #
  results <- rbind(results, data.frame(threshold = threshold,  
                                        pair_count = pair_count,  
                                        correlation_coefficient = correlation_coefficient,  
                                        p_value = p_value,  
                                        ci_lower = ci_lower,  
                                        ci_upper = ci_upper))  
}  
#
write.csv(results, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/asSNP_HepG2_SuRE_activity_correlation_results_by_threshold.csv", row.names = FALSE)

##########################################################
###3.2.nonasSNP(in corresponding peak) in SuRE database###
##########################################################
library(purrr)
##########################non-asSNP(in corresponding peak)
ASCA=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR")
nonASCA=ASCA %>% filter(C0.BBINOM.FDR>=0.1)
nonASCA=nonASCA[,c(3,6,10,11)]
nonASCA$ALT_REF_ratio_log2 <- log2(nonASCA$ALL.AF / (1 - nonASCA$ALL.AF)) 
nonASCA <- nonASCA %>%  
  separate(RSID, into = c("chr_nonasSNP", "pos_nonasSNP", "ref_nonasSNP", "alt_nonasSNP"), sep = ":", remove = F)
nonASCA$chrbp=paste0(nonASCA$chr_nonasSNP,":",nonASCA$pos_nonasSNP) 
##
nonASCA=inner_join(nonASCA,df,by=c("chrbp"="chrbp"))
nonASCA <- nonASCA %>%  
  mutate(hepg2.SuRE_activity_new = case_when(  
    ref_nonasSNP == ref & alt_nonasSNP == alt ~ hepg2.SuRE_activity,  
    ref_nonasSNP == alt & alt_nonasSNP == ref ~ -hepg2.SuRE_activity,
    TRUE ~ NA_real_ 
  ))
##
any(is.na(nonASCA)) 
#[1] FALSE
length(unique(nonASCA$RSID))
#[1] 64513
write.csv(nonASCA, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/nonasSNP_SuRE_intersected.csv", row.names = FALSE)

###########################cor
#######################HepG2
#
thresholds <- c(1, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7)
#
midpoints <- sqrt(thresholds[-length(thresholds)] * thresholds[-1])
#
new_thresholds <- sort(c(thresholds, midpoints), decreasing = TRUE)
#
results <- data.frame(threshold = integer(),  
                      pair_count = integer(),  
                      correlation_coefficient = numeric(),  
                      p_value = numeric(),  
                      ci_lower = numeric(),  
                      ci_upper = numeric(),  
                      stringsAsFactors = FALSE)  

#
for (threshold in new_thresholds) {  
  # 
  filtered_data <- nonASCA %>%  
    filter(hepg2.wilcox.p.value <= threshold)  
  
  #
  pair_count <- nrow(filtered_data)  
  
  # 
  if (pair_count > 1) { 
    cor_test_result <- cor.test(filtered_data$ALT_REF_ratio_log2,
                                filtered_data$hepg2.SuRE_activity_new,
                                use = "complete.obs")
  
    #
    correlation_coefficient <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
  
    #
    n <- sum(complete.cases(filtered_data$ALT_REF_ratio_log2,
                            filtered_data$hepg2.SuRE_activity_new))
  
    #
    sem <- (1 - correlation_coefficient^2) / sqrt(n - 2)
  
    #
    ci_lower <- correlation_coefficient - sem
    ci_upper <- correlation_coefficient + sem
  
  } else {
    correlation_coefficient <- NA
    p_value <- NA
    ci_lower <- NA
    ci_upper <- NA
  }
  
  # 
  results <- rbind(results, data.frame(threshold = threshold,  
                                        pair_count = pair_count,  
                                        correlation_coefficient = correlation_coefficient,  
                                        p_value = p_value,  
                                        ci_lower = ci_lower,  
                                        ci_upper = ci_upper))  
}  
#
write.csv(results, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/nonasSNP_HepG2_SuRE_activity_correlation_results_by_threshold.csv", row.names = FALSE)

###################
###4.scatterplot###
###################
library(ggplot2)  
setwd("/data1/gy/ATAC_for_review/Figure3E&FigureS7B/output")
################################################
###i.caSNP & non-caSNP (in corresponding OCR)###
################################################
caSNP_result=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/caSNP_in_corresponing_peak_HepG2_SuRE_activity_correlation_results_by_threshold.csv")
caSNP_result$group="caSNP"
noncaSNP_result=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/noncaSNP_in_corresponing_peak_HepG2_SuRE_activity_correlation_results_by_threshold.csv")
noncaSNP_result$group="non-caSNP"
result=rbind(caSNP_result,noncaSNP_result)
result=subset(result,threshold>1e-7)
#
result$threshold <- as.factor(-log10(result$threshold))  
##
p1 <- ggplot(result, aes(x = threshold, 
                         y = correlation_coefficient, 
                         color = group)) +

  geom_point(size = 3, alpha = 0.75) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, alpha = 0.75) + 
  geom_text(aes(label = pair_count), hjust = -0.3, vjust = 2, size = 2.5, color = "black") + 

  labs(title = "",
       x = "Maximal SuRE SNP P value",
       y = "Correlation coefficient\nbetween SuRE SNP ΔExpression (ALT-REF) values\nand caSNP (in corresponding OCR) caQTL effect size") +

  #
  scale_color_manual(values = c("caSNP" = "#094D92", "non-caSNP" = "grey")) +

  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_blank(),
    legend.title = element_blank() 
  )
ggsave(p1,filename="caSNP_noncaSNP_in_corresponing_peak_HepG2_SuRE_activity_correlation_results_by_threshold.Scatterplot.pdf",width=8,height=7) ##Figure3E

#################################################
###ii.asSNP & non-asSNP (in corresponding OCR)###
#################################################
asSNP_result=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/asSNP_HepG2_SuRE_activity_correlation_results_by_threshold.csv")
asSNP_result$group="asSNP"
nonasSNP_result=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SuRE_cor/nonasSNP_HepG2_SuRE_activity_correlation_results_by_threshold.csv")
nonasSNP_result$group="non-asSNP"
result=rbind(asSNP_result,nonasSNP_result)
result=subset(result,threshold>1e-7)
#
result$threshold <- as.factor(-log10(result$threshold))  
##
p2 <- ggplot(result, aes(x = threshold, 
                         y = correlation_coefficient, 
                         color = group)) + 

  geom_point(size = 3, alpha = 0.75) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, alpha = 0.75) + 
  geom_text(aes(label = pair_count), hjust = -0.3, vjust = 2, size = 2.5, color = "black") +

  labs(title = "",
       x = "Maximal SuRE SNP P value",
       y = "Correlation coefficient\nbetween SuRE SNP ΔExpression(ALT-REF) values\nand asSNP ASCA log2(ALT/REF)") +

  #
  scale_color_manual(values = c("asSNP" = "#238444", "non-asSNP" = "#7A7A7A")) +

  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_blank(),
    legend.title = element_blank() 
  )
ggsave(p2,filename="asSNP_nonasSNP_HepG2_SuRE_activity_correlation_results_by_threshold.Scatterplot.pdf",width=8,height=7) ##FigureS7B

