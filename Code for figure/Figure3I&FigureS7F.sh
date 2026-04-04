########################################
###########Figure3I&FigureS7F###########
########################################
######################
###1.SNP-SELEX data###
######################
##http://renlab.sdsc.edu/GVATdb/index.html
##PBS negative value means that the TF prefers alternative allele
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
#############batch1
df1=fread("/data1/gy/public/SNP-SELEX/GVATdb.csv")
df1$chr <- sapply(str_split(df1$oligo, ":"), function(x) head(x, n = 1)) 
df1$pos <- sapply(str_split(df1$oligo, "-"), function(x) tail(x, n = 1)) 
df1$pos=as.numeric(df1$pos)-20
df1$snp=paste(df1$chr,df1$pos,df1$ref,df1$alt,sep=":")
df1=df1[,c("TF","snp","oligo_auc","oligo_pval","pbs","pval")]
df1$batch="batch1"
#############batch2
df2=fread("/data1/gy/public/SNP-SELEX/GVATdb.novel_batch.csv")
df2$snp <- gsub("_", ":", df2$snp) 
df2$batch="batch2"
##
df=rbind(df1,df2)
dim(df)
#[1] 2660658       7
length(unique(df$snp))
#[1] 151289
df <- df %>%  
  separate(snp, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = F)
df$chrbp=paste0(df$chr,":",df$pos) 
fwrite(df,"/data1/gy/public/SNP-SELEX/GVATdb.merged_2batch.final.txt",col.names=T,row.names=F,sep="\t",quote=F)

##########################################################
###2.caSNP(in corresponding peak) in SNP-SELEX database###
##########################################################
library(purrr)
##########################caSNP
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt")
caQTL=caQTL %>% filter(pheno_var_dist == 0)
caQTL=caQTL[,c(1,8,12,13)]
caQTL <- caQTL %>%  
  separate(var_id, into = c("chr_caSNP", "pos_caSNP", "ref_caSNP", "alt_caSNP"), sep = ":", remove = F)
caQTL$chrbp=paste0(caQTL$chr_caSNP,":",caQTL$pos_caSNP) 
##In merging the SNP-SELEX results, as PBS = refscore - altscore, a negative value signifies stronger binding for the ALT allele. Consequently, PBS values must be adjusted to reflect a positive binding strength.
caQTL=inner_join(caQTL,df,by=c("chrbp"="chrbp"))
caQTL <- caQTL %>%  
  mutate(pbs_new = case_when(  
    ref_caSNP == ref & alt_caSNP == alt ~ -pbs, 
    ref_caSNP == alt & alt_caSNP == ref ~ pbs, 
    TRUE ~ NA_real_    
  ))
##
any(is.na(caQTL)) 
#[1] FALSE
dim(caQTL)
#[1] 12330    20
length(unique(caQTL$var_id))
#[1] 861
write.csv(caQTL, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SNP-SELEX_cor/caSNP_in_corresponing_peak_SNP-SELEX_intersected.csv", row.names = FALSE)

###########################correlation analysis
#
thresholds <- 0:17  
#
results <- data.frame(threshold = integer(),  
                      pair_count = integer(),  
                      correlation_coefficient = numeric(),  
                      p_value = numeric(),  
                      ci_lower = numeric(),  
                      ci_upper = numeric(),  
                      stringsAsFactors = FALSE)  

#
for (threshold in thresholds) {  
  #
  filtered_data <- caQTL %>%  
    filter(abs(pbs_new) > threshold)  
  
  #
  pair_count <- nrow(filtered_data)  
  
  #
  if (pair_count > 1) { 
    cor_test_result <- cor.test(filtered_data$beta, filtered_data$pbs_new, use = "complete.obs")  
    
    #
    correlation_coefficient <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
  
    #
    n <- sum(complete.cases(filtered_data$beta,
                            filtered_data$pbs_new))
  
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
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SNP-SELEX_cor") 
write.csv(results, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SNP-SELEX_cor/caSNP_in_corresponing_peak_SNP-SELEX_correlation_results_by_threshold.csv", row.names = FALSE)

##########################################################
###3.asSNP(in corresponding peak) in SNP-SELEX database###
##########################################################
library(purrr)
##########################asSNP
ASCA=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1")
ASCA=ASCA[,c(3,6,10,11)]
ASCA$ALT_REF_ratio_log2 <- log2(ASCA$ALL.AF / (1 - ASCA$ALL.AF)) 
ASCA <- ASCA %>%  
  separate(RSID, into = c("chr_asSNP", "pos_asSNP", "ref_asSNP", "alt_asSNP"), sep = ":", remove = F)
ASCA$chrbp=paste0(ASCA$chr_asSNP,":",ASCA$pos_asSNP) 
##
ASCA=inner_join(ASCA,df,by=c("chrbp"="chrbp"))
ASCA <- ASCA %>%  
  mutate(pbs_new = case_when(  
    ref_asSNP == ref & alt_asSNP == alt ~ -pbs, 
    ref_asSNP == alt & alt_asSNP == ref ~ pbs,
    TRUE ~ NA_real_ 
  ))
##
any(is.na(ASCA)) 
#[1] FALSE
dim(ASCA)
#[1] 25908    22
length(unique(ASCA$RSID))
#[1] 1759
write.csv(ASCA, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SNP-SELEX_cor/asSNP_SNP-SELEX_intersected.csv", row.names = FALSE)

###########################correlation analysis
#
thresholds <- 0:17  
#
results <- data.frame(threshold = integer(),  
                      pair_count = integer(),  
                      correlation_coefficient = numeric(),  
                      p_value = numeric(),  
                      ci_lower = numeric(),  
                      ci_upper = numeric(),  
                      stringsAsFactors = FALSE)  

#
for (threshold in thresholds) {  
  #
  filtered_data <- ASCA %>%  
    filter(abs(pbs_new) > threshold)  
  
  # 
  pair_count <- nrow(filtered_data)  
  
  #
  if (pair_count > 1) { 
    cor_test_result <- cor.test(filtered_data$ALT_REF_ratio_log2, filtered_data$pbs_new, use = "complete.obs")  
    
    #
    correlation_coefficient <- cor_test_result$estimate
    p_value <- cor_test_result$p.value
  
    #
    n <- sum(complete.cases(filtered_data$ALT_REF_ratio_log2,
                            filtered_data$pbs_new))
  
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
write.csv(results, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SNP-SELEX_cor/asSNP_SNP-SELEX_correlation_results_by_threshold.csv", row.names = FALSE)

###################
###4.scatterplot###
###################
setwd("/data1/gy/ATAC_for_review/Figure3I&FigureS7F/output")
library(ggplot2)  
####################
###i.caSNP inpeak###
####################
caSNP_result=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SNP-SELEX_cor/caSNP_in_corresponing_peak_SNP-SELEX_correlation_results_by_threshold.csv")
# 
caSNP_result$threshold <- as.factor(results$threshold)  
##
p1 <- ggplot(caSNP_result, aes(x = threshold, y = correlation_coefficient)) +  
    geom_point(size = 3,alpha=0.75, color = "#094d92") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "#094d92",alpha=0.75) +
    geom_text(aes(label = pair_count), hjust = -0.3, vjust = 2, size = 2.5) + 
    scale_x_continuous(breaks = seq(0, max(caSNP_result$threshold) + 2, by = 2)) + 
    labs(title = "",  
         x = "Minimal absolute SNP-SELEX SNP-TF pair PBS value",  
         y = "Pearson correlation coefficient between\nSNP-SELEX SNP-TF pair PBS values and\ncaSNP (in corresponding peak) caQTL effect") +  
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
        axis.line = element_blank()
    )   
ggsave(p1,filename="caSNP_in_corresponing_peak_SNP-SELEX_correlation_results_by_threshold.Scatterplot.pdf") ##Figure3I

##############
###ii.asSNP###
##############
asSNP_result=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/SNP-SELEX_cor/asSNP_SNP-SELEX_correlation_results_by_threshold.csv")
#
asSNP_result$threshold <- as.factor(results$threshold)  
##
p2 <- ggplot(asSNP_result, aes(x = threshold, y = correlation_coefficient)) +  
    geom_point(size = 3,alpha=0.75, color = "#238444") + 
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "#238444",alpha=0.75) +  
    geom_text(aes(label = pair_count), hjust = -0.3, vjust = 2, size = 2.5) +
    scale_x_continuous(breaks = seq(0, max(caSNP_result$threshold) + 2, by = 2)) + 
    labs(title = "",  
         x = "Minimal absolute SNP-SELEX SNP-TF pair PBS value",  
         y = "Pearson correlation coefficient between\nSNP-SELEX SNP-TF pair PBS values and\nasSNP ASCA log2(ALT/REF)") +  
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
        axis.line = element_blank() 
    )   
ggsave(p2,filename="asSNP_SNP-SELEX_correlation_results_by_threshold.Scatterplot.pdf")  ##FigureS7F
