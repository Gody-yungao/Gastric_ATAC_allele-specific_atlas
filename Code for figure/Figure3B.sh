#########################
#########Figure3B########
#########################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(ggplot2)
library(ggpubr)  
##
ASCA_sig=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1")
ASCA_sig=ASCA_sig[,c(1,3,4,5,6,12,13)]
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
caQTL_inpeak=subset(caQTL,caQTL$pheno_var_dist==0)
caQTL_inpeak=caQTL_inpeak[,c(1,7,8,9,12,13)]

#####
m1=merge(caQTL_inpeak,ASCA_sig,by.x="var_id",by.y="RSID")

#
library(Rmpfr)
#res <- cor.test(m1$C0.AF, m1$beta, method = "pearson")
precise_cor_p <- function(x, y, prec = 256) {
  res <- cor.test(x, y, method = "pearson")
  t_mpfr <- mpfr(abs(res$statistic), precBits = prec)
  #
  p_mpfr <- 2 * pnorm(-t_mpfr, lower.tail = TRUE)
  list(r = res$estimate, t = res$statistic, df = res$parameter,
       p_precise = p_mpfr)
}
res=precise_cor_p(m1$C0.AF, m1$beta)
format(res$p_precise, format = "e", digits = 3)
#[1] "7.35e-4187"

#
setwd("/data1/gy/ATAC_for_review/Figure3B/output")
p <- ggplot(m1, aes(x=C0.AF,y=beta)) +  
  geom_point(color = "#238444",alpha=0.75) + 
  geom_smooth(method = 'lm', se = TRUE) + 
  theme_bw() +  
  theme(  
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black") 
  ) +  
  labs(x = "asSNP allele fraction", y = "caQTL effect size (beta)") +  
  stat_cor(data = m1, method = "pearson", label.sep = "\n",   
           r.digits = 3, p.digits = 3, position = "identity",   
           label.x = 0.2, label.y = 1.4)
ggsave(p,file="asSNP_af.vs.caQTL_beta.corplot.pdf", height=5, width=5) ##Figure3B
