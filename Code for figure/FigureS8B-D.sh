#########################
#######FigureS8B-D#######
#########################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
####result
coloc_anno_result=read.csv("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.all_pair.distal.TAD_loop_ABC_anno.csv")
coloc_anno_result$pair=paste0(coloc_anno_result$Peak,"~",coloc_anno_result$Geneid)

####input
caQTL_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_eQTL_coloc.txt")
eQTL_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/262sample_ciseQTL_1MB_result.cis_qtl_pairs.chr1_22.with_symbol.for_caQTL_eQTL_coloc.txt")

####################caQTL
# 
caQTL_coloc=subset(caQTL_coloc,pheno_id %in% coloc_anno_result$Peak)
#
caQTL_coloc[, c("REF", "ALT") := tstrsplit(var_id, ":", fixed = TRUE)[3:4]]  
#
caQTL_coloc=caQTL_coloc[,c(1,2,9,10,12,13,16,17,18)]
caQTL_coloc$chrbp=paste0(caQTL_coloc$var_chr,":",caQTL_coloc$var_start)
caQTL_coloc=caQTL_coloc[,-c(3,4)]
colnames(caQTL_coloc)=c("caQTL_SNPid","pheno_id","caQTL_Pval","caQTL_beta","caQTL_MAF","caQTL_REF","caQTL_ALT","chrbp")

####################eQTL
eQTL_coloc$phenotype_id=substr(eQTL_coloc$phenotype_id,1,15)
eQTL_coloc=subset(eQTL_coloc,phenotype_id %in% coloc_anno_result$Geneid) 
#
eQTL_coloc$chrbp <- stringr::str_extract(eQTL_coloc$variant_id, "^[^:]+:[^:]+")
eQTL_coloc$chrbp = paste0("chr", eQTL_coloc$chrbp)
eQTL_coloc=eQTL_coloc[,c(1,2,3,4,6,9,10,12,14)]   
colnames(eQTL_coloc)=c("phenotype_id","eQTL_SNPid","eQTL_REF","eQTL_ALT","eQTL_MAF","eQTL_Pval","eQTL_beta","gene_name","chrbp")

#####################coloc
coloc_anno_result=left_join(coloc_anno_result,caQTL_coloc,by=c("Peak"="pheno_id","coloc_topSNP_caqtlID"="caQTL_SNPid"))
coloc_anno_result=left_join(coloc_anno_result,eQTL_coloc,by=c("Geneid"="phenotype_id","chrbp"="chrbp"))
dim(coloc_anno_result)
#[1] 1752   30
#
coloc_anno_result <- coloc_anno_result[complete.cases(coloc_anno_result), ]  
dim(coloc_anno_result)
#[1] 1752   30

######################
complement <- function(base) {  
  sapply(base, function(b) {  
    switch(b,  
           "A" = "T",  
           "T" = "A",  
           "C" = "G",  
           "G" = "C",  
           b) 
  })  
}
#
coloc_anno_result=as.data.table(coloc_anno_result)
coloc_anno_result[, eQTL_beta_new := ifelse(  
  caQTL_REF == eQTL_REF & caQTL_ALT == eQTL_ALT,   
  eQTL_beta,   
  ifelse(caQTL_ALT == eQTL_REF & caQTL_REF == eQTL_ALT,   
         -eQTL_beta,   
         ifelse(complement(caQTL_REF) == eQTL_REF & complement(caQTL_ALT) == eQTL_ALT,   
                eQTL_beta,   
                ifelse(complement(caQTL_ALT) == eQTL_REF & complement(caQTL_REF) == eQTL_ALT,   
                       -eQTL_beta,   
                       NA))))]
#
library(dplyr)
coloc_anno_result <- coloc_anno_result %>% 
  mutate(
    sign_concordance = sign(caQTL_beta) == sign(eQTL_beta_new)
  )
  
#########################scatterplot
library(ggplot2)  
library(ggpubr)  
pph_threshold=0.5
############1.all distal coloc
prop_true <- mean(coloc_anno_result$sign_concordance, na.rm = TRUE)
prop_true  
#[1] 0.6957763
p1 <- ggplot(coloc_anno_result, aes(x = caQTL_beta, y = eQTL_beta_new)) +  
    geom_point(color = "#FE9929") +    
    geom_smooth(method = 'lm', se = FALSE) +  # Linear model without confidence interval  
    theme_bw() +  
    theme(  
      panel.grid = element_blank(),  
      panel.border = element_blank(),  
      axis.line = element_blank()  
    ) +  
    labs(x = "caQTL effect size of coloc top SNP", y = "eQTL effect size of coloc top SNP") +  
    stat_cor(method = "pearson", label.sep = "\n", r.digits = 3, p.digits = 3,   
             position = "identity", label.x = 1.0, label.y = 1.45) +  
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +  
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
    coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6))+
    ggtitle("All colocalized distal caOCR-eGene pairs")  # Title for the plot  

#cor
library(Rmpfr)
precise_cor_p <- function(x, y, prec = 256) {
  res <- cor.test(x, y, method = "pearson")
  t_mpfr <- mpfr(abs(res$statistic), precBits = prec)
  #
  p_mpfr <- 2 * pnorm(-t_mpfr, lower.tail = TRUE)
  list(r = res$estimate, t = res$statistic, df = res$parameter,
       p_precise = p_mpfr)
}
res=precise_cor_p(coloc_anno_result$caQTL_beta, coloc_anno_result$eQTL_beta_new)
format(res$p_precise, format = "e", digits = 3)
#[1] "1.33e-46"

##################2.loop anno
prop_true <- mean(coloc_anno_result[loop_anno=="TRUE"]$sign_concordance, na.rm = TRUE)
prop_true  
#[1] 0.8382979
p2 <- ggplot(coloc_anno_result[loop_anno=="TRUE"], aes(x = caQTL_beta, y = eQTL_beta_new)) +  
    geom_point(color = "#FDD0A2") +    
    geom_smooth(method = 'lm', se = FALSE) +  # Linear model without confidence interval  
    theme_bw() +  
    theme(  
      panel.grid = element_blank(),  
      panel.border = element_blank(),  
      axis.line = element_blank()  
    ) +  
    labs(x = "caQTL effect size of coloc top SNP", y = "eQTL effect size of coloc top SNP") +  
    stat_cor(method = "pearson", label.sep = "\n", r.digits = 3, p.digits = 3,   
             position = "identity", label.x = 1.0, label.y = 1.45) +  
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +  
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
    coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6)) +
    ggtitle("Colocalized pairs within chromatin loops")  # Title for the plot  

#cor
library(Rmpfr)
res=precise_cor_p(coloc_anno_result[loop_anno=="TRUE"]$caQTL_beta, coloc_anno_result[loop_anno=="TRUE"]$eQTL_beta_new)
format(res$p_precise, format = "e", digits = 3)
#[1] "1.99e-21"

###########3.TAD anno
prop_true <- mean(coloc_anno_result[TAD_anno=="TRUE"]$sign_concordance, na.rm = TRUE)
prop_true  
#[1] 0.7578521
p3 <- ggplot(coloc_anno_result[TAD_anno=="TRUE"], aes(x = caQTL_beta, y = eQTL_beta_new)) +  
    geom_point(color = "#FDAE6B") +    
    geom_smooth(method = 'lm', se = FALSE) +  # Linear model without confidence interval  
    theme_bw() +  
    theme(  
      panel.grid = element_blank(),  
      panel.border = element_blank(),  
      axis.line = element_blank()  
    ) +  
    labs(x = "caQTL effect size of coloc top SNP", y = "eQTL effect size of coloc top SNP") +  
    stat_cor(method = "pearson", label.sep = "\n", r.digits = 3, p.digits = 3,   
             position = "identity", label.x = 1.0, label.y = 1.45) +  
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +  
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
    coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6)) +
    ggtitle("Colocalized pairs within chromatin TAD")  # Title for the plot  

#cor
library(Rmpfr)
res=precise_cor_p(coloc_anno_result[TAD_anno=="TRUE"]$caQTL_beta, coloc_anno_result[TAD_anno=="TRUE"]$eQTL_beta_new)
format(res$p_precise, format = "e", digits = 3)
#[1] "1.87e-50"

###########4.ABC model anno
prop_true <- mean(coloc_anno_result[ABC_anno=="TRUE"]$sign_concordance, na.rm = TRUE)
prop_true  
#[1] 0.8017241
p4 <- ggplot(coloc_anno_result[ABC_anno=="TRUE"], aes(x = caQTL_beta, y = eQTL_beta_new)) +  
    geom_point(color = "#FD8D3C") +    
    geom_smooth(method = 'lm', se = FALSE) +  # Linear model without confidence interval  
    theme_bw() +  
    theme(  
      panel.grid = element_blank(),  
      panel.border = element_blank(),  
      axis.line = element_blank()  
    ) +  
    labs(x = "caQTL effect size of coloc top SNP", y = "eQTL effect size of coloc top SNP") +  
    stat_cor(method = "pearson", label.sep = "\n", r.digits = 3, p.digits = 3,   
             position = "identity", label.x = 1.0, label.y = 1.45) +  
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +  
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  
    coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6)) +
    ggtitle("Colocalized pairs with annotation of ABC model")  # Title for the plot  

#cor
library(Rmpfr)
res=precise_cor_p(coloc_anno_result[ABC_anno=="TRUE"]$caQTL_beta, coloc_anno_result[ABC_anno=="TRUE"]$eQTL_beta_new)
format(res$p_precise, format = "e", digits = 3)
#[1] "4.33e-36"

###########combine
setwd("/data1/gy/ATAC_for_review/FigureS8B-D/output")
library(cowplot)
combined <- plot_grid(p1, p2, p3, p4,
                      ncol = 2, nrow = 2)
ggsave("caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.all_pair.distal.TAD_loop_ABC_anno.coloc_top_SNP.caQTL_eQTL_beta.corplot.pdf",
       combined,
       width = 12, height = 12) ##FigureS8B-D
