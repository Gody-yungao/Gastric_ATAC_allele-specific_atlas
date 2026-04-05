###############################
###########FigureS8A###########
###############################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
####input
caQTL_eQTL_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv")
caQTL_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_eQTL_coloc.txt")
eQTL_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/262sample_ciseQTL_1MB_result.cis_qtl_pairs.chr1_22.with_symbol.for_caQTL_eQTL_coloc.txt")
####################caQTL
#
caQTL_coloc=subset(caQTL_coloc,pheno_id %in% caQTL_eQTL_coloc$Peak)
#
caQTL_coloc[, c("REF", "ALT") := tstrsplit(var_id, ":", fixed = TRUE)[3:4]]  
#
caQTL_coloc=caQTL_coloc[,c(1,2,9,10,12,13,16,17,18)]
caQTL_coloc$chrbp=paste0(caQTL_coloc$var_chr,":",caQTL_coloc$var_start)
caQTL_coloc=caQTL_coloc[,-c(3,4)]
colnames(caQTL_coloc)=c("caQTL_SNPid","pheno_id","caQTL_Pval","caQTL_beta","caQTL_MAF","caQTL_REF","caQTL_ALT","chrbp")

####################eQTL
eQTL_coloc=subset(eQTL_coloc,phenotype_id %in% caQTL_eQTL_coloc$Geneid) 
#
eQTL_coloc$chrbp <- stringr::str_extract(eQTL_coloc$variant_id, "^[^:]+:[^:]+")
eQTL_coloc$chrbp = paste0("chr", eQTL_coloc$chrbp)
eQTL_coloc=eQTL_coloc[,c(1,2,3,4,6,9,10,12,14)]   
colnames(eQTL_coloc)=c("phenotype_id","eQTL_SNPid","eQTL_REF","eQTL_ALT","eQTL_MAF","eQTL_Pval","eQTL_beta","gene_name","chrbp")

#####################coloc
caQTL_eQTL_coloc=left_join(caQTL_eQTL_coloc,caQTL_coloc,by=c("Peak"="pheno_id","coloc_topSNP_caqtlID"="caQTL_SNPid"))
caQTL_eQTL_coloc=left_join(caQTL_eQTL_coloc,eQTL_coloc,by=c("Geneid"="phenotype_id","chrbp"="chrbp"))
dim(caQTL_eQTL_coloc)
#[1] 2053   26
caQTL_eQTL_coloc$pair=paste0(caQTL_eQTL_coloc$Peak,"~",caQTL_eQTL_coloc$Geneid)
dim(caQTL_eQTL_coloc)
#[1] 2053   27
#
caQTL_eQTL_coloc <- caQTL_eQTL_coloc[complete.cases(caQTL_eQTL_coloc), ]  
dim(caQTL_eQTL_coloc)
#[1] 2053   27

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
caQTL_eQTL_coloc[, eQTL_beta_new := ifelse(  
  caQTL_REF == eQTL_REF & caQTL_ALT == eQTL_ALT,   
  eQTL_beta,   
  ifelse(caQTL_ALT == eQTL_REF & caQTL_REF == eQTL_ALT,   
         -eQTL_beta,   
         ifelse(complement(caQTL_REF) == eQTL_REF & complement(caQTL_ALT) == eQTL_ALT,   
                eQTL_beta,   
                ifelse(complement(caQTL_ALT) == eQTL_REF & complement(caQTL_REF) == eQTL_ALT,   
                       -eQTL_beta,   
                       NA))))]

########correlation analysis
library(Rmpfr)
precise_cor_p <- function(x, y, prec = 256) {
  res <- cor.test(x, y, method = "pearson")
  t_mpfr <- mpfr(abs(res$statistic), precBits = prec)
  #
  p_mpfr <- 2 * pnorm(-t_mpfr, lower.tail = TRUE)
  list(r = res$estimate, t = res$statistic, df = res$parameter,
       p_precise = p_mpfr)
}
res=precise_cor_p(caQTL_eQTL_coloc$caQTL_beta, caQTL_eQTL_coloc$eQTL_beta_new)
format(res$p_precise, format = "e", digits = 3)
#[1] "2.99e-78"

#########################scatterplot
setwd("/data1/gy/ATAC_for_review/FigureS8A/output")
library(ggplot2)  
library(ggpubr)  
pph_threshold=0.5
p <- ggplot(caQTL_eQTL_coloc, aes(x = caQTL_beta, y = eQTL_beta_new)) +  
    geom_point(color = "#D53E4F") +
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
    coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6))
    #ggtitle(paste("PPH4 Threshold:", pph_threshold))  # Title for the plot  
  
#
ggsave(filename = paste0("caQTL_eQTL_coloc_PPH4_", pph_threshold, ".coloc_top_SNP.caQTL_eQTL_beta.corplot.pdf"),
       plot = p, width = 5, height = 5)  ## FigureS8A
