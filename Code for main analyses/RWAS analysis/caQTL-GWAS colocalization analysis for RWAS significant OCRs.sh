##############################################################################
#########caQTL-GWAS colocalization analysis for RWAS significant OCRs#########
##############################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(coloc)
library(dplyr)
library(data.table)
# input
caqtl=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/caQTL_GWAS_coloc/input/95sample_caQTL_200kb_nominal_result.in_EAS.for_coloc")
gwas=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/caQTL_GWAS_coloc/input/GWAS_QCandFilter.metaResult_new.in_EAS.for_coloc")
Peak=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/caQTL_GWAS_coloc/input/Peak_file.for_coloc")
sig=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter/sig/RWAS.fdr0.1.sig.result.txt")
Peak=subset(Peak,Peak$Geneid %in% sig$ID)

#
smy_list <- list()  
#
for (i in 1:nrow(sig)) {  
  #
  current_peakid <- Peak$Geneid[i]  
  current_chr <- Peak$X.Chr[i]  
  current_start <- Peak$Start[i]+1  
  current_end <- Peak$End[i]  
  #
  caqtl_summary <- caqtl[Peakid == current_peakid]  
  #
  gwas_summary <- gwas[chr == current_chr & bp >= current_start - 210000 & bp <= current_end + 210000]  
  #
  summary <- left_join(caqtl_summary, gwas_summary, by=c('caqtlID'='caqtlID'))
  #
  if (nrow(summary) > 0) {  
    result <- coloc.abf(  
      dataset1 = list(snp = summary$caqtlID, pvalues = summary$Pgwas, type = "cc", s = 0.075, N = 223476, MAF = summary$MAF_gwas),  
      dataset2 = list(snp = summary$caqtlID, pvalues = summary$Pcaqtl, type = "quant", N = 95, MAF = summary$MAF_caqtl)  
    )  
    #
    max_snp_row <- result$results[which.max(result$results$SNP.PP.H4), ]
    max_snp <- max_snp_row$snp
    max_snp.PPH4 <- max_snp_row$SNP.PP.H4
    #
    colocout <- data.frame(  
      Peak = current_peakid,  
      Start = current_start,  
      End = current_end,  
      PPH4 = result$summary[['PP.H4.abf']],  
      PPH3 = result$summary[['PP.H3.abf']],  
      PPH2 = result$summary[['PP.H2.abf']],  
      PPH1 = result$summary[['PP.H1.abf']],  
      HitSNP = length(unique(summary$SNP)),  
      minGWASP = min(summary$Pgwas),
      mincaQTLP = min(summary$Pcaqtl),
      coloc_topSNP_caqtlID = max_snp,
      coloc_topSNP_PPH4 = max_snp.PPH4
    )  
    
    smy_list[[i]] <- colocout
  }  
}  

#
smy_all <- do.call(rbind, smy_list)
smy_all_df=as.data.frame(smy_all)
smy_all_df_0.5=subset(smy_all_df,smy_all_df$PPH4>=0.5)
###
write.csv(smy_all_df,'/data1/gy/ATACseq_RWAS/RWAS_STITCH/caQTL_GWAS_coloc/coloc.abf_result/RWAS_sig_Peak.caQTL_GWAS_coloc_200kb_withRegion_PPH4.csv',quote=F,row.names=F)
write.csv(smy_all_df_0.5,'/data1/gy/ATACseq_RWAS/RWAS_STITCH/caQTL_GWAS_coloc/coloc.abf_result/RWAS_sig_Peak.caQTL_GWAS_coloc_200kb_withRegion_PPH4_0.5.csv',quote=F,row.names=F)
q()
