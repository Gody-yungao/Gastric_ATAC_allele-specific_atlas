################################
############Figure4B############
################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
all_data=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4.chr1_22.csv")
all_data_0.5=subset(all_data,all_data$PPH4>=0.5)

##########caOCR num
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt")
length(unique(caQTL$pheno_id))
#[1] 13665

##########eGene num
eQTL=fread("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR0.05.with_symbol.all.chr.txt")
length(unique(eQTL$phenotype_id))
#[1] 11232


##########colocalizd caOCR-eGene pair
dim(unique(all_data_0.5))
#[1] 2053   13
##########colocalizd caOCR num
length(unique(all_data_0.5$Peak))
#[1] 1350
##########colocalizd eGene num
length(unique(all_data_0.5$Geneid))
#[1] 1205

########Adobe Illustrator plot venn plot
