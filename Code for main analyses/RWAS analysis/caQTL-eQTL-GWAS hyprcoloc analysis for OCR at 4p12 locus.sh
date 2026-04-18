##############################################################################
###########caQTL-eQTL-GWAS hyprcoloc analysis for OCR at 4p12 locus###########
##############################################################################
#######################################################################################################
###########1.4p12 locus chr4:48072992-48073492 caQTL-NFXL1 bulk eQTL-GWAS hyprcoloc analysis###########
#######################################################################################################
mkdir -p /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/caQTL_262eQTL_GWAS_hyprcoloc
conda activate R_base
library(data.table)
library(hyprcoloc)
library(dplyr)
###############################################################
#############i.caQTL eQTL GWAS intersected SNPlist#############
###############################################################
snplist=data.table(fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/caQTL_eQTL_intersected_SNPlist.txt"))
snplist=snplist[,c("chrbp","eqtlID","caqtlID","MAF_eqtl","MAF_caqtl")]
GWAS=fread("/data1/gy/EAS_GWAS_meta/GWAS_QCandFilter.new.metaResult")
GWAS$chrbp=paste0("chr",GWAS$chr,":",GWAS$bp)
GWAS$MAF_gwas <- pmin(GWAS$Freq1, 1 - GWAS$Freq1)  
GWAS=GWAS[,c("RSID","chrbp","Allele1","Allele2","MAF_gwas")]
snplist=merge(snplist,GWAS,by="chrbp")

#
snplist$REF <- sub(".*:([A-Z])(:[A-Z])$", "\\1", snplist$caqtlID)  
snplist$ALT <- sub(".*:([A-Z]):([A-Z])$", "\\2", snplist$caqtlID)  
eqtl_geno <- fread(
  "/data1/gy/262_eQTL_final_result/geno/262geno.all.chr.Rsq03.MAF01.G05.M05.Hwe6e.geno.withheader.txt",
  select = 1:5
)
eqtl_geno=eqtl_geno[,3:5]
colnames(eqtl_geno)[2:3]=c("ref","alt")
snplist=left_join(snplist,eqtl_geno,by=c("eqtlID"="ID"))

#######
snplist1 = snplist[which((snplist$Allele1 == snplist$REF & snplist$Allele2 == snplist$ALT) | (snplist$Allele1 == snplist$ALT & snplist$Allele2 == snplist$REF)),]
dim(snplist1)
#[1] 4161830      13
#######
snplist2 = snplist[!which((snplist$Allele1 == snplist$REF & snplist$Allele2 == snplist$ALT) | (snplist$Allele1 == snplist$ALT & snplist$Allele2 == snplist$REF)),]
dim(snplist2)
#[1] 268  13
#######
snplist2$nchar1 = nchar(snplist2$Allele1) + nchar(snplist2$Allele2)
snplist2$nchar2 = nchar(snplist2$REF) + nchar(snplist2$ALT)
snplist2 = snplist2[which(snplist2$nchar1 == snplist2$nchar2),]
snplist2_1 = snplist2[which(snplist2$nchar1 == 2),]
snplist2_2 = snplist2[which(snplist2$nchar1 > 2),]
snplist2_1 = snplist2_1[which(snplist2_1$Allele1 != snplist2_1$REF & snplist2_1$Allele2 != snplist2_1$ALT & 
                        snplist2_1$Allele1 != snplist2_1$ALT & snplist2_1$Allele2 != snplist2_1$REF),]
snplist2_1 = snplist2_1[,c(1:13)] 
dim(snplist2_1)
#[1] 0 13
#######
snplist = rbind(snplist1,snplist2_1)

##########################################################
#############ii.caQTL seQTL GWAS beta se file#############
##########################################################
#######caqtl
peakID="chr4:48072992-48073492"
caqtl=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/caQTL_input_per_chr/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_eQTL_coloc.chr4.txt")
caqtl=subset(caqtl,pheno_id %in% peakID)
caqtl=caqtl[,c("pheno_id","var_id","beta","p_nominal")]
caqtl$se=sqrt(((caqtl$beta)^2)/qchisq(caqtl$p_nominal,1,lower.tail=F))
colnames(caqtl)=c("peakID","caqtlID","caqtlBeta","caqtlP","caqtlSE")
#######262 bulk eqtl
symbol="NFXL1"
eqtl=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/eQTL_input_per_chr/262sample_ciseQTL_1MB_result.cis_qtl_pairs.chr1_22.with_symbol.for_caQTL_eQTL_coloc.chr4.txt")
eqtl=subset(eqtl,gene_name %in% symbol)
eqtl=eqtl[,c("phenotype_id","gene_name","variant_id","slope","pval_nominal","slope_se")]
colnames(eqtl)[1]=c("gene_id")
colnames(eqtl)[3:6]=c("eqtlID","eqtlBeta","eqtlP","eqtlSE")
#######GWAS
GWAS=fread("/data1/gy/EAS_GWAS_meta/GWAS_QCandFilter.new.metaResult")
GWAS=GWAS[,c("RSID","Effect","P-value","StdErr")]
colnames(GWAS)[2:4]=c("GWASBeta","GWASP","GWASSE")
#######merge
snplist=merge(snplist,caqtl,by="caqtlID")
snplist=merge(snplist,eqtl,by="eqtlID")
snplist=merge(snplist,GWAS,by="RSID")
dim(snplist)
#[1] 840  25

######################################################
#############iii.correct allele direction#############
######################################################
complementary <- function(bases) {  
  sapply(bases, function(base) {  
    switch(base,  
           C = "G",  
           G = "C",  
           A = "T",  
           T = "A",  
           base)
  })  
}  

# eQTL beta
snplist$eqtlBeta_new <- ifelse(  
  (snplist$ref == snplist$REF & snplist$alt == snplist$ALT) |  
  (snplist$ref == complementary(snplist$REF) & snplist$alt == complementary(snplist$ALT)),  
  snplist$eqtlBeta,
  -snplist$eqtlBeta
)  
# GWAS beta
snplist$GWASBeta_new <- ifelse(  
  (snplist$Allele2 == snplist$REF & snplist$Allele1 == snplist$ALT) |  
  (snplist$Allele2 == complementary(snplist$REF) & snplist$Allele1 == complementary(snplist$ALT)),  
  snplist$GWASBeta,
  -snplist$GWASBeta
)   

###########################################
#############iv.hyprcoloc test#############
###########################################
# hypercoloc input
##beta
betas=as.data.frame(snplist[,c("chrbp","caqtlBeta","eqtlBeta_new","GWASBeta_new")])
rownames(betas)=betas$chrbp
betas=betas[,-1]
colnames(betas)=c("caQTL","eQTL","GWAS")
##se
ses=as.data.frame(snplist[,c("chrbp","caqtlSE","eqtlSE","GWASSE")])
rownames(ses)=ses$chrbp
ses=ses[,-1]
colnames(ses)=c("caQTL","eQTL","GWAS")
##
traits <- c("caQTL","eQTL","GWAS")
rsid <- rownames(betas)
##
result=hyprcoloc(effect.est=as.matrix(betas), effect.se=as.matrix(ses), trait.names = traits, snp.id = rsid,bb.alg =F )

#
hyprcoloc_out <- data.frame(  
      Peak = peakID,
      Geneid = snplist$gene_id[1],
      Symbol = symbol,
      posterior_prob = result$results[['posterior_prob']], 
      regional_prob =  result$results[['regional_prob']],
      candidate_snp = result$results[['candidate_snp']],
      posterior_explained_by_snp = result$results[['posterior_explained_by_snp']],
      HitSNP = length(unique(snplist$chrbp)),  
      mincaQTLP = min(snplist$caqtlP),
      mineQTLP = min(snplist$eqtlP),
      minGWASP = min(snplist$GWASP)
    )
hyprcoloc_out
#                    Peak             Geneid Symbol posterior_prob regional_prob
#1 chr4:48072992-48073492 ENSG00000170448.11  NFXL1         0.4692        0.8019
#  candidate_snp posterior_explained_by_snp HitSNP   mincaQTLP     mineQTLP
#1 chr4:48076486                     0.1205    840 2.02282e-10 1.769404e-05
#   minGWASP
#1 1.232e-05

##
outdir="/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/caQTL_262eQTL_GWAS_hyprcoloc/"
write.csv(hyprcoloc_out,file=paste0(outdir,'caQTL_eQTL_GWAS_hyprcoloc_200kb_1MB_PP.result.csv'),quote=F,row.names=F)
q()

##############################################################################################################
###########2.4p12 locus chr4:48072992-48073492 caQTL-NIPAL1 epithelium eQTL-GWAS hyprcoloc analysis###########
##############################################################################################################
mkdir -p /data1/gy/ATACseq_RWAS/GC_candidate_region_final_262_STITCH/chr4:48072992-48073492/caQTL_seQTL_GWAS_hyprcoloc/input
conda activate R_base
R
library(data.table)
library(hyprcoloc)
library(dplyr)
###############################################################
#############i.caQTL eQTL GWAS intersected SNPlist#############
###############################################################
snplist=data.table(fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/epi/input/caQTL_sceQTL_intersected_SNPlist.txt"))
snplist=snplist[,c("chrbp","eqtlID","caqtlID","MAF_eqtl","MAF_caqtl")]
GWAS=fread("/data1/gy/EAS_GWAS_meta/GWAS_QCandFilter.new.metaResult")
GWAS$chrbp=paste0("chr",GWAS$chr,":",GWAS$bp)
GWAS$MAF_gwas <- pmin(GWAS$Freq1, 1 - GWAS$Freq1)  
GWAS=GWAS[,c("RSID","chrbp","Allele1","Allele2","MAF_gwas")]
snplist=merge(snplist,GWAS,by="chrbp")

# 
snplist$REF <- sub(".*:([A-Z])(:[A-Z])$", "\\1", snplist$caqtlID)  
snplist$ALT <- sub(".*:([A-Z]):([A-Z])$", "\\2", snplist$caqtlID)  
snplist$ref <- sub(".*:([A-Z])(:[A-Z])$", "\\1", snplist$eqtlID)  
snplist$alt <- sub(".*:([A-Z]):([A-Z])$", "\\2", snplist$eqtlID) 

#######
snplist1 = snplist[which((snplist$Allele1 == snplist$REF & snplist$Allele2 == snplist$ALT) | (snplist$Allele1 == snplist$ALT & snplist$Allele2 == snplist$REF)),]
dim(snplist1)
#[1] 3844531      13
#######
snplist2 = snplist[!which((snplist$Allele1 == snplist$REF & snplist$Allele2 == snplist$ALT) | (snplist$Allele1 == snplist$ALT & snplist$Allele2 == snplist$REF)),]
dim(snplist2)
#[1] 513  13
#######
snplist2$nchar1 = nchar(snplist2$Allele1) + nchar(snplist2$Allele2)
snplist2$nchar2 = nchar(snplist2$REF) + nchar(snplist2$ALT)
snplist2 = snplist2[which(snplist2$nchar1 == snplist2$nchar2),]
snplist2_1 = snplist2[which(snplist2$nchar1 == 2),]
snplist2_2 = snplist2[which(snplist2$nchar1 > 2),]
snplist2_1 = snplist2_1[which(snplist2_1$Allele1 != snplist2_1$REF & snplist2_1$Allele2 != snplist2_1$ALT & 
                        snplist2_1$Allele1 != snplist2_1$ALT & snplist2_1$Allele2 != snplist2_1$REF),]
snplist2_1 = snplist2_1[,c(1:13)] 
dim(snplist2_1)
#[1] 0 13
#######
snplist = rbind(snplist1,snplist2_1)

##########################################################
#############ii.caQTL seQTL GWAS beta se file#############
##########################################################
#######caqtl
peakID="chr4:48072992-48073492"
caqtl=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/epi/input/caQTL_input_per_chr/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_sceQTL_coloc.chr4.txt")
caqtl=subset(caqtl,pheno_id %in% peakID)
caqtl=caqtl[,c("pheno_id","var_id","beta","p_nominal")]
caqtl$se=sqrt(((caqtl$beta)^2)/qchisq(caqtl$p_nominal,1,lower.tail=F))
colnames(caqtl)=c("peakID","caqtlID","caqtlBeta","caqtlP","caqtlSE")
#######sceqtl
symbol="NIPAL1"
sceqtl=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/epi/input/eQTL_input_per_chr/epi_sceQTL_result.allpairs.withsymbol.withFDR.for_caQTL_sceQTL_coloc.chr4.txt")
sceqtl=subset(sceqtl,gene_name %in% symbol)
sceqtl=sceqtl[,c("gene_id","gene_name","variant_id","slope","pval_nominal","slope_se")]
colnames(sceqtl)[3:6]=c("eqtlID","eqtlBeta","eqtlP","eqtlSE")
#######GWAS
GWAS=fread("/data1/gy/EAS_GWAS_meta/GWAS_QCandFilter.new.metaResult")
GWAS=GWAS[,c("RSID","Effect","P-value","StdErr")]
colnames(GWAS)[2:4]=c("GWASBeta","GWASP","GWASSE")
#######merge
snplist=merge(snplist,caqtl,by="caqtlID")
snplist=merge(snplist,sceqtl,by="eqtlID")
snplist=merge(snplist,GWAS,by="RSID")
dim(snplist)
#[1] 842  25

######################################################
#############iii.correct allele direction#############
######################################################
complementary <- function(bases) {  
  sapply(bases, function(base) {  
    switch(base,  
           C = "G",  
           G = "C",  
           A = "T",  
           T = "A",  
           base)
  })  
}  

# eQTL bera
snplist$eqtlBeta_new <- ifelse(  
  (snplist$ref == snplist$REF & snplist$alt == snplist$ALT) |  
  (snplist$ref == complementary(snplist$REF) & snplist$alt == complementary(snplist$ALT)),  
  snplist$eqtlBeta,
  -snplist$eqtlBeta
)  
# GWAS beta
snplist$GWASBeta_new <- ifelse(  
  (snplist$Allele2 == snplist$REF & snplist$Allele1 == snplist$ALT) |  
  (snplist$Allele2 == complementary(snplist$REF) & snplist$Allele1 == complementary(snplist$ALT)),  
  snplist$GWASBeta,
  -snplist$GWASBeta
)   

###########################################
#############iv.hyprcoloc test#############
###########################################
# hypercoloc input
##beta
betas=as.data.frame(snplist[,c("chrbp","caqtlBeta","eqtlBeta_new","GWASBeta_new")])
rownames(betas)=betas$chrbp
betas=betas[,-1]
colnames(betas)=c("caQTL","eQTL","GWAS")
##se
ses=as.data.frame(snplist[,c("chrbp","caqtlSE","eqtlSE","GWASSE")])
rownames(ses)=ses$chrbp
ses=ses[,-1]
colnames(ses)=c("caQTL","eQTL","GWAS")
##
traits <- c("caQTL","eQTL","GWAS")
rsid <- rownames(betas)
##
result=hyprcoloc(effect.est=as.matrix(betas), effect.se=as.matrix(ses), trait.names = traits, snp.id = rsid,bb.alg =F )

#
hyprcoloc_out <- data.frame(  
      Peak = peakID,
      Geneid = snplist$gene_id[1],
      Symbol = symbol,
      posterior_prob = result$results[['posterior_prob']], 
      regional_prob =  result$results[['regional_prob']],
      candidate_snp = result$results[['candidate_snp']],
      posterior_explained_by_snp = result$results[['posterior_explained_by_snp']],
      HitSNP = length(unique(snplist$chrbp)),  
      mincaQTLP = min(snplist$caqtlP),
      minsceQTLP = min(snplist$eqtlP),
      minGWASP = min(snplist$GWASP)
    )  

hyprcoloc_out
#                    Peak          Geneid Symbol posterior_prob regional_prob
#1 chr4:48072992-48073492 ENSG00000163293 NIPAL1         0.6922        0.9531
#  candidate_snp posterior_explained_by_snp HitSNP   mincaQTLP  minsceQTLP
#1 chr4:48090040                      0.471    842 2.02282e-10 9.58249e-18
#   minGWASP
#1 1.232e-05

##
outdir="/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/caQTL_seQTL_GWAS_hyprcoloc/"
write.csv(hyprcoloc_out,file=paste0(outdir,'caQTL_sceQTL_GWAS_hyprcoloc_200kb_1MB_PP.result.csv'),quote=F,row.names=F)
q()
