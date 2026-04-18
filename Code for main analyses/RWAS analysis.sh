###################################
###########RWAS analysis###########
###################################
##ref: https://github.com/scbaca/cwas/tree/master PMID:36071171
####################################################################
#################1.input file RWAS weight calculation###############
####################################################################
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/input/ASCA_matrix_perchr
cd /data1/gy/ATACseq_RWAS/RWAS_STITCH/input/ASCA_matrix_perchr
##########input ASCA count matrix
VCF=/data1/gy/ATACseq_RWAS/ASCA_STITCH/input/3.stratas_prep_ase_vcf/95sample_BYES_merged.sort.vcf.gz
zcat $VCF \
| grep -v '#' \
| cut -f 1-5,9- \
| awk 'BEGIN{OFS="\t"} {for(i=1;i<=5;i++) printf "%s ",$i; for(i=6;i<=NF;i++) {gsub(":"," ",$i); printf "%s\t",$i}; print ""}' \
| tr ',' '\t' \
| sed 's/[|,]/ /g' \
| tr '\t' ' ' \
| sed 's/GT AS//g' \
| while read -r line; do
    chr=$(echo "$line" | awk '{print $1}')
    if [[ "$chr" =~ ^chr([1-9]|1[0-9]|2[0-2])$ ]]; then
        echo "$line" >> "95sample_BYES_merged.MAT.${chr}.txt"
    fi
done

##########total OCR accessibility TMM_qnorm matrix
cd /data1/gy/ATACseq_RWAS/RWAS_STITCH/input
zcat /data1/gy/ATACseq_RWAS/caQTL_STITCH/exp/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.bed.gz | cut -f4,7- > 95sample_IterativeOverlapPeakSet.final.TMM_qnorm.total.mat
##covariate
cp /data1/gy/ATACseq_RWAS/caQTL_STITCH/multiple_model_covar_all/TMM/15.95sample_sex_age_Hp_PC1_5_peer1_7.covar ./
##position file
awk -F'\t' '{print $2}' \
/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/plink/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.bim \
> 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.pos


#########################################################
#################2.RWAS weight calculation###############
#########################################################
conda activate RWAS
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/predict/WEIGHTS
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/predict/results
cd /data1/gy/ATACseq_RWAS/RWAS_STITCH/predict
###########±200kb of OCR (±200.25kb of OCR center)
for i in {1..22}; do
Rscript /Public/gaoyun/ASE_ATACseq/stratAS_code/stratas.R \
--input /data1/gy/ATACseq_RWAS/RWAS_STITCH/input/ASCA_matrix_perchr/95sample_BYES_merged.MAT.chr${i}.txt \
--seed 123 \
--max_rho 0.2 \
--window 200250 \
--min_cov 1 \
--covar /data1/gy/ATACseq_RWAS/caQTL_STITCH/multiple_model_covar_all/TMM/15.95sample_sex_age_Hp_PC1_5_peer1_7.covar \
--total_matrix /data1/gy/ATACseq_RWAS/RWAS_STITCH/input/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.total.mat \
--fill_cnv FALSE \
--samples /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/baseline/95sample_BYES.AS.PHE \
--peaks /data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.for_stratas_input.sort.bed \
--global_param /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/5.stratas_input_params_file/95sample_BYES_withheader.global.params \
--predict_snps /data1/gy/ATACseq_RWAS/RWAS_STITCH/input/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.pos \
--predict --predict_only \
> ./results/predict.95sample_BYES_output_200kb_chr${i}.profile
done
##########combine
cat /data1/gy/software/cwas-master/misc/predict.head ./results/predict.95sample_BYES_output_200kb_chr*.profile > ./results/predict.95sample_BYES_output_200kb_all_chr.profile 

#################################################
#########3.prepare input files for FUSION#########
#################################################
mkdir -p ../fusion/input
ls ./WEIGHTS/ | grep RDat > ../fusion/input/weights.tmp
sed 's/chr//' ../fusion/input/weights.tmp | \
sed 's/.wgt.RDat//' | sed 's/:/\t/' | \
sed 's/-/\t/' > ../fusion/input/coords.tmp
paste ../fusion/input/weights.tmp \
<(sed 's/.wgt.RDat//' ../fusion/input/weights.tmp) \
../fusion/input/coords.tmp > ../fusion/input/comb.tmp
printf "WGT\tID\tCHR\tP0\tP1\n" > ../fusion/input/weights.pos
cat ../fusion/input/comb.tmp >> ../fusion/input/weights.pos
rm ../fusion/input/coords.tmp ../fusion/input/weights.tmp ../fusion/input/comb.tmp

################################################
#########4.prepare GWAS file for FUSION#########
################################################
################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
common=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/input/shared_snp_list.bed")
gwas=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/input/gwas/GWAS_QCandFilter.new.with_SNPid_RWAS.metaResult")
gwas=subset(gwas,gwas$SNPid_RWAS %in% common$V4)
gwas$REF <- sapply(strsplit(gwas$SNPid_RWAS,':'),'[',3)
gwas$ALT <- sapply(strsplit(gwas$SNPid_RWAS,':'),'[',4)

#######
gwas1 = gwas[which((gwas$REF == gwas$Allele1 & gwas$ALT == gwas$Allele2) | (gwas$REF == gwas$Allele2 & gwas$ALT == gwas$Allele1)),]
dim(gwas1)
#[1] 4697331      22
#######
gwas2 = gwas[!gwas$SNPid_RWAS %in% gwas1$SNPid_RWAS, ]
dim(gwas2)
#[1]  0 22

########
gwas$Effect <- ifelse(gwas$REF == gwas$Allele1 & gwas$ALT == gwas$Allele2,   
                      gwas$Effect * -1,   
                      gwas$Effect)  

########
gwas$Freq <- ifelse(gwas$REF == gwas$Allele1 & nwas$ALT == gwas$Allele2,   
                      1-gwas$Freq1,   
                      gwas$Freq1)  
########
gwas <- gwas[,c('SNPid_RWAS','ALT','REF','Freq','Effect','StdErr','P-value','N')]
fwrite(gwas,'/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/input/gwas/GWAS_QCandFilter.new.with_SNPid_RWAS.for_RWAS.metaResult',row.names=F,col.names=T,quote=F,sep='\t') 
q()

##################################convert to Z-SCORE
conda activate ldsc
##############
##Note that reversing the A1 and A2 alleles does not affect the RWAS results, as long as the SNP effect allele direction is consistent between the GWAS and RWAS model data.
python /data1/gy/software/ldsc/munge_sumstats.py \
--sumstats /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/input/gwas/GWAS_QCandFilter.new.with_SNPid_RWAS.for_RWAS.metaResult \
--snp SNPid_RWAS \
--a1 REF \
--a2 ALT \
--frq Freq \
--p P-value \
--N 223476 \
--N-cas 16817 \
--N-con 206659 \
--N-col N \
--out /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/input/gwas/GWAS_QCandFilter.new.with_SNPid_RWAS.without_ambiguous.for_RWAS.QC 

#######################################################################
#########5.run FUSION (accessibility model - GWAS association)#########
#######################################################################
conda activate fusion
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/log
#####
for mod in lasso lasso.as lasso.combined top1.as top1 top1.combined
do  
    for chr in $(seq 1 22)
    do  
        Rscript /data1/gy/software/cwas-master/fusion/FUSION.assoc_test.R \
        --sumstats /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/input/gwas/GWAS_QCandFilter.new.with_SNPid_RWAS.without_ambiguous.for_RWAS.QC.sumstats.gz \
        --out /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.${mod}.${chr}.dat \
        --weights /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/input/weights.pos \
        --weights_dir /data1/gy/ATACseq_RWAS/RWAS_STITCH/predict/WEIGHTS/ \
        --ref_ld_chr /data1/gy/ATACseq_RWAS/RWAS_STITCH/ldref/per_chr/1000G.EAS. \
        --chr ${chr} \
        --force_model ${mod} |tee /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/log/${mod}.${chr}.log 
    done  
done

#############################################
#########6.merge results from FUSION#########
#############################################
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged
HEADER=/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged/header.tmp
head -n 1 $(echo $(for s in {1..22}; do echo "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.lasso.$s.dat"; done) | sed 's/ .*//') > "$HEADER"  

#
output_lasso="/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged/fusion.lasso.txt"  
output_lasso_as="/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged/fusion.lasso.as.txt"  
output_lasso_combined="/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged/fusion.lasso.combined.txt"  
output_top1_as="/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged/fusion.top1.as.txt"  
output_top1_combined="/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged/fusion.top1.combined.txt"  
output_top1="/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged/fusion.top1.txt"  
cp "$HEADER" "$output_lasso"  
cp "$HEADER" "$output_lasso_as"  
cp "$HEADER" "$output_lasso_combined"  
cp "$HEADER" "$output_top1_as"  
cp "$HEADER" "$output_top1_combined"  
cp "$HEADER" "$output_top1"  

#chr1-22(except MHC region)
for s in {1..22}; do  
    cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.lasso.$s.dat" | grep -v PANEL >> "$output_lasso"  
    cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.lasso.as.$s.dat" | grep -v PANEL >> "$output_lasso_as"  
    cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.lasso.combined.$s.dat" | grep -v PANEL >> "$output_lasso_combined"  
    cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.top1.as.$s.dat" | grep -v PANEL >> "$output_top1_as"  
    cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.top1.combined.$s.dat" | grep -v PANEL >> "$output_top1_combined"  
    cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.top1.$s.dat" | grep -v PANEL >> "$output_top1"  
done  

#MHC region
cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.lasso.6.dat.MHC" | grep -v PANEL >> "$output_lasso"  
cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.lasso.as.6.dat.MHC" | grep -v PANEL >> "$output_lasso_as"  
cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.lasso.combined.6.dat.MHC" | grep -v PANEL >> "$output_lasso_combined"  
cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.top1.as.6.dat.MHC" | grep -v PANEL >> "$output_top1_as"  
cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.top1.combined.6.dat.MHC" | grep -v PANEL >> "$output_top1_combined"  
cat "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/meta.top1.6.dat.MHC" | grep -v PANEL >> "$output_top1"  

#
rm "$HEADER"  

#####################################################
#########7.filter merged results from FUSION#########
#####################################################
## included all caOCRs for analyses using the total-accessibility models (lasso.total and top1.total),
## all asOCRs for analyses using the allele-specific models (lasso.allelic and top1.allelic),
## and only OCRs defined as both caOCRs and asOCRs for the combined models (lasso.combined and top1.combined).
######
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged")
library(data.table)
library(qvalue)
library(ggplot2)
library(stringr)
##
caQTL_sig=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt")
caQTL_sig_peak=unique(caQTL_sig$pheno_id)
ASCA_sig=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1")
ASCA_sig_peak=unique(ASCA_sig$NAME)
intersect=intersect(caQTL_sig_peak,ASCA_sig_peak)
##
length(caQTL_sig_peak)
#[1] 13665
length(ASCA_sig_peak)
#[1] 5658
length(intersect)
#[1] 2924

############filter
##lasso.total
data_lasso=fread("fusion.lasso.txt")
data_lasso=subset(data_lasso,ID %in% caQTL_sig_peak)
fwrite(data_lasso,"../merged_filter/fusion.lasso.filter.txt",col.names=T,row.names=F,quote=F,sep="\t")
dim(data_lasso)
#[1] 9105   20
##lasso.allelic
data_lasso.as=fread("fusion.lasso.as.txt")
data_lasso.as=subset(data_lasso.as,ID %in% ASCA_sig_peak)
fwrite(data_lasso.as,"../merged_filter/fusion.lasso.as.filter.txt",col.names=T,row.names=F,quote=F,sep="\t")
dim(data_lasso.as)
#[1] 5658   20
##lasso.combined
data_lasso.combined=fread("fusion.lasso.combined.txt")
data_lasso.combined=subset(data_lasso.combined,ID %in% intersect)
fwrite(data_lasso.combined,"../merged_filter/fusion.lasso.combined.filter.txt",col.names=T,row.names=F,quote=F,sep="\t")
dim(data_lasso.combined)
#[1] 2924   20
##top1.total
data_top1=fread("fusion.top1.txt")
data_top1=subset(data_top1,ID %in% caQTL_sig_peak)
fwrite(data_top1,"../merged_filter/fusion.top1.filter.txt",col.names=T,row.names=F,quote=F,sep="\t")
dim(data_top1)
#[1] 9105   20
##top1.allelic
data_top1.as=fread("fusion.top1.as.txt")
data_top1.as=subset(data_top1.as,ID %in% ASCA_sig_peak)
fwrite(data_top1.as,"../merged_filter/fusion.top1.as.filter.txt",col.names=T,row.names=F,quote=F,sep="\t")
dim(data_top1.as)
#[1] 5658   20
##top1.combined
data_top1.combined=fread("fusion.top1.combined.txt")
data_top1.combined=subset(data_top1.combined,ID %in% intersect)
fwrite(data_top1.combined,"../merged_filter/fusion.top1.combined.filter.txt",col.names=T,row.names=F,quote=F,sep="\t")
dim(data_top1.combined)
#[1] 2924   20
q()

########################################################################################
#########8.combine RWAS results from different models and significant filtering#########
########################################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(dplyr)
library(qvalue)
library(ggplot2)
library(stringr)
##
fusion.path="/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter"
##
dir.create(paste0(fusion.path,"/sig"))                                         
dir.create(paste0(fusion.path,"/cv"))                                         

##combine results of all tested models
files=list.files(fusion.path, pattern="fusion.*.txt", full.name=T)
m = do.call(rbind, lapply(files, function(x) read.table(x, stringsAsFactors=F,header=T, sep="\t",
   colClasses=c(rep("character",3), rep("numeric",4), "character", "numeric", "character", rep("numeric", 5), rep("character",3), rep("numeric",2)))))

##replace each MODELCV.PV with the minimum of the two values it contains
p.tmp = matrix(as.numeric(unlist(str_split(m$MODELCV.PV, ","))), nrow=nrow(m), ncol=2, byrow=T)
p.best=suppressWarnings(apply(p.tmp, 1, FUN=function(x) min(x,na.rm=T))) #introduces Inf where both rows are NA
p.best[!is.finite(p.best)] = NA
m$MODELCV.PV=p.best
r2.tmp = matrix(as.numeric(unlist(str_split(m$MODELCV.R2, ","))), nrow=nrow(m), ncol=2, byrow=T)
m$MODELCV.R2=sapply(seq(1:nrow(r2.tmp)),FUN=function(x) ifelse(!is.na(p.tmp[x,1]) & p.best[x]==p.tmp[x,1], r2.tmp[x,1], r2.tmp[x,2]))

##remove models with no cross-validation or no GWAS SNP
m = subset(m, !is.na(MODELCV.PV) & !is.na(BEST.GWAS.Z)) 

##retain only models with CV p-value < 0.05, and for each OCR, keep only the model with the lowest CV p-value
m.cvsig = subset(m,MODELCV.PV < 0.05)
m.cvsig = m.cvsig %>% group_by(ID) %>% slice(which.min(MODELCV.PV))

##FDR correction and significant RWAS OCR (FDR<0.1)
m.cvsig$TWAS.fdr=p.adjust(m.cvsig$TWAS.P,method="BH")
m.cvsig.fdr.sig = subset(m.cvsig, TWAS.fdr < 0.1)
cat(paste0("significant RWAS from these models after fdr correction: ", nrow(m.cvsig.fdr.sig), "\n\n"))
#significant RWAS from these models after fdr correction: 58

##write
m.cvsig.fdr.sig$PANEL="RWAS"
write.table(m.cvsig, paste0(fusion.path,"/cv/RWAS.withFDR.result.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) 
write.table(m.cvsig.fdr.sig, paste0(fusion.path,"/sig/RWAS.fdr0.1.sig.result.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

##add caOCR/asOCR anno
library(data.table)
m.cvsig.fdr.sig=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter/sig/RWAS.fdr0.1.sig.result.txt")
caOCR=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.bed")
asOCR=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1.bed")
ca_ids <- caOCR$V4
as_ids <- asOCR$V4
m.cvsig.fdr.sig[, caOCR := ifelse(ID %chin% ca_ids, "Yes", "No")]
m.cvsig.fdr.sig[, asOCR := ifelse(ID %chin% as_ids, "Yes", "No")]
write.table(m.cvsig.fdr.sig,
            "/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter/sig/RWAS.fdr0.1.sig.result.with_caOCR_asOCR_anno.txt",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

###################################################################################
#########9.caQTL-GWAS colocalization analysis for 58 RWAS significant OCRs#########
###################################################################################
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
