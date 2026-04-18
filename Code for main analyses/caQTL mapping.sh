###########################################
###############caQTL mapping###############
###########################################
################################################
##1.OCR accessibility bed.gz file for QTLtools##
################################################
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_STITCH/exp
#####
/Public/gaoyun/software/R-4.2.0/bin/R
TMM <- read.delim("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.TMM_qnorm",header = T)
TMM_bed=TMM[,c(2:4,1,5:101)]
TMM_bed$Start=TMM_bed$Start-1
colnames(TMM_bed)[1]="#Chr"
write.table(TMM_bed,"/data1/gy/ATACseq_RWAS/caQTL_STITCH/exp/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.bed",sep="\t",col.names=T,row.names=F,quote=F)
q()

#bgzip
cd /data1/gy/ATACseq_RWAS/caQTL_STITCH/exp
bgzip 95sample_IterativeOverlapPeakSet.final.TMM_qnorm.bed
tabix 95sample_IterativeOverlapPeakSet.final.TMM_qnorm.bed.gz

#################################
##2.genotype PCA for individual##
#################################
conda activate /Public/gaoyun/miniconda3/envs/ATACseq_ASE
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
#####
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_STITCH/genotype_pca_from_qtltools
cd /data1/gy/ATACseq_RWAS/caQTL_STITCH/genotype_pca_from_qtltools
/data1/gy/software/qtltools_v1.2/QTLtools_1.2_source/bin/QTLtools pca \
--vcf /data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.vcf.gz \
--scale --center \
--maf 0.05 --distance 50000 \
--out 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.50kb_dis.95percent
######PC1-PC5
head -n 6 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.50kb_dis.95percent.pca > 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.50kb_dis.first_5percent.pca

#############################
##3.PEER factor calculation##
##############################
###
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/covar
/Public/gaoyun/software/R-4.2.0/bin/R 
library(openxlsx)
library(data.table)
baseline=read.xlsx("/data1/gy/ATACseq_RWAS/ATACseq/baseline/95sample_baseline.xlsx")
rownames(baseline)=baseline[,1]
baseline=baseline[,c(3,2,4)]
baseline_trans=as.data.frame(t(baseline))
baseline_trans$SampleID=rownames(baseline_trans)
baseline_trans=baseline_trans[,c(96,1:95)]
geno_PCA=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/genotype_pca_from_qtltools/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.50kb_dis.first_5percent.pca")
covar=rbind(baseline_trans,geno_PCA)
rownames(covar)=covar$SampleID
covar_trans=t(covar[,-1])
write.table(covar,"/data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/covar/95sample_sex_age_Hp_PC1_5.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(cbind(SampleID=rownames(covar_trans),covar_trans),"/data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/covar/95sample_sex_age_Hp_PC1_5.transpose.txt",sep="\t",col.names=T,row.names=F,quote=F)
q()

###
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/exp
/Public/gaoyun/software/R-4.2.0/bin/R 
library(data.table)
TMM=read.delim("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.TMM_qnorm",header=T)
TMM=TMM[,-c(2:6)]
write.table(TMM,"/data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/exp/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.txt",sep = "\t",quote = F,col.names = T,row.names = F)
q()

###
##https://github.com/PMBio/peer/wiki/Tutorial
##https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R
##https://www.biostars.org/p/9463502/
conda activate /Public/gaoyun/miniconda3/envs/peer
#######
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/peer_output_TMM_qnorm
for peer_num in $(seq 1 15); do
  echo "Running PEER with ${peer_num} factors..."
  Rscript \
  /data1/gy/code/R_script/run_PEER.R \
  /data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/exp/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.txt \
  --covariates /data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/covar/95sample_sex_age_Hp_PC1_5.transpose.txt \
  95sample_IterativeOverlapPeakSet.${peer_num}peer.final.TMM_qnorm \
  ${peer_num} \
  --output_dir /data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/peer_output_TMM_qnorm
done

###################
##4.caQTL mapping##
###################
#####################################################################################################################
###########i.using covariates including sex_age_Hp_PC1_5_peer1_15 for caQTL mapping(1MB permutation-based)###########
#####################################################################################################################
conda activate /Public/gaoyun/miniconda3/envs/ATACseq_ASE
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
##################23.sex_age_Hp_PC1_5_peer1_15
covar=23.95sample_sex_age_Hp_PC1_5_peer1_15.covar
model_type=23.sex_age_Hp_PC1_5_peer1_15
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/${model_type}/1MB_permutation
/data1/gy/software/qtltools_v1.2/QTLtools_1.2_source/bin/QTLtools cis \
--permute 1000 \
--vcf /data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.vcf.gz \
--bed /data1/gy/ATACseq_RWAS/caQTL_STITCH/exp/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.bed.gz \
--cov /data1/gy/ATACseq_RWAS/caQTL_STITCH/multiple_model_covar_all/TMM/${covar} \
--seed 123 \
--window 1000000 \
--out /data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/${model_type}/1MB_permutation/95sample_caQTL_1MB_permutation_result.txt
################permutation result
cd /data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/23.sex_age_Hp_PC1_5_peer1_15/1MB_permutation
echo "pheno_id pheno_chr pheno_start pheno_end pheno_strand n_proximal_var pheno_var_dist var_id var_chr var_start var_end p_degree_freedom dummy beta_dist_par1 beta_dist_par2 p_nominal beta p_empirical p_adjust_beta_dist" \
| sed s/" "/"\t"/g | gzip -c > 95sample_caQTL_1MB_permutation_result.tsv.gz
cat 95sample_caQTL_1MB_permutation_result.txt | sed s/" "/"\t"/g | gzip -c >> 95sample_caQTL_1MB_permutation_result.tsv.gz

#########################################################
###########ii.1MB permutation result (FDR<0.1)###########
#########################################################
cd /data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/23.sex_age_Hp_PC1_5_peer1_15/1MB_permutation
####I.FDR adjustment of permutation-based p-values for caQTLs
/Public/gaoyun/software/R-4.2.0/bin/Rscript /data1/gy/code/cwas_script/qtltools-runFDR_cis.R \
95sample_caQTL_1MB_permutation_result.tsv.gz \
0.1 \
95sample_caQTL_1MB_permutation_result

####II.check the beta approximation used by qtltools
/Public/gaoyun/software/R-4.2.0/bin/Rscript /data1/gy/code/cwas_script/qtltools-check_beta_approx.R \
95sample_caQTL_1MB_permutation_result.tsv.gz \
95sample_caQTL_1MB_permutation_result

##################################################################################################################
###########iii.Statistics of the distances (kb) between top SNPs and OCR center (1MB permutation-based)###########
##################################################################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
t=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/23.sex_age_Hp_PC1_5_peer1_15/1MB_permutation/95sample_caQTL_1MB_permutation_result-significant.tsv.gz")

####peak center
t$mid=(t$pheno_start+t$pheno_end)/2
##########distances between top SNPs and OCR center (kb)
t$dist_center=(t$var_start-t$mid)/1000

##########################Statistics of the distances between top SNPs and OCR center (kb)
t$inpeak = t$pheno_var_dist==0
message("in peak: ", sum(t$inpeak)/nrow(t))
#in peak: 0.165191740412979
t$tenkbpeak = abs(t$dist_center) < 10
message("within 10kb: ", sum(t$tenkbpeak)/nrow(t))
#within 10kb: 0.557522123893805
t$twentyfivekbpeak = abs(t$dist_center) < 25
message("within 25kb: ", sum(t$twentyfivekbpeak)/nrow(t))
#within 25kb: 0.709070796460177
t$onehundredkbpeak = abs(t$dist_center) < 100
message("within 100kb: ", sum(t$onehundredkbpeak)/nrow(t))
#within 100kb: 0.870575221238938
t$twohundredkbpeak = abs(t$dist_center) < 200
message("within 200kb: ", sum(t$twohundredkbpeak)/nrow(t))
#within 200kb: 0.918510324483776

# More than 90% of top SNPs were located within 200 kb of the OCR center,
# so a 200 kb window was used for subsequent nominal caQTL mapping.
q()

############################################################################################################################################
###########iv.Loop over covariate sets including fixed sex_age_Hp_PC1_5 and 0 to 15 PEER factors for caQTL mapping(200kb nominal)###########
############################################################################################################################################
########directory
base_dir=/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp
##
mkdir -p ${base_dir}/8.sex_age_Hp_PC1_5/200kb_nominal
##
for i in $(seq 1 15); do
  dir_num=$((8 + i))
  peer_suffix=$(seq -s "_" 1 ${i})
  mkdir -p ${base_dir}/${dir_num}.sex_age_Hp_PC1_5_peer${peer_suffix}/200kb_nominal
done

#########covar
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_STITCH/multiple_model_covar_all/TMM
cd /data1/gy/ATACseq_RWAS/caQTL_STITCH/multiple_model_covar_all
##
cp /data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/covar/95sample_sex_age_Hp_PC1_5.txt \
  ./TMM/8.95sample_sex_age_Hp_PC1_5.covar
##
for peer_num in $(seq 1 15); do
  out_num=$((8 + peer_num))
  cp /data1/gy/ATACseq_RWAS/caQTL_STITCH/peer/peer_output_TMM_qnorm/95sample_IterativeOverlapPeakSet.${peer_num}peer.final.TMM_qnorm.PEER_covariates.txt \
    ./TMM/${out_num}.95sample_sex_age_Hp_PC1_5_peer1_${peer_num}.covar
done

#########caQTL mapping using QTLtools 
conda activate /Public/gaoyun/miniconda3/envs/ATACseq_ASE
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
###
vcf=/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.vcf.gz
bed=/data1/gy/ATACseq_RWAS/caQTL_STITCH/exp/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.bed.gz
cov_dir=/data1/gy/ATACseq_RWAS/caQTL_STITCH/multiple_model_covar_all/TMM
out_root=/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp
qtltools=/data1/gy/software/qtltools_v1.2/QTLtools_1.2_source/bin/QTLtools
# sex_age_Hp_PC1_5 + 0 to 15 PEER factors
for peer_num in $(seq 0 15); do
  if [ "$peer_num" -eq 0 ]; then
    idx=8
    covar=8.95sample_sex_age_Hp_PC1_5.covar
    model_type=8.sex_age_Hp_PC1_5
  elif [ "$peer_num" -eq 1 ]; then
    idx=9
    covar=9.95sample_sex_age_Hp_PC1_5_peer1.covar
    model_type=9.sex_age_Hp_PC1_5_peer1
  else
    idx=$((8 + peer_num))
    covar=${idx}.95sample_sex_age_Hp_PC1_5_peer1_${peer_num}.covar
    model_type=${idx}.sex_age_Hp_PC1_5_peer1_${peer_num}
  fi

  mkdir -p ${out_root}/${model_type}/200kb_nominal

  echo "Running ${model_type} ..."

  ${qtltools} cis \
    --vcf ${vcf} \
    --bed ${bed} \
    --cov ${cov_dir}/${covar} \
    --nominal 1 \
    --window 200000 \
    --out ${out_root}/${model_type}/200kb_nominal/95sample_caQTL_200kb_nominal_result.txt
done

########################################################################################################################################
###########v.Perform FDR correction and summarize statistics for 200 kb caQTL results across different covariate combinations###########
########################################################################################################################################
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/all_model_result_statistics
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
#
model_types <- c("8.sex_age_Hp_PC1_5", "9.sex_age_Hp_PC1_5_peer1",
                 "10.sex_age_Hp_PC1_5_peer1_2", "11.sex_age_Hp_PC1_5_peer1_3", "12.sex_age_Hp_PC1_5_peer1_4",
                 "13.sex_age_Hp_PC1_5_peer1_5", "14.sex_age_Hp_PC1_5_peer1_6", "15.sex_age_Hp_PC1_5_peer1_7",
                 "16.sex_age_Hp_PC1_5_peer1_8", "17.sex_age_Hp_PC1_5_peer1_9", "18.sex_age_Hp_PC1_5_peer1_10",
                 "19.sex_age_Hp_PC1_5_peer1_11", "20.sex_age_Hp_PC1_5_peer1_12", "21.sex_age_Hp_PC1_5_peer1_13",
                 "22.sex_age_Hp_PC1_5_peer1_14", "23.sex_age_Hp_PC1_5_peer1_15")

#
statistic_list <- list()
#
for (model_type in model_types) {
  ###
  caQTL <- fread(paste0("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/", model_type, "/200kb_nominal/95sample_caQTL_200kb_nominal_result.txt"), header = FALSE)
  colnames(caQTL)=c("pheno_id","pheno_chr","pheno_start","pheno_end","pheno_strand","n_proximal_var","pheno_var_dist","var_id","var_chr","var_start","var_end","p_nominal","beta","top_proximal_var")
  caQTL$FDR=p.adjust(caQTL$p_nominal,"BH")
  fwrite(caQTL,paste0("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/", model_type, "/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt"),sep="\t",col.names=T,row.names=F,quote=F)
  caQTL_0.1=subset(caQTL,caQTL$FDR<0.1)
  fwrite(caQTL_0.1,paste0("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/", model_type, "/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt"),sep="\t",col.names=T,row.names=F,quote=F)
  caQTL_0.1_dis0=subset(caQTL_0.1,caQTL_0.1$pheno_var_dist==0)
  caQTL_0.1_dis1kb=subset(caQTL_0.1,abs(caQTL_0.1$pheno_var_dist)<=1000)
  ###
  num_sigSNPpeak_pair=length(caQTL_0.1$var_id)
  num_sigSNPpeak_pair_dis1kb=length(caQTL_0.1_dis1kb$var_id)
  num_sigSNPpeak_pair_dis0=length(caQTL_0.1_dis0$var_id)
  num_uniquesigSNP=length(unique(caQTL_0.1$var_id))
  num_uniquesigSNP_dis1kb=length(unique(caQTL_0.1_dis1kb$var_id))
  num_uniquesigSNP_dis0=length(unique(caQTL_0.1_dis0$var_id))
  num_uniquesigpeak=length(unique(caQTL_0.1$pheno_id))
  num_uniquesigpeak_dis1kb=length(unique(caQTL_0.1_dis1kb$pheno_id))
  num_uniquesigpeak_dis0=length(unique(caQTL_0.1_dis0$pheno_id))
  ratio_sigSNPpeak_pair_dis1kb_all=num_sigSNPpeak_pair_dis1kb/num_sigSNPpeak_pair
  ratio_sigSNPpeak_pair_dis0_all=num_sigSNPpeak_pair_dis0/num_sigSNPpeak_pair
  ratio_uniquesigSNP_dis1kb_all=num_uniquesigSNP_dis1kb/num_uniquesigSNP
  ratio_uniquesigSNP_dis0_all=num_uniquesigSNP_dis0/num_uniquesigSNP
  ratio_uniquesigpeak_dis1kb_all=num_uniquesigpeak_dis1kb/num_uniquesigpeak
  ratio_uniquesigpeak_dis0_all=num_uniquesigpeak_dis0/num_uniquesigpeak
  #
  statistic_list[[model_type]]=list(model_type=model_type,
            num_sigSNPpeak_pair=num_sigSNPpeak_pair,num_sigSNPpeak_pair_dis1kb=num_sigSNPpeak_pair_dis1kb,num_sigSNPpeak_pair_dis0=num_sigSNPpeak_pair_dis0,
            num_uniquesigSNP=num_uniquesigSNP,num_uniquesigSNP_dis1kb=num_uniquesigSNP_dis1kb,num_uniquesigSNP_dis0=num_uniquesigSNP_dis0,
            num_uniquesigpeak=num_uniquesigpeak,num_uniquesigpeak_dis1kb=num_uniquesigpeak_dis1kb,num_uniquesigpeak_dis0=num_uniquesigpeak_dis0,
            ratio_sigSNPpeak_pair_dis1kb_all=ratio_sigSNPpeak_pair_dis1kb_all,ratio_sigSNPpeak_pair_dis0_all=ratio_sigSNPpeak_pair_dis0_all,
            ratio_uniquesigSNP_dis1kb_all=ratio_uniquesigSNP_dis1kb_all,ratio_uniquesigSNP_dis0_all=ratio_uniquesigSNP_dis0_all,
            ratio_uniquesigpeak_dis1kb_all=ratio_uniquesigpeak_dis1kb_all,ratio_uniquesigpeak_dis0_all=ratio_uniquesigpeak_dis0_all)
}

##
statistic_list_final <- do.call(rbind, statistic_list)
fwrite(statistic_list_final,"/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/all_model_result_statistics/allmodel_type_95sample_caQTL_nominal_200kb_result.statistics",sep = "\t", col.names = T, row.names = FALSE, quote = FALSE)
q()

# The covariate combination including 7 PEER factors yielded the largest number of caOCRs, (see FigureS6C and code of FigureS6C-E.sh)
# so this model was used for downstream analyses.
