############################
########ASCA analysis#######
############################
##############################################################################
#########1.Extract SNP sites where at least one sample is heterozygous########
##############################################################################
conda activate RWAS
cd /data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis
bcftools view -i 'COUNT(GT="het")>0' \
95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.vcf.gz \
-o 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.het.for_analysis.withchr.phased.vcf
bgzip 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.het.for_analysis.withchr.phased.vcf
tabix 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.het.for_analysis.withchr.phased.vcf.gz
bcftools view -H 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.vcf.gz|wc -l
#4993520
bcftools view -H 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.het.for_analysis.withchr.phased.vcf.gz|wc -l
#4993520
###with no SNPs removed

##########################################################################
#########2.Add readgroup information to the WASP-filtered BAM file########
##########################################################################
mkdir -p /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/1.bam_withRG
###
cd /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/1.bam_withRG
ls /data1/gy/ATACseq_RWAS/ATACseq/WASP_STITCH/7.dedup_bam/*.WASPfilter.sort.bam > 1
ls /data1/gy/ATACseq_RWAS/ATACseq/WASP_STITCH/7.dedup_bam/*.WASPfilter.sort.bam |cut -d"/" -f 8|cut -d"." -f 1  > 0
paste 0 1  > config_wasp_bam
##
cat config_wasp_bam |while read id;
do echo $id
arr=($id)
bam_input=${arr[1]}
sampleid=${arr[0]}
bam_output=/data1/gy/ATACseq_RWAS/ASCA_STITCH/input/1.bam_withRG/${sampleid}.WASPfilter.red.sort.bam 
java -jar /Public/gaoyun/software/picard.jar AddOrReplaceReadGroups \
 I=$bam_input \
 O=$bam_output \
 RGLB=lib \
 RGPL=illumina \
 RGPU=run \
 RGSM=$sampleid
samtools index $bam_output 
done

######################################
#########3.gatk ASEReadCounter########
######################################
cd /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/1.bam_withRG
cp -r /data1/gy/NO_tissue/ATACseq_ASE/ref /data1/gy/ATACseq_RWAS/ASCA_STITCH/input
mkdir ../2.allele_counts
##
cat config_wasp_bam |while read id;
do echo $id
arr=($id)
sampleid=${arr[0]}
bam_RG=/data1/gy/ATACseq_RWAS/ASCA_STITCH/input/1.bam_withRG/${sampleid}.WASPfilter.red.sort.bam 
vcf=/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.het.for_analysis.withchr.phased.vcf.gz
gatk ASEReadCounter \
   -R /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/ref/hg19.fa \
   -O /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/2.allele_counts/${sampleid}.csv \
   -I $bam_RG \
   -V $vcf \
   -min-depth 1
done

##################################################
#########4.add allele-specific count to vcf########
##################################################
mkdir ../phased_vcf_persample
mkdir ../3.stratas_prep_ase_vcf
cd ../3.stratas_prep_ase_vcf
##
# Extract a single sample's genotype from a multi-sample phased VCF,
# write the sample's ref/alt read counts from a .csv file to the VCF FORMAT/AS field,
# then sort and index the resulting VCF.
cat ../1.bam_withRG/config_wasp_bam|while read id;
do echo $id
arr=($id)
sampleid=${arr[0]}

bcftools view -s ${sampleid} \
/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.het.for_analysis.withchr.phased.vcf.gz | bgzip \
> ../phased_vcf_persample/ATACseq_${sampleid}.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.het.for_analysis.withchr.phased.vcf.gz

# Input (vcf.gz) with phased=genotypes for this individual
VCF=../phased_vcf_persample/ATACseq_${sampleid}.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.het.for_analysis.withchr.phased.vcf.gz
# Input (.csv) with alleles for this individual
CSV=../2.allele_counts/${sampleid}.csv
# Output prefix
OUT=${sampleid}
# Sample identifier to put into the VCF header
SAMPLE=${sampleid}

#################################vcf with count
# Construct the header
zcat $VCF | awk -vs=$SAMPLE 'BEGIN{ OFS="\t" } { if(substr($1,1,1) != "#") exit 1; if($1 == "#CHROM") { print "##FORMAT=<ID=AS,Number=2,Type=Integer,Description=\"Read counts for the ref and alt alleles\">"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,s; } else print $0; }' > $OUT.vcf
# Join the body
join -a1 -1 1 -2 1 \
<(zcat $VCF | awk 'substr($1,1,1) != "#" { printf $1":"$2" "; print $0 }' | sed 's/chr//' | sort -k1,1) \
<(cat $CSV | tail -n+2 | awk '{ printf $1":"$2" "; printf "%d %d\n",$6,$7 }' | sed 's/chr//' | sort -k1,1) \
| cut -d ' ' -f2- \
| awk 'BEGIN{ OFS="\t" } { if (NF != 10) print $1,$2,$3,$4,$5,$6,$7,$8,$9":AS",$10":"$11","$12; else print $1,$2,$3,$4,$5,$6,$7,$8,$9":AS",$10":0,0"; }' \
>> $OUT.vcf
# Sort and index
gatk SortVcf -I $OUT.vcf -O $OUT.sort.vcf 
bgzip -c $OUT.sort.vcf > $OUT.sort.vcf.gz
rm $OUT.sort.vcf
rm $OUT.vcf 
tabix $OUT.sort.vcf.gz
done

##############################combine vcf with allele-specific count of 95 individuals
vcf_files=$(ls *.sort.vcf.gz)
bcftools merge ${vcf_files} | bgzip -c > 95sample_BYES_merged.sort.vcf.gz

##############################################
#########5.Estimating prior parameters########
##############################################
####
cd ../3.stratas_prep_ase_vcf
mkdir ../4.stratas_prep_ase_counts_file
mkdir ../5.stratas_input_params_file
cat ../1.bam_withRG/config_wasp_bam_part_09 |while read id;
do echo $id
arr=($id)
sampleid=${arr[0]}
zcat ${sampleid}.sort.vcf.gz \
| grep -v '#' \
| cut -f 1,2,10 \
| tr ':' '\t' \
| awk '$3 != "1|1" && $3 != "0|0"' \
| tr ',' '\t' \
| sed 's/chr//' \
| awk 'BEGIN { print "CHR POS HAP REF.READS ALT.READS" } $4 + $5 > 0' \
 > ../4.stratas_prep_ase_counts_file/${sampleid}.counts
##
sed '1s/ /\t/g' ../4.stratas_prep_ase_counts_file/${sampleid}.counts | bgzip -c > ../4.stratas_prep_ase_counts_file/${sampleid}.AS.counts.gz

######Estimating prior parameters
###params.R from
/Public/gaoyun/software/R-4.2.0/bin/Rscript /data1/gy/code/R_script/RWAS/params.R \
--inp_counts ../4.stratas_prep_ase_counts_file/${sampleid}.AS.counts.gz \
--out ../5.stratas_input_params_file/${sampleid} \
--min_cov 5 \
--id ${sampleid}
done

###################################combine prior parameters of 95 individuals
cd ../5.stratas_input_params_file
output_file="95sample_BYES.global.params"
touch $output_file
#
for file in *.params
do
    tail -n +2 $file >> $output_file
done
#
head -n 1 BYES02022.global.params > header 
cat header $output_file >> 95sample_BYES_withheader.global.params

#######################################
#########6.peak file for stratas#######
#######################################
cd /data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample
/Public/gaoyun/software/R-4.2.0/bin/R 
library(data.table)
data=as.data.frame(fread("95sample_IterativeOverlapPeakSet.sort.bed"))
data$V5=(data$V2+1+data$V3)/2
data$V4=paste0(data$V1,":",data$V2+1,"-",data$V3)
colnames(data)=c("CHR","P0","P1","NAME","CENTER")
dim(data)
#[1] 109187      5
write.table(data,"95sample_IterativeOverlapPeakSet.for_stratas_input.sort.bed",sep="\t",quote=F,col.names=T,row.names=F)
q()

################################################################
#########7 Testing for population allele specificity########
################################################################
###################################################
#########i.input for allele specificity test#######
###################################################
cd /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/5.stratas_input_params_file
mkdir ../6.stratas_input_MAT_file
mkdir ../baseline
####A vcf file containing all individuals is converted to counts and split into batches as follows:
VCF=../3.stratas_prep_ase_vcf/95sample_BYES_merged.sort.vcf.gz
OUT=../6.stratas_input_MAT_file/95sample_BYES_merged
zcat $VCF\
| grep -v '#' \
| cut -f 1-5,9- \
| awk 'BEGIN{OFS="\t"} {for(i=1;i<=5;i++) printf "%s ",$i; for(i=6;i<=NF;i++) {gsub(":"," ",$i); printf "%s\t",$i}; print ""}' \
| tr ',' '\t' \
| sed 's/[|,]/ /g' \
| tr '\t' ' ' \
| sed 's/GT AS//g' \
| split -d -l 20000 - $OUT.MAT.split.
####sample file
awk 'BEGIN {print "ID\tCONDITION"} {print $1 "\t0"}' ../5.stratas_input_params_file/95sample_BYES.global.params \
> ../baseline/95sample_BYES.AS.PHE

##########################################
#########ii.allele specificity test#######
##########################################
conda activate /Public/gaoyun/miniconda3/envs/ATACseq_ASE
cd ../6.stratas_input_MAT_file
mkdir -p ../7.stratas_output_file/stratas_output_250bp
ls /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/6.stratas_input_MAT_file/95sample_BYES_merged.MAT.split.* > 1
ls /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/6.stratas_input_MAT_file/95sample_BYES_merged.MAT.split.* |cut -d"/" -f 8|cut -d"." -f 4 > 0
paste 0 1 > config_MAT

####ASCA analysis
##--indiv TRUE output individual-level read count information
##--window -1 inpeak
cat config_MAT|while read id;
do echo $id
arr=($id)
i=${arr[0]}
/Public/gaoyun/software/R-4.2.0/bin/Rscript /data1/gy/software/stratAS-master/stratas.R \
--input /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/6.stratas_input_MAT_file/95sample_BYES_merged.MAT.split.${i} \
--seed 123 \
--max_rho 0.2 \
--window -1 \
--min_cov 1 \
--fill_cnv FALSE \
--indiv TRUE \
--samples /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/baseline/95sample_BYES.AS.PHE \
--peaks /data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.for_stratas_input.sort.bed \
--global_param /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/5.stratas_input_params_file/95sample_BYES_withheader.global.params \
> /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/7.stratas_output_file/stratas_output_250bp/95sample_BYES_output_250bp_${i}.profile
done

########combine results
cd /data1/gy/ATACseq_RWAS/ASCA_STITCH/input
mkdir -p ../ASCA_result/region_250bp
cp ./7.stratas_output_file/stratas_output_250bp/$(ls /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/7.stratas_output_file/stratas_output_250bp | head -n 1) ../ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined
for filename in $(ls /data1/gy/ATACseq_RWAS/ASCA_STITCH/input/7.stratas_output_file/stratas_output_250bp | tail -n +2); do
   tail -n +2 ./7.stratas_output_file/stratas_output_250bp/$filename >> ../ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined
done

#######FDR correction
cd ../ASCA_result/region_250bp
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
data=as.data.frame(fread("95sample_stratas_output_files_250bp_combined"))
data$'C0.BBINOM.FDR'=p.adjust(data$'C0.BBINOM.P',"BH")
fwrite(data,"95sample_stratas_output_files_250bp_combined_withFDR",sep="\t",quote=F,col.names=T,row.names=F)
data_0.1=subset(data,data$'C0.BBINOM.FDR'<0.1)
length(unique(data_0.1$NAME))
#[1] 5658
fwrite(data_0.1,"95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1",sep="\t",quote=F,col.names=T,row.names=F)
q()
