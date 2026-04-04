######################################################
###########Figure3G-H&FigureS7C&FigureS7E#############
######################################################
##############################################################################################
##########1.Check if all SNP reference alleles correspond to the hg19 forward strand##########
##############################################################################################
conda activate RWAS
bcftools norm \
--check-ref -s \
-f /data1/gy/public/genome/hg19/hg19.fa \
-Oz \
/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.vcf.gz \
-o /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.norm.vcf.gz
#Lines   total/split/realigned/skipped:	4993520/0/0/0
#REF/ALT total/modified/added:  	4993520/0/0
rm /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.norm.vcf.gz

#############################################################
################2.caSNP_inpeak,asSNP bed file################
#############################################################
conda activate motifbreakR
R
library(data.table)
##
caSNP_inpeak=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.vcf")
asSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_ASCA_FDR0.1.SNP.vcf")
caSNP_inpeak$start=caSNP_inpeak$POS-1
caSNP_inpeak$V5=0
caSNP_inpeak$V6="+"
asSNP$start=asSNP$POS-1
asSNP$V5=0
asSNP$V6="+"
caSNP_inpeak=caSNP_inpeak[,c(1,6,2,3,7,8)]
asSNP=asSNP[,c(1,6,2,3,7,8)]

####
write.table(caSNP_inpeak,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(asSNP,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_ASCA_FDR0.1.SNP.bed",sep="\t",quote=F,col.names=F,row.names=F)

###########################################################################################################
################3.control SNP set（SNP_inpeak_with_nocaQTL,SNP_inpeak_with_noASCA bed file)################
###########################################################################################################
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
ASCA=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR")
caQTL_inpeak_nosig=subset(caQTL,pheno_var_dist==0 & FDR>=0.1)
ASCA_nosig=subset(ASCA,C0.BBINOM.FDR>=0.1)
caQTL_inpeak_nosig$start=caQTL_inpeak_nosig$var_start-1
caQTL_inpeak_nosig$V5=0
caQTL_inpeak_nosig$V6="+"
caQTL_inpeak_nosig=caQTL_inpeak_nosig[,c(9,16,10,8,17,18)]
ASCA_nosig$start=ASCA_nosig$POS-1
ASCA_nosig$V5=0
ASCA_nosig$V6="+"
ASCA_nosig=ASCA_nosig[,c(1,24,2,3,25,26)]

##
write.table(caQTL_inpeak_nosig,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_caQTL_nosig.SNP.in_corresponding_peak.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(ASCA_nosig,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_ASCA_nosig.SNP.bed",sep="\t",quote=F,col.names=F,row.names=F)
q()

########################################################################################################################
################4.motifbreakR anno and calculate TF disrupt score for caSNP/asSNP (in corresponding OCR)################
########################################################################################################################
####motifbreakR format
####name is defined as chromosome:start:REF:ALT
R
library(motifbreakR)
library(BSgenome)
####
available.SNPs()
available.genomes()
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
caSNP_inpeak.snps.mb <- snps.from.file("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.bed",
 search.genome = BSgenome.Hsapiens.UCSC.hg19,
 format = "bed")
asSNP.snps.mb <- snps.from.file("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_ASCA_FDR0.1.SNP.bed",
 search.genome = BSgenome.Hsapiens.UCSC.hg19,
 format = "bed")
 
####
library(MotifDb)
MotifDb
table(mcols(MotifDb)$organism, mcols(MotifDb)$dataSource)
#
# Step 1: Filter for human (Hsapiens) motifs
hsapiens_motifs <- subset(MotifDb, organism == "Hsapiens")
# Step 2: Remove motifs from the 'stamlab' source
filtered_motifs <- hsapiens_motifs[!grepl("stamlab", names(hsapiens_motifs))]

#
jaspar2016_motifs <- MotifDb[grep("jaspar2016", values(MotifDb)$dataSource)]  
#
jaspar2016_hsapiens <- jaspar2016_motifs[grep("Hsapiens", values(jaspar2016_motifs)$organism)]  
#
jaspar2016_hsapiens  

#
jaspar2018_motifs <- MotifDb[grep("jaspar2018", values(MotifDb)$dataSource)]  
#
jaspar2018_hsapiens <- jaspar2018_motifs[grep("Hsapiens", values(jaspar2018_motifs)$organism)]  
#
jaspar2018_hsapiens  

#
jaspar2022_motifs <- MotifDb[grep("jaspar2022", values(MotifDb)$dataSource)]  
#
jaspar2022_hsapiens <- jaspar2022_motifs[grep("Hsapiens", values(jaspar2022_motifs)$organism)]  
#
jaspar2022_hsapiens  

###############################################################################combind
##
#
all_raw_names <- c(names(jaspar2016_hsapiens), 
                   names(jaspar2018_hsapiens), 
                   names(jaspar2022_hsapiens))

#
temp_names <- sub("^[^-]+-[^-]+-", "", all_raw_names)

#
clean_names <- sub("-MA[0-9]+\\.[0-9]+$", "", temp_names)

#
final_tfs <- unique(toupper(clean_names))


#
tfs_no_parentheses <- sub("\\(.*\\)", "", final_tfs)

#
tfs_no_parentheses[tfs_no_parentheses == "T"] <- "TBXT"

#
final_clean_tfs <- unique(tfs_no_parentheses)

# 
cat("Final Unique TF num:", length(final_clean_tfs), "\n")
#########################Final Unique TF num: 676

####
caSNP_inpeak_results <- motifbreakR(snpList = caSNP_inpeak.snps.mb, filterp = TRUE,
                       pwmList = filtered_motifs ,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10))
asSNP_results <- motifbreakR(snpList = asSNP.snps.mb, filterp = TRUE,
                       pwmList = filtered_motifs,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 20))          
####
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/total_database")
names(caSNP_inpeak_results) <- NULL   
caSNP_inpeak_results_df=as.data.frame(caSNP_inpeak_results)
names(asSNP_results) <- NULL   
asSNP_results_df=as.data.frame(asSNP_results)
#
caSNP_inpeak_results_df$motifPos <- sapply(caSNP_inpeak_results_df$motifPos, function(x) {  
  if (is.null(x)) {  
    return(NA) 
  } else {  
    paste(x, collapse = ",") 
  }  
}) 
asSNP_results_df$motifPos <- sapply(asSNP_results_df$motifPos, function(x) {  
  if (is.null(x)) {  
    return(NA) 
  } else {  
    paste(x, collapse = ",")
  }  
}) 
#  
library(openxlsx)
write.xlsx(caSNP_inpeak_results_df[caSNP_inpeak_results_df$effect == "strong",],"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/total_database/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
write.xlsx(asSNP_results_df[asSNP_results_df$effect == "strong",],"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/total_database/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.xlsx")
q()

################################################################################################################################
################5.motifbreakR anno and calculate TF disrupt score for non-caSNP/non-asSNP (in corresponding OCR)################
################################################################################################################################
####
caSNP_inpeak_nosig_results <- motifbreakR(snpList = caSNP_inpeak_nosig.snps.mb, filterp = TRUE,
                       pwmList = filtered_motifs ,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 20))
asSNP_nosig_results <- motifbreakR(snpList = asSNP_nosig.snps.mb, filterp = TRUE,
                       pwmList = filtered_motifs,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 20))          
####
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/total_database")
names(caSNP_inpeak_nosig_results) <- NULL   
caSNP_inpeak_nosig_results_df=as.data.frame(caSNP_inpeak_nosig_results)
names(asSNP_nosig_results) <- NULL   
asSNP_nosig_results_df=as.data.frame(asSNP_nosig_results)
# 
caSNP_inpeak_nosig_results_df$motifPos <- sapply(caSNP_inpeak_nosig_results_df$motifPos, function(x) {  
  if (is.null(x)) {  
    return(NA) 
  } else {  
    paste(x, collapse = ",") 
  }  
}) 
asSNP_nosig_results_df$motifPos <- sapply(asSNP_nosig_results_df$motifPos, function(x) {  
  if (is.null(x)) {  
    return(NA) 
  } else {  
    paste(x, collapse = ",") 
  }  
}) 
#  
library(openxlsx)
write.xlsx(caSNP_inpeak_nosig_results_df[caSNP_inpeak_nosig_results_df$effect == "strong",],"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/total_database/95sample_caQTL_nosig.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
write.xlsx(asSNP_nosig_results_df[asSNP_nosig_results_df$effect == "strong",],"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/total_database/95sample_ASCA_nosig.SNP.motifbreakR_anno.strong.xlsx")
q()

#################################################################
################6.keep jaspar2016,2018,2022 motif################
#################################################################
R
library(openxlsx)
caSNP_inpeak_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/total_database/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
asSNP_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/total_database/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.xlsx")
caSNP_inpeak_nosig_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/total_database/95sample_caQTL_nosig.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
asSNP_nosig_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/total_database/95sample_ASCA_nosig.SNP.motifbreakR_anno.strong.xlsx")
##
library(dplyr)
caSNP_inpeak_results_df <- caSNP_inpeak_results_df %>%
  filter(grepl("jaspar20", dataSource, ignore.case = TRUE))
asSNP_results_df <- asSNP_results_df %>%
  filter(grepl("jaspar20", dataSource, ignore.case = TRUE))
caSNP_inpeak_nosig_results_df <- caSNP_inpeak_nosig_results_df %>%
  filter(grepl("jaspar20", dataSource, ignore.case = TRUE))
asSNP_nosig_results_df <- asSNP_nosig_results_df %>%
  filter(grepl("jaspar20", dataSource, ignore.case = TRUE))
#
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar")
write.xlsx(caSNP_inpeak_results_df,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
write.xlsx(asSNP_results_df,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.xlsx")
#
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar")
write.xlsx(caSNP_inpeak_nosig_results_df,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar/95sample_caQTL_nosig.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
write.xlsx(asSNP_nosig_results_df,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar/95sample_ASCA_nosig.SNP.motifbreakR_anno.strong.xlsx")
q()

######################################################################
################7.keep TFs expressed in 262sample RNAseq################
######################################################################
R
library(openxlsx)
caSNP_inpeak_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
asSNP_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.xlsx")
caSNP_inpeak_nosig_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar/95sample_caQTL_nosig.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
asSNP_nosig_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar/95sample_ASCA_nosig.SNP.motifbreakR_anno.strong.xlsx")

#
clean_gene_symbol <- function(x) {
  #
  x <- gsub("[\\(（][^\\)）]*[\\)）]", "", x)
  #
  x <- gsub("\\s+", "", x)
  return(x)
}

#
caSNP_inpeak_results_df$Symbol <- clean_gene_symbol(caSNP_inpeak_results_df$geneSymbol)
asSNP_results_df$Symbol <- clean_gene_symbol(asSNP_results_df$geneSymbol)
caSNP_inpeak_nosig_results_df$Symbol <- clean_gene_symbol(caSNP_inpeak_nosig_results_df$geneSymbol)
asSNP_nosig_results_df$Symbol <- clean_gene_symbol(asSNP_nosig_results_df$geneSymbol)

##
caSNP_inpeak_results_df$Symbol <- replace(
  caSNP_inpeak_results_df$Symbol,
  caSNP_inpeak_results_df$Symbol == "Pax6",
  "PAX6"
)
asSNP_results_df$Symbol <- replace(
  asSNP_results_df$Symbol,
  asSNP_results_df$Symbol == "Pax6",
  "PAX6"
)
caSNP_inpeak_nosig_results_df$Symbol <- replace(
  caSNP_inpeak_nosig_results_df$Symbol,
  caSNP_inpeak_nosig_results_df$Symbol == "Pax6",
  "PAX6"
)
asSNP_nosig_results_df$Symbol <- replace(
  asSNP_nosig_results_df$Symbol,
  asSNP_nosig_results_df$Symbol == "Pax6",
  "PAX6"
)

##
caSNP_inpeak_results_df$Symbol <- replace(
  caSNP_inpeak_results_df$Symbol,
  caSNP_inpeak_results_df$Symbol == "T",
  "TBXT"
)
asSNP_results_df$Symbol <- replace(
  asSNP_results_df$Symbol,
  asSNP_results_df$Symbol == "T",
  "TBXT"
)
caSNP_inpeak_nosig_results_df$Symbol <- replace(
  caSNP_inpeak_nosig_results_df$Symbol,
  caSNP_inpeak_nosig_results_df$Symbol == "T",
  "TBXT"
)
asSNP_nosig_results_df$Symbol <- replace(
  asSNP_nosig_results_df$Symbol,
  asSNP_nosig_results_df$Symbol == "T",
  "TBXT"
)

##
TF_list=unique(c(caSNP_inpeak_results_df$Symbol,asSNP_results_df$Symbol,caSNP_inpeak_nosig_results_df$Symbol,asSNP_nosig_results_df$Symbol))
length(TF_list)
#[1] 623

##
gtf <- rtracklayer::import('/data1/gy/public/gtf/gencode.v29lift37.annotation.gtf')
gtf <- as.data.frame(gtf)
gtf<- dplyr::select(gtf,c(gene_name,gene_id))#,gene_biotype
gtf<-unique(gtf)
gtf$gene_id=substr(gtf$gene_id,1,15)
##
TF_list1 <- TF_list[TF_list %in% gtf$gene_name] ##TF complex
TF_list2 <- TF_list[!TF_list %in% gtf$gene_name] ##non-TF complex

##keep
library(data.table)
exp=fread("/data1/RNAseq/262samples/eQTL/262sample_Stomach.expression.bed.gz")
exp$gene_id=substr(exp$gene_id,1,15)
gtf_exp=subset(gtf,gene_id %in% exp$gene_id)

####
TF_list1_exp <- TF_list1[TF_list1 %in% unique(gtf_exp$gene_name)]
#
split_tf <- function(x) unlist(strsplit(x, split = "[:\\-]{2}|-"))
#
keep_idx <- sapply(TF_list2, function(tfs) {
  genes <- unlist(strsplit(tfs, split = "[:]{2}|-"))
  all(genes %in% gtf_exp$gene_name)
})
#
TF_list2_exp <- TF_list2[keep_idx]
##
TF_list_exp=c(TF_list1_exp,TF_list2_exp)
length(TF_list_exp)
#474


##
caSNP_inpeak_results_df_exp = subset(caSNP_inpeak_results_df,Symbol %in% TF_list_exp)
asSNP_results_df_exp = subset(asSNP_results_df,Symbol %in% TF_list_exp)
caSNP_inpeak_nosig_results_df_exp = subset(caSNP_inpeak_nosig_results_df,Symbol %in% TF_list_exp)
asSNP_nosig_results_df_exp = subset(asSNP_nosig_results_df,Symbol %in% TF_list_exp)

####
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar_262exp")
write.xlsx(caSNP_inpeak_results_df_exp,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar_262exp/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
write.xlsx(asSNP_results_df_exp,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar_262exp/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.xlsx")
#
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar_262exp")
write.xlsx(caSNP_inpeak_nosig_results_df_exp,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar_262exp/95sample_caQTL_nosig.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
write.xlsx(asSNP_nosig_results_df_exp,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar_262exp/95sample_ASCA_nosig.SNP.motifbreakR_anno.strong.xlsx")
q()

#######################################################################################################################
################8.Enrichment results (per TF) for caSNP_inpeak and asSNP from motifbreakR TF annotation################
#######################################################################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(openxlsx)
#####
caSNP_inpeak=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.bed")
asSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_ASCA_FDR0.1.SNP.bed")
caSNP_inpeak_nosig=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_caQTL_nosig.SNP.in_corresponding_peak.bed")
asSNP_nosig=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_ASCA_nosig.SNP.bed")
#####
caSNP_inpeak_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar_262exp/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
asSNP_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar_262exp/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.xlsx")
caSNP_inpeak_nosig_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar_262exp/95sample_caQTL_nosig.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
asSNP_nosig_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output_control/jaspar_262exp/95sample_ASCA_nosig.SNP.motifbreakR_anno.strong.xlsx")
#####
length(unique(caSNP_inpeak_results_df$Symbol))
#[1] 473
length(unique(asSNP_results_df$Symbol))
#[1] 474
length(unique(caSNP_inpeak_nosig_results_df$Symbol))
#[1] 474
length(unique(asSNP_nosig_results_df$Symbol))
#[1] 474
TF_list=unique(c(caSNP_inpeak_results_df$Symbol,asSNP_results_df$Symbol,caSNP_inpeak_nosig_results_df$Symbol,asSNP_nosig_results_df$Symbol))
length(TF_list)
#[1] 474
##
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/enrichment_result")
#####################################################
####################i.caSNP_inpeak###################
#####################################################
#
results <- data.frame(FEATURE = character(),  
                      Odds_Ratio = numeric(),  
                      P_Value = numeric(),  
                      CI_Lower = numeric(),  
                      CI_Upper = numeric(), 
                      SNP_num = numeric(), 
                      stringsAsFactors = FALSE)  

#
for (i in 1:length(TF_list)) {  
  TF <- TF_list[i]  
  
  #
  TF_snps_caSNP_inpeak <- subset(caSNP_inpeak_results_df, geneSymbol == TF)  
  TF_snps_caSNP_inpeak_count <- length(unique(TF_snps_caSNP_inpeak$SNP_id))  
  TF_snps_caSNP_inpeak_nosig <- subset(caSNP_inpeak_nosig_results_df, geneSymbol == TF)  
  TF_snps_caSNP_inpeak_nosig_count <- length(unique(TF_snps_caSNP_inpeak_nosig$SNP_id))  

  #
  a <- TF_snps_caSNP_inpeak_count      
  b <- nrow(caSNP_inpeak) - TF_snps_caSNP_inpeak_count   
  c <- TF_snps_caSNP_inpeak_nosig_count   
  d <- nrow(caSNP_inpeak_nosig) - TF_snps_caSNP_inpeak_nosig_count  

  #
  if (a > 0 && b > 0 && c > 0 && d > 0) {  
    #
    contingency_table <- matrix(c(a, b, c, d), nrow = 2)  

    #
    fisher_test <- fisher.test(contingency_table, alternative = "two.sided", conf.int = TRUE)  
  
    #
    results <- rbind(results, data.frame(FEATURE = TF,  
                                        Odds_Ratio = fisher_test$estimate,  
                                        P_Value = fisher_test$p.value,  
                                        CI_Lower = fisher_test$conf.int[1],  
                                        CI_Upper = fisher_test$conf.int[2],
                                        SNP_num = TF_snps_caSNP_inpeak_count))  
  } else {  
    #
    results <- rbind(results, data.frame(FEATURE = TF,  
                                        Odds_Ratio = NA,  
                                        P_Value = NA,  
                                        CI_Lower = NA,  
                                        CI_Upper = NA,
                                        SNP_num = NA))
  }    
}
results$FDR=p.adjust(results$P_Value,method="BH")
write.csv(results,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/enrichment_result/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.enrichment_result.csv",row.names=F)

###############################################
####################ii.asSNP###################
###############################################
# 
results <- data.frame(FEATURE = character(),  
                      Odds_Ratio = numeric(),  
                      P_Value = numeric(),  
                      CI_Lower = numeric(),  
                      CI_Upper = numeric(),  
                      SNP_num = numeric(), 
                      stringsAsFactors = FALSE)  

#
for (i in 1:length(TF_list)) {  
  TF <- TF_list[i]  
  
  #
  TF_snps_asSNP <- subset(asSNP_results_df, geneSymbol == TF)  
  TF_snps_asSNP_count <- length(unique(TF_snps_asSNP$SNP_id))  
  TF_snps_asSNP_nosig <- subset(asSNP_nosig_results_df, geneSymbol == TF)  
  TF_snps_asSNP_nosig_count <- length(unique(TF_snps_asSNP_nosig$SNP_id))  

  #
  a <- TF_snps_asSNP_count      
  b <- nrow(asSNP) - TF_snps_asSNP_count   
  c <- TF_snps_asSNP_nosig_count   
  d <- nrow(asSNP_nosig) - TF_snps_asSNP_nosig_count  

  #
  if (a > 0 && b > 0 && c > 0 && d > 0) {  
    #
    contingency_table <- matrix(c(a, b, c, d), nrow = 2)  

    #
    fisher_test <- fisher.test(contingency_table, alternative = "two.sided", conf.int = TRUE)  
  
    #
    results <- rbind(results, data.frame(FEATURE = TF,  
                                        Odds_Ratio = fisher_test$estimate,  
                                        P_Value = fisher_test$p.value,  
                                        CI_Lower = fisher_test$conf.int[1],  
                                        CI_Upper = fisher_test$conf.int[2],
                                        SNP_num = TF_snps_asSNP_count))  
  } else {  
    #
    results <- rbind(results, data.frame(FEATURE = TF,  
                                        Odds_Ratio = NA,  
                                        P_Value = NA,  
                                        CI_Lower = NA,  
                                        CI_Upper = NA,
                                        SNP_num = NA))  
  }    
}
results$FDR=p.adjust(results$P_Value,method="BH")
write.csv(results,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/enrichment_result/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.enrichment_result.csv",row.names=F)

###################################################################
################9.Enrichment result dotplot (Top20)################
###################################################################
library(ggplot2)  
library(dplyr)  
###########################################
#####i.caSNP_inpeak enrichment dotplot#####
###########################################
setwd("/data1/gy/ATAC_for_review/Figure3G-H&FigureS7C&FigureS7E/output/caSNP_inOCR")
results=read.csv("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/enrichment_result/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.enrichment_result.csv")
###
library(dplyr)  
#
results <- results %>%  
  arrange(FDR)  
#
results <- results %>%  
  mutate(FEATURE = factor(FEATURE, levels = FEATURE))  

###
library(ggplot2)
library(viridis)
p1 = ggplot(results[1:20,], aes(x = FEATURE, y = -log10(P_Value))) +
  geom_point(aes(size=Odds_Ratio,color=SNP_num)) +
  scale_color_viridis_c(option = "magma",
                        name = "# of variants\naltering motifs") +  
  labs(x = "",
       y = "-log10(P)") +
  theme_minimal() +
  theme(
    # 
    axis.line.x = element_line(linewidth = 0.3, color = "black"),  
    axis.line.y = element_blank(),                              
    axis.title = element_text(size = 7),
    axis.text.x = element_text(angle = 30, hjust = 1, size=7),
    axis.text.y = element_text(size=7),
    plot.title = element_text(size = 7, hjust = 0.5),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 5),
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank(), 
    #
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray50", linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.06, "cm")
  )
ggsave(p1,
       filename="95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.enrichment_result.top20.dotplot.pdf",
       width=5,height=4) ##Figure3G(top)

######RNAseq expression
exp=fread("/data1/RNAseq/262samples/eQTL/262sample.expression.tpm.txt")
exp=as.data.frame(exp)
rownames(exp)=exp$ID
exp=exp[,-1]
exp=as.data.frame(t(exp))
exp$gene_id=rownames(exp)
gtf <- rtracklayer::import('/data1/gy/public/gtf/gencode.v29lift37.annotation.gtf')
gtf <- as.data.frame(gtf)
gtf<- dplyr::select(gtf,c(gene_name,gene_id))#,gene_biotype
gtf<-unique(gtf)
gtf$gene_id=substr(gtf$gene_id,1,15)
exp=merge(gtf,exp,by="gene_id")
exp$median_expression <- apply(exp[, 3:ncol(exp)], 1, median, na.rm = TRUE)
exp=exp[,c(1,2,265)]

#
results_top20 <- results[1:20, ]
#
results_top20$FEATURE <- as.character(results_top20$FEATURE)
#
results_top20$expression_value <- NA
#
for (i in 1:nrow(results_top20)) {
  feature <- results_top20$FEATURE[i]
  
  #
  if (grepl("::", feature)) {
    # 拆分复合物
    genes <- unlist(strsplit(feature, "::"))
    gene1 <- genes[1]
    gene2 <- genes[2]
    
    #
    exp_value1 <- exp$median_expression[exp$gene_name == gene1]
    exp_value2 <- exp$median_expression[exp$gene_name == gene2]
    
    #
    if (length(exp_value1) == 0) exp_value1 <- NA
    if (length(exp_value2) == 0) exp_value2 <- NA
    
    #
    avg_exp <- mean(c(exp_value1, exp_value2), na.rm = TRUE)
    #
    if (is.nan(avg_exp)) avg_exp <- NA
    
    results_top20$expression_value[i] <- avg_exp
  } else {
    #
    exp_value <- exp$median_expression[exp$gene_name == feature]
    if (length(exp_value) == 0) exp_value <- NA
    results_top20$expression_value[i] <- exp_value
  }
}

################################TF expression heatmap
library(ComplexHeatmap)  
library(circlize)  
#
mat <- results_top20[, c("FEATURE", "expression_value")]
rownames(mat) <- mat$FEATURE  
mat$expression_value=log2(1+mat$expression_value)

#
col_exp <- colorRamp2(c(0, 6), c("white", "#E31A1C"))  

#
ht_exp <- Heatmap(
  t(mat[, "expression_value", drop = FALSE]),  
  name = "Expression",
  col = col_exp,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  width = unit(20, "cm"),
  height = unit(1, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(
  direction = "horizontal",
  legend_width = unit(4, "cm")
  )
)
pdf("95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.enrichment_result.expression.top20.heatmap.pdf",height=3,width=8) ##Figure3G(bottom)
draw(ht_exp, heatmap_legend_side = "bottom")
dev.off()

#####################################
#####ii.asSNP enrichment dotplot#####
#####################################
setwd("/data1/gy/ATAC_for_review/Figure3G-H&FigureS7C&FigureS7E/output/asSNP")
results=read.csv("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/enrichment_result/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.enrichment_result.csv")
###
library(dplyr)  
#
results <- results %>%  
  arrange(FDR)  
#
results <- results %>%  
  mutate(FEATURE = factor(FEATURE, levels = FEATURE))  

###
library(ggplot2)
library(viridis)
p1 = ggplot(results[1:20,], aes(x = FEATURE, y = -log10(P_Value))) +
  geom_point(aes(size=Odds_Ratio,color=SNP_num)) +
  #scale_size_continuous(range = c(1, 7)) +  
  scale_color_viridis_c(option = "magma",
                        limits = c(50,150),
                        oob = scales::squish,
                        name = "# of variants\naltering motifs") +  
  labs(x = "",
       y = "-log10(P)") +
  #scale_y_continuous(
  #  breaks = seq(0, ceiling(max(-log10(results$P_Value))), by = 5) 
  #) +
  scale_y_continuous(breaks = c(16, 20, 24)) +
  theme_minimal() +
  theme(
    #
    axis.line.x = element_line(linewidth = 0.3, color = "black"),
    axis.line.y = element_blank(),          
    axis.title = element_text(size = 7),
    axis.text.x = element_text(angle = 30, hjust = 1, size=7),
    axis.text.y = element_text(size=7),
    plot.title = element_text(size = 7, hjust = 0.5),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 5),
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank(),
    #
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray50", linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.06, "cm")
  )
ggsave(p1,
       filename="95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.enrichment_result.top20.dotplot.pdf",
       width=5,height=4) ##FigureS7C(top)

######RNAseq expression
exp=fread("/data1/RNAseq/262samples/eQTL/262sample.expression.tpm.txt")
exp=as.data.frame(exp)
rownames(exp)=exp$ID
exp=exp[,-1]
exp=as.data.frame(t(exp))
exp$gene_id=rownames(exp)
gtf <- rtracklayer::import('/data1/gy/public/gtf/gencode.v29lift37.annotation.gtf')
gtf <- as.data.frame(gtf)
gtf<- dplyr::select(gtf,c(gene_name,gene_id))#,gene_biotype
gtf<-unique(gtf)
gtf$gene_id=substr(gtf$gene_id,1,15)
exp=merge(gtf,exp,by="gene_id")
exp$median_expression <- apply(exp[, 3:ncol(exp)], 1, median, na.rm = TRUE)
exp=exp[,c(1,2,265)]

#
results_top20 <- results[1:20, ]
#
results_top20$FEATURE <- as.character(results_top20$FEATURE)
#
results_top20$expression_value <- NA
#
for (i in 1:nrow(results_top20)) {
  feature <- results_top20$FEATURE[i]
  
  #
  if (grepl("::", feature)) {
    #
    genes <- unlist(strsplit(feature, "::"))
    gene1 <- genes[1]
    gene2 <- genes[2]
    
    #
    exp_value1 <- exp$median_expression[exp$gene_name == gene1]
    exp_value2 <- exp$median_expression[exp$gene_name == gene2]
    
    #
    if (length(exp_value1) == 0) exp_value1 <- NA
    if (length(exp_value2) == 0) exp_value2 <- NA
    
    #
    avg_exp <- mean(c(exp_value1, exp_value2), na.rm = TRUE)
    #
    if (is.nan(avg_exp)) avg_exp <- NA
    
    results_top20$expression_value[i] <- avg_exp
  } else {
    #
    exp_value <- exp$median_expression[exp$gene_name == feature]
    if (length(exp_value) == 0) exp_value <- NA
    results_top20$expression_value[i] <- exp_value
  }
}

################################TF expression heatmap
library(ComplexHeatmap)  
library(circlize)  
#
mat <- results_top20[, c("FEATURE", "expression_value")]
rownames(mat) <- mat$FEATURE  
mat$expression_value=log2(1+mat$expression_value)

#
col_exp <- colorRamp2(c(0, 6), c("white", "#E31A1C"))  

#
ht_exp <- Heatmap(
  t(mat[, "expression_value", drop = FALSE]),  
  name = "Expression",
  col = col_exp,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  width = unit(20, "cm"),
  height = unit(1, "cm"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(
  direction = "horizontal", 
  legend_width = unit(4, "cm") 
  )
)
pdf("5sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.enrichment_result.expression.top20.heatmap.pdf",height=3,width=8)  ##FigureS7C(bottom)
draw(ht_exp, heatmap_legend_side = "bottom")
dev.off()

#####################################################################################################################
################10.association between motif-disrupting allele relative entropy and caQTL/ASCA effect################
#####################################################################################################################
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect")
library(dplyr)
library(ggplot2)
##############################
########i.caSNP_inpeak########
##############################
setwd("/data1/gy/ATAC_for_review/Figure3G-H&FigureS7C&FigureS7E/output/caSNP_inOCR")
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/scatter_plot_caQTL")
caSNP_inpeak_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar_262exp/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")
caSNP_inpeak_results_df=caSNP_inpeak_results_df[,c(6,7,8,27,29)]
#
caSNP_inpeak_results_df <- aggregate(alleleDiff ~ SNP_id + Symbol + REF + ALT, 
                   data = caSNP_inpeak_results_df, 
                   FUN = function(x) mean(x, na.rm = TRUE))
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/input")
write.xlsx(caSNP_inpeak_results_df,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.xlsx")

##
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt")
caQTL=subset(caQTL,pheno_var_dist==0)
caQTL=caQTL[,c(8,13)]
caSNP_inpeak_results_df=left_join(caSNP_inpeak_results_df,caQTL,by=c("SNP_id"="var_id"))
caSNP_inpeak_results_df$disrupt_diff = -caSNP_inpeak_results_df$alleleDiff
##
caSNP_inpeak_results_df <- caSNP_inpeak_results_df %>%  
  mutate(  
    disrupt_diff_new = ifelse(disrupt_diff < 0, -disrupt_diff, disrupt_diff),  
    beta_new = ifelse(disrupt_diff < 0, -beta, beta) 
  )  
## lm~beta+0
# TFlist
TF_list=unique(caSNP_inpeak_results_df$Symbol)

# 
results <- data.frame(  
  TF = character(), 
  slope = numeric(), 
  p_value = numeric(), 
  stringsAsFactors = FALSE  
)  

# 
for (TF in TF_list) {  
  #
  caSNP_inpeak_results_df_sub <- subset(caSNP_inpeak_results_df, Symbol == TF)  
  
  #
  if (nrow(caSNP_inpeak_results_df_sub) > 1) {  
    #
    lm_result <- tryCatch(  
      lm(beta_new ~ disrupt_diff_new +0, data = caSNP_inpeak_results_df_sub),  
      error = function(e) NA  
    )  
    
    #
    if (inherits(lm_result, "lm")) { 
      slope <- coef(lm_result)[1]
      p_value <- summary(lm_result)$coefficients[1, 4]
      #
      siglabel <- ifelse(p_value >= 0.05, "noPsig",   
                   ifelse(slope < 0, "decreased_Psig", "increased_Psig"))  

      #scatterplot
      y_max <- max(abs(caSNP_inpeak_results_df_sub$beta_new)) * 1.1  

      p <- ggplot(caSNP_inpeak_results_df_sub, aes(x = disrupt_diff_new, y = beta_new)) +  
           geom_point(color = "black", size = 2) +
           geom_abline(  
            slope = slope,   
            intercept = 0,   
            color = ifelse(slope < 0, "#EF6548",  "#4EB3D3"),
            linetype = "solid"  
          ) + 
           geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
           geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
           labs(  
           title = paste("TF:", TF),  
           subtitle = paste0("Slope: ", round(slope, 3), ", P-value: ", format(p_value, scientific = TRUE)),  
           x = "Relative entropy\n(noneffect allele - effect allele)",  
           y = "caQTL effect size"  
          ) +  
           scale_y_continuous(expand = c(0, 0)) +
           coord_cartesian(ylim = c(-y_max, y_max)) +
           theme_minimal() + 
           theme(  
           panel.grid = element_blank(),
           panel.border = element_rect(color = "black", fill = NA, size = 1), 
           axis.ticks = element_line(color = "black"), 
           axis.text = element_text(size = 10, color = "black"),
           plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
           plot.subtitle = element_text(hjust = 0.5, size = 12)
          )
      
      #
      ggsave(filename = paste0("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/scatter_plot_caQTL/scatterplot_", TF,"_",siglabel, ".pdf"),
              plot = p, width = 5.5, height = 6)  ##Figure3H (right)
        
    } else {  
      slope <- NA  
      p_value <- NA  
    }  
  } else {  
    #
    slope <- NA  
    p_value <- NA  
  }  
  
  #
  results <- rbind(results, data.frame(TF = TF, slope = slope, p_value = p_value))  
}  
##FDR
results$FDR=p.adjust(results$p_value,method="BH")
#
results <- results %>%  
  mutate(  
    sig = case_when(  
      # 
      !is.na(FDR) & FDR < 0.05 & slope < 0 ~ "decreased",  
      
      # 
      !is.na(FDR) & FDR < 0.05 & slope > 0 ~ "increased",  
      
      #
      TRUE ~ "non-significant"   
    )  
  )
table(results$sig)
#      decreased       increased non-significant 
#            118               5             350
#
write.csv(results, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.lm_result.csv", row.names = FALSE)


###barplot
#
total_N <- sum(sig_counts$count)
increased_N <- sig_counts$count[2]
decreased_N <- sig_counts$count[3]
sig_N <- increased_N+decreased_N
nonsig_N <- total_N - sig_N 

#
plot_data1 <- data.frame(
  category = c("non-sig", "sig"),
  count = c(nonsig_N, sig_N)
)
plot_data2 <- data.frame(
  category = c("increased", "decreased"),
  count = c(increased_N, decreased_N)
)

#
plot_data1$category <- factor(plot_data1$category, levels = c("non-sig", "sig"))
plot_data2$category <- factor(plot_data2$category, levels = c("increased", "decreased"))

#
p1_1 <- ggplot(plot_data1, aes(x = "caSNP-associated", y = count, fill = category)) + 
  
  #
  geom_col(position = "fill", width = 0.5, alpha = 1) +
  
  #
  geom_text(aes(label = paste0(count, "\n(", percent(count/total_N, accuracy = 0.1), ")")), 
            position = position_fill(vjust = 0.5), 
            size = 5, color = "black") +
  
  #
  scale_fill_manual(values = c("sig" = "#fcbba1", "non-sig" = "grey80")) +
  
  #
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage") +
  
  #
  theme_classic() +
  labs(x = NULL) + 
  
  theme(
    #
    axis.text.x = element_text(size = 12, color = "black"), 
    legend.position = "none"
  )

p1_2 <- ggplot(plot_data2, aes(x = "significant", y = count, fill = category)) + 
  
  #
  geom_col(position = "fill", width = 0.5, alpha = 1) +
  
  #
  geom_text(aes(label = paste0(count, "\n(", percent(count/sig_N, accuracy = 0.1), ")")), 
            position = position_fill(vjust = 0.5), 
            size = 5, color = "black") +
  
  #
  scale_fill_manual(values = c("increased" = "#4eb3d3", "decreased" = "#ef6548")) +
  
  #
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage") +
  
  #
  theme_classic() +
  labs(x = NULL) + 
  
  theme(
    #
    axis.text.x = element_text(size = 12, color = "black"),
    legend.position = "none"
  )

#
library(cowplot)
#
combined_plot1 <- plot_grid(p1_1, p1_2, 
                           align = "h", 
                           axis = "tb", 
                           nrow = 1, 
                           rel_widths = c(1, 1)) 
ggsave(combined_plot1,file="95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.motifbreakR_anno.strong.lm_result.withNA.barplot.pdf",width=5,height=6) ##Figure3H (left)

########################
########ii.asSNP########
########################
setwd("/data1/gy/ATAC_for_review/Figure3G-H&FigureS7C&FigureS7E/output/asSNP")
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/scatter_plot_ASCA")
asSNP_results_df=read.xlsx("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/output/jaspar_262exp/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.xlsx")
asSNP_results_df=asSNP_results_df[,c(6,7,8,27,29)]
#
asSNP_results_df <- aggregate(alleleDiff ~ SNP_id + Symbol + REF + ALT, 
                   data = asSNP_results_df, 
                   FUN = function(x) mean(x, na.rm = TRUE))
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/input")
write.xlsx(asSNP_results_df,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/input/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.xlsx")

##
ASCA=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1")
ASCA$ALT_REF_ratio_log2 <- log2(ASCA$ALL.AF / (1 - ASCA$ALL.AF)) 
ASCA=ASCA[,c(3,24)]
asSNP_results_df=left_join(asSNP_results_df,ASCA,by=c("SNP_id"="RSID"))
asSNP_results_df$disrupt_diff=-asSNP_results_df$alleleDiff
#
asSNP_results_df <- asSNP_results_df %>%  
  mutate(  
    disrupt_diff_new = ifelse(disrupt_diff < 0, -disrupt_diff, disrupt_diff),
    ALT_REF_ratio_log2_new = ifelse(disrupt_diff < 0, -ALT_REF_ratio_log2, ALT_REF_ratio_log2)
  )  
# TFlist
TF_list=unique(asSNP_results_df$Symbol)

#
results <- data.frame(  
  TF = character(),    
  slope = numeric(),    
  p_value = numeric(),  
  stringsAsFactors = FALSE  
)  

#
for (TF in TF_list) {  
  #
  asSNP_results_df_sub <- subset(asSNP_results_df, Symbol == TF)  
  
  #
  if (nrow(asSNP_results_df_sub) > 1) {  
    #
    lm_result <- tryCatch(  
      lm(ALT_REF_ratio_log2_new ~ disrupt_diff_new + 0, data = asSNP_results_df_sub),  
      error = function(e) NA  
    )  
    
    #
    if (inherits(lm_result, "lm")) {
      slope <- coef(lm_result)[1]
      p_value <- summary(lm_result)$coefficients[1, 4]
      #
      siglabel <- ifelse(p_value >= 0.05, "noPsig",   
                   ifelse(slope < 0, "decreased_Psig", "increased_Psig"))  

      #
      #
      y_max <- max(abs(asSNP_results_df_sub$ALT_REF_ratio_log2_new)) * 1.1  

      p <- ggplot(asSNP_results_df_sub, aes(x = disrupt_diff_new, y = ALT_REF_ratio_log2_new)) +  
           geom_point(color = "black", size = 2) + 
           geom_abline(  
            slope = slope,   
            intercept = 0,   
            color = ifelse(slope < 0, "#EF6548",  "#4EB3D3"),
            linetype = "solid"  
          ) + 
           geom_vline(xintercept = 0, color = "black", linetype = "dotted") + 
           geom_hline(yintercept = 0, color = "black", linetype = "dotted") + 
           labs(  
           title = paste("TF:", TF),  
           subtitle = paste0("Slope: ", round(slope, 3), ", P-value: ", format(p_value, scientific = TRUE)),  
           x = "Relative entropy\n(noneffect allele - effect allele)",
           y = expression("ASCA log"[2]*"(ALT/REF)") 
          ) +  
           scale_y_continuous(expand = c(0, 0)) + 
           coord_cartesian(ylim = c(-y_max, y_max)) + 
           theme_minimal() + 
           theme(  
           panel.grid = element_blank(),
           panel.border = element_rect(color = "black", fill = NA, size = 1),
           axis.ticks = element_line(color = "black"), 
           axis.text = element_text(size = 10, color = "black"),
           plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
           plot.subtitle = element_text(hjust = 0.5, size = 12) 
          )
      
      #
      ggsave(filename = paste0("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/scatter_plot_ASCA/scatterplot_", TF,"_",siglabel, ".pdf"),
              plot = p, width = 5.5, height = 6)  ##FigureS7E (right)

    } else {  
      slope <- NA  
      p_value <- NA  
    }  
  } else {  
    #
    slope <- NA  
    p_value <- NA  
  }  
  
  #
  results <- rbind(results, data.frame(TF = TF, slope = slope, p_value = p_value))  
}  
##FDR
results$FDR=p.adjust(results$p_value,method="BH")
#
results <- results %>%  
  mutate(  
    sig = case_when(  
      #
      !is.na(FDR) & FDR < 0.05 & slope < 0 ~ "decreased",  
      
      #
      !is.na(FDR) & FDR < 0.05 & slope > 0 ~ "increased",  
      
      #
      TRUE ~ "non-significant"   
    )  
  )
table(results$sig)
#      decreased       increased non-significant 
#            107              11             356
# 
write.csv(results, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/motif-disrupting_allele_vs_ATAC_effect/95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.lm_result.csv", row.names = FALSE)

###barplot
#
total_N <- sum(sig_counts$count)
increased_N <- sig_counts$count[2]
decreased_N <- sig_counts$count[3]
sig_N <- increased_N+decreased_N
nonsig_N <- total_N - sig_N 

#
plot_data1 <- data.frame(
  category = c("non-sig", "sig"),
  count = c(nonsig_N, sig_N)
)
plot_data2 <- data.frame(
  category = c("increased", "decreased"),
  count = c(increased_N, decreased_N)
)

#
plot_data1$category <- factor(plot_data1$category, levels = c("non-sig", "sig"))
plot_data2$category <- factor(plot_data2$category, levels = c("increased", "decreased"))

#
p2_1 <- ggplot(plot_data1, aes(x = "asSNP-associated", y = count, fill = category)) + 
  
  #
  geom_col(position = "fill", width = 0.5, alpha = 1) +
  
  #
  geom_text(aes(label = paste0(count, "\n(", percent(count/total_N, accuracy = 0.1), ")")), 
            position = position_fill(vjust = 0.5), 
            size = 5, color = "black") +
  
  #
  scale_fill_manual(values = c("sig" = "#fcbba1", "non-sig" = "grey80")) +
  
  #
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage") +
  
  #
  theme_classic() +
  labs(x = NULL) + 
  
  theme(
    #
    axis.text.x = element_text(size = 12, color = "black"), 
    legend.position = "none"
  )
p2_2 <- ggplot(plot_data2, aes(x = "significant", y = count, fill = category)) + 
  
  #
  geom_col(position = "fill", width = 0.5, alpha = 1) +
  
  #
  geom_text(aes(label = paste0(count, "\n(", percent(count/sig_N, accuracy = 0.1), ")")), 
            position = position_fill(vjust = 0.5), 
            size = 5, color = "black") +
  
  #
  scale_fill_manual(values = c("increased" = "#4eb3d3", "decreased" = "#ef6548")) +
  
  #
  scale_y_continuous(labels = scales::percent_format(), name = "Percentage") +
  
  #
  theme_classic() +
  labs(x = NULL) + 
  
  theme(
    #
    axis.text.x = element_text(size = 12, color = "black"),
    legend.position = "none"
  )

##
library(cowplot)
#
combined_plot2 <- plot_grid(p2_1, p2_2, 
                           align = "h", 
                           axis = "tb", 
                           nrow = 1, 
                           rel_widths = c(1, 1)) 
ggsave(combined_plot2,file="95sample_ASCA_FDR0.1.SNP.motifbreakR_anno.strong.lm_result.withNA.barplot.pdf",width=5,height=6) ##FigureS7E (left)
