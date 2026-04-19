##############################################################
############caQTL-epithelium_sceQTL colocalization############
##############################################################
type=epi
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/$type/input
##################################################################
######1.Extract SNPs shared by the sceQTL and caQTL datasets######
##################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(tidyfst)
#########################################eQTL
type="epi"
##eqtl
eQTL_input=paste0("/data1/gy/sceQTL/eQTL_CHN100K/",type,"/",type,"/",type,"_sceQTL_result.allpairs.withsymbol.withFDR.txt")
eqtl=fread(eQTL_input)
eqtl=eqtl[,c("variant_id","maf")]
eqtl=unique(eqtl)
dim(eqtl)
#[1] 4594724       2
##dosage
dosage_input=paste0("/data1/gy/sceQTL/genotype/",type,"/",type,"_CHN100K_filter05.dosage.txt")
dosage=fread(dosage_input)
print(paste0("sample num:",ncol(dosage)-9))
#[1] "sample num:203"
dosage=dosage[,1:5]
dosage$chrbp=paste0("chr",dosage$'#CHROM',":",dosage$POS)
dosage$ID=paste(dosage$chrbp,dosage$REF,dosage$ALT,sep=":")
dosage=subset(dosage,nchar(REF)==1 & nchar(ALT)==1)
dosage=dosage[,3:6]
colnames(dosage)[1]=c("eqtlID")
##merge
eqtl=merge(dosage,eqtl,by.x="eqtlID",by.y="variant_id")
colnames(eqtl)=c("eqtlID","Allele1","Allele2","chrbp","MAF_eqtl")
dim(eqtl)
#[1] 4578958 
eqtl = eqtl[,c('chrbp','Allele1','Allele2','MAF_eqtl','eqtlID')]

###########################################caQTL freq
freq = fread('/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/plink/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.frq')
freq$caqtlID = freq$SNP
freq[, c("chr", "bp", "a1", "a2") := tstrsplit(SNP, ":", fixed = TRUE)]  
colnames(freq)[5]="MAF_caqtl"
freq = freq[,c('chr','bp','a1','a2','MAF_caqtl','caqtlID')]
###
freq$chrbp = paste(freq$chr,freq$bp,sep=':')
freq = freq[,c('chrbp','a1','a2','MAF_caqtl','caqtlID')]

###########################################merge
eqtl = merge(eqtl,freq,by='chrbp')
dim(eqtl)
#[1] 4039080       9
#######
eqtl1 = eqtl[which((eqtl$Allele1 == eqtl$a1 & eqtl$Allele2 == eqtl$a2) | (eqtl$Allele1 == eqtl$a2 & eqtl$Allele2 == eqtl$a1)),]
dim(eqtl1)
#[1] 4038968       9
#######
eqtl2 = eqtl[!which((eqtl$Allele1 == eqtl$a1 & eqtl$Allele2 == eqtl$a2) | (eqtl$Allele1 == eqtl$a2 & eqtl$Allele2 == eqtl$a1)),]
dim(eqtl2)
#[1] 112   9
#######
eqtl2$nchar1 = nchar(eqtl2$Allele1) + nchar(eqtl2$Allele2)
eqtl2$nchar2 = nchar(eqtl2$a1) + nchar(eqtl2$a2)
eqtl2 = eqtl2[which(eqtl2$nchar1 == eqtl2$nchar2),]
eqtl2_1 = eqtl2[which(eqtl2$nchar1 == 2),]
eqtl2_2 = eqtl2[which(eqtl2$nchar1 > 2),]
eqtl2_1 = eqtl2_1[which(eqtl2_1$Allele1 != eqtl2_1$a1 & eqtl2_1$Allele2 != eqtl2_1$a2 & 
                        eqtl2_1$Allele1 != eqtl2_1$a2 & eqtl2_1$Allele2 != eqtl2_1$a1),]
eqtl2_1 = eqtl2_1[,c(1:9)] 
dim(eqtl2_1 )
#[1] 0 9
#######
eqtl = rbind(eqtl1,eqtl2_1)
dim(eqtl)
#[1] 4038968       9
snplist_outdir=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/caQTL_sceQTL_intersected_SNPlist.txt")
fwrite(eqtl,snplist_outdir,col.names=T,row.names=F,quote=F)

################################################################################
######2.Prepare the caQTL and sceQTL coloc input files for the shared SNPs######
################################################################################
snplist=fread(snplist_outdir)
##caqtl
caqtl=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
caqtl=subset(caqtl,caqtl$var_id %in% snplist$caqtlID)
snplist_sub=snplist[,c("MAF_caqtl","caqtlID")]
caqtl=merge(caqtl,snplist_sub,by.x="var_id",by.y="caqtlID")
caqtl_outdir=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_sceQTL_coloc.txt")
fwrite(caqtl,caqtl_outdir,col.names=T,row.names=F,quote=F)
##eqtl
eQTL_input=paste0("/data1/gy/sceQTL/eQTL_CHN100K/",type,"/",type,"/",type,"_sceQTL_result.allpairs.withsymbol.withFDR.txt")
eqtl=fread(eQTL_input)
eqtl=subset(eqtl,eqtl$variant_id %in% snplist$eqtlID)
eqtl_outdir=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/",type,"_sceQTL_result.allpairs.withsymbol.withFDR.for_caQTL_sceQTL_coloc.txt")
fwrite(eqtl[,1:11],eqtl_outdir,col.names=T,row.names=F,quote=F)
q()

##################################################################################
###3.Prepare the caQTL and sceQTL coloc input files per chr for the shared SNPs###
##################################################################################
######caQTL
type="epi"
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/$type/input/caQTL_input_per_chr
cd /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/$type/input/caQTL_input_per_chr
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
type="epi"
caqtl_input=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_sceQTL_coloc.txt")
caqtl=fread(caqtl_input)
#
for (i in 1:22) {  
  #
  subset_data <- subset(caqtl, pheno_chr == paste0("chr",i))  
  #
  file_name <- paste0("95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_sceQTL_coloc.chr", i,".txt")
  #
  fwrite(subset_data, file_name,sep="\t",col.names=T,row.names=F,quote=F)  
}  
q()

######eQTL
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/$type/input/eQTL_input_per_chr
cd /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/$type/input/eQTL_input_per_chr
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
type="epi"
eqtl_input=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/",type,"_sceQTL_result.allpairs.withsymbol.withFDR.for_caQTL_sceQTL_coloc.txt")
eqtl=fread(eqtl_input)
eqtl[, c("chr", "bp") := tstrsplit(variant_id, ":", fixed = TRUE, keep = 1:2)]  
#
for (i in 1:22) {  
  #
  subset_data <- subset(eqtl, chr == paste0("chr",i))
  #
  file_name <- paste0(type,"_sceQTL_result.allpairs.withsymbol.withFDR.for_caQTL_sceQTL_coloc.chr", i,".txt")
  #
  fwrite(subset_data, file_name,sep="\t",col.names=T,row.names=F,quote=F)  
}  
q()

###########################################################################################################
###4.Generate OCR-gene pairs for colocalization analysis (caOCR boundaries ±200kb within eGene TSS ±1Mb)###
###########################################################################################################
type="epi"
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/$type/input/peak_gene_pair_per_chr
################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
############caqtl
caqtl=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
caqtl=subset(caqtl,caqtl$FDR<0.1)
caqtl=dplyr::select(caqtl,c(pheno_id,pheno_chr,pheno_start,pheno_end))
caqtl=unique(caqtl)
dim(caqtl)
#[1] 13665     4

#############eqtl 
type="epi"
eqtl_input=paste0("/data1/gy/sceQTL/eQTL_CHN100K/",type,"/",type,"/",type,"_sceQTL_result.allpairs.withsymbol.withFDR.FDR0.05.txt")
eqtl=fread(eqtl_input)
eqtl=dplyr::select(eqtl,c(gene_id,gene_name))
eqtl=unique(eqtl)
dim(eqtl)
#[1] 2959    2
gtf <- rtracklayer::import('/data1/gy/public/gtf/gencode.v29lift37.annotation.gtf')
gtf <- as.data.frame(gtf)
gtf=subset(gtf,gtf$type=="gene")
gtf<- dplyr::select(gtf,c(seqnames,start,end,strand,gene_id))#,gene_biotype
gtf<-unique(gtf)
gtf$TSS=gtf$start
gtf$TSS[gtf$strand=="-"]=gtf$end[gtf$strand=="-"]
gtf<-unique(gtf)
gtf$gene_id <- substr(gtf$gene_id,1,15)
gtf<- dplyr::select(gtf,c(seqnames,TSS,gene_id))
eqtl <- dplyr::left_join(eqtl, gtf, by = c("gene_id"="gene_id")) 
dim(eqtl)
#[1] 2959    4

##################
library(dplyr)  
result_list <- lapply(1:nrow(eqtl), function(i) {  
  eqtl_row <- eqtl[i,]  
  #
  matching_caqtl <- caqtl %>%  
  filter(  
    pheno_chr == as.character(eqtl_row$seqnames),  
    (pheno_start - 200000 >= (eqtl_row$TSS - 1000000)) & 
      (pheno_end + 200000 <= (eqtl_row$TSS + 1000000))
  )  
  #
  if (nrow(matching_caqtl) > 0) {  
    return(bind_cols(eqtl_row, matching_caqtl))  
  } else {  
    return(NULL)  
  }  
}) 
#
filtered_eqtl_caqtl <- do.call(bind_rows, result_list[!sapply(result_list, is.null)])  
dim(filtered_eqtl_caqtl)
#[1] 54321     8
###
type="epi"
output=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/peak_gene_pair.for_caQTL_sceQTL_coloc.txt")
fwrite(filtered_eqtl_caqtl, output,sep="\t",col.names=T,row.names=F,quote=F)

#############per chr
for (i in 1:22) {  
  #
  subset_data <- subset(filtered_eqtl_caqtl, pheno_chr == paste0("chr",i))  
  #
  outdir=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/peak_gene_pair_per_chr/")
  file_name <- paste0(outdir,"peak_gene_pair.for_caQTL_sceQTL_coloc.chr", i,".txt")
  #
  fwrite(subset_data, file_name,sep="\t",col.names=T,row.names=F,quote=F)  
}  
q()

################################################
###5.caQTL-epi_sceQTL colocalization analysis###
################################################
type="epi"
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/$type/coloc_result/per_chr
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(coloc)
type="epi"
#############
for (chr_num in 1:22) {  
chr <- paste0("chr", chr_num) 
####
inputdir=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/")
caqtl_file=paste0(inputdir,"caQTL_input_per_chr/","95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_sceQTL_coloc.",chr,".txt")
eqtl_file=paste0(inputdir,"eQTL_input_per_chr/",type,"_sceQTL_result.allpairs.withsymbol.withFDR.for_caQTL_sceQTL_coloc.",chr,".txt")
peak_gene_file=paste0(inputdir,"peak_gene_pair_per_chr/","peak_gene_pair.for_caQTL_sceQTL_coloc.",chr,".txt")
##
caqtl=as.data.table(fread(caqtl_file))
caqtl=caqtl[,c("var_id","pheno_id","p_nominal")]
colnames(caqtl)=c("caqtlID","pheno_id","Pcaqtl")
##
eqtl=as.data.table(fread(eqtl_file))
eqtl=eqtl[,c("variant_id","gene_id","pval_nominal")]
colnames(eqtl)=c("eqtlID","gene_id","Peqtl")
##
peak_gene=as.data.table(fread(peak_gene_file))
#########
snplist_input=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/input/caQTL_sceQTL_intersected_SNPlist.txt")
snplist=data.table(fread(snplist_input))
snplist=snplist[,c("eqtlID","caqtlID","MAF_eqtl","MAF_caqtl")]
caqtl=merge(caqtl,snplist,by="caqtlID")
#
library(dplyr)
peak_gene <- peak_gene %>%  
  mutate(  
    # 计算距离  
    distance = case_when(  
      TSS >= pheno_start & TSS <= pheno_end ~ 0,
      TSS < pheno_start ~ pheno_start - TSS,
      TSS > pheno_end ~  pheno_end - TSS
  )  
)  
#
smy_list <- list()  
# 
for (i in 1:nrow(peak_gene)) {  
  #
  print(paste("Processing Peak-Gene index:", i))  
  #
  current_peakid <- peak_gene$pheno_id[i]  
  current_geneid <- peak_gene$gene_id[i]  
  current_genename <- peak_gene$gene_name[i]  
  current_chr <- peak_gene$pheno_chr[i]  
  current_dis <- peak_gene$distance[i] 
  
  #
  caqtl_summary <- caqtl[pheno_id == current_peakid]  
  #
  eqtl_summary <- eqtl[gene_id == current_geneid]  
  #
  summary <- merge(caqtl_summary, eqtl_summary, by='eqtlID')
  #
  if (nrow(summary) > 0) {  
    result <- coloc.abf(  
      dataset1 = list(snp = summary$caqtlID, pvalues = summary$Peqtl, type = "quant", N = 203, MAF = summary$MAF_eqtl),  
      dataset2 = list(snp = summary$caqtlID, pvalues = summary$Pcaqtl, type = "quant", N = 95, MAF = summary$MAF_caqtl)  
    )  
    #
    max_snp_row <- result$results[which.max(result$results$SNP.PP.H4), ]
    max_snp <- max_snp_row$snp
    max_snp.PPH4 <- max_snp_row$SNP.PP.H4
    #
    colocout <- data.frame(  
      Peak = current_peakid,
      Geneid = current_geneid,
      Symbol = current_genename,
      Distance = current_dis,
      PPH4 = result$summary[['PP.H4.abf']],  
      PPH3 = result$summary[['PP.H3.abf']],  
      PPH2 = result$summary[['PP.H2.abf']],  
      PPH1 = result$summary[['PP.H1.abf']],  
      HitSNP = length(unique(summary$caqtlID)),  
      mineqtlP = min(summary$Peqtl),
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
outdir=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/coloc_result/per_chr/")
write.csv(smy_all_df,file=paste0(outdir,'caQTL_sceQTL_coloc_200kb_1MB_PPH4.',chr,'.csv'),quote=F,row.names=F)
}
#####################all chr
dir <- paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/coloc_result/per_chr/")
csv_files <- list.files(path = dir, pattern = "\\.csv$", full.names = TRUE)  
#
all_data <- do.call(rbind, lapply(csv_files, function(file) {  
  read.csv(file, stringsAsFactors = FALSE)  
}))  
##
dim(all_data)
#[1] 54295    13
all_data_0.5=subset(all_data,all_data$PPH4>=0.5)
##########colocalized caOCR-epi_eGene pair
dim(unique(all_data_0.5))
#[1] 676  13
##########colocalized caOCR num
length(unique(all_data_0.5$Peak))
#[1] 577
##########colocalized epi_eGene num
length(unique(all_data_0.5$Geneid))
#[1] 377
outdir=paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",type,"/coloc_result/")
write.csv(all_data, paste0(outdir,"caQTL_seQTL_coloc_200kb_1MB_PPH4.chr1_22.csv"), row.names = F,quote=F) 
write.csv(all_data_0.5, paste0(outdir,"caQTL_seQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv"), row.names = F,quote=F) 
q()

#########################################################################################
########6.Retain caOCR-epi_eGene pairs colocalized with epithelium-specific OCRs#########
#########################################################################################
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(GenomicRanges)  
############OCR cell-type anno
OCR_annotated=fread("/data1/gy/ATACseq_RWAS/ATACseq/OCR_celltype_anno/macs2_type_peak_anno/95sample_IterativeOverlapPeakSet.type_anno.txt")
OCR_single_type <- OCR_annotated %>%
  filter(type_specific != "" & type_specific != FALSE)
OCR_single_type_sub <- OCR_single_type %>%
  filter(type == "Epithelium")

#############epithelum-specific caOCR
##epithelum-specific OCR num
length(OCR_single_type_sub$OCR)
#[1] 18031
caOCR=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.bed")
OCR_single_type_sub_caOCR=subset(OCR_single_type_sub,OCR %in% caOCR$V4)
##epithelum-specific OCR num
length(OCR_single_type_sub_caOCR$OCR)
#[1] 2825

###############
sceQTL_coloc_path <- paste0(
  "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc/",
  "epi",
  "/coloc_result/caQTL_sceQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv"
)
sceQTL_coloc0.5 <- read.csv(sceQTL_coloc_path)
##retain epithelium-specific OCR
sceQTL_coloc0.5_OCR_single_type_sub <- OCR_single_type_sub %>%
  left_join(sceQTL_coloc0.5, by = c("OCR" = "Peak")) %>%
  rename(Peak = OCR)

###############merge with bulk coloc
##caQTL-bulk_eQTL coloc
eQTL_coloc0.5 <- read.csv(
  "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv"
)
eQTL_coloc0.5$Geneid <- substr(eQTL_coloc0.5$Geneid, 1, 15)
eQTL_coloc0.5_OCR_single_type_sub <- OCR_single_type_sub %>%
  left_join(eQTL_coloc0.5, by = c("OCR" = "Peak")) %>%
  rename(Peak = OCR)
##merge
merged_data <- full_join(
  sceQTL_coloc0.5_OCR_single_type_sub,
  eQTL_coloc0.5_OCR_single_type_sub,
  by = c("Peak", "Geneid"),
  suffix = c("_sceQTL", "_eQTL")
)

#############Epithelium sceQTL coloc num for epithelium-specific OCR
sceQTL_coloc_num <- sum(!is.na(merged_data$Distance_sceQTL))
message(paste0("Epithelium sceQTL coloc num = ", sceQTL_coloc_num))
#Epithelium sceQTL coloc num = 149

#############
output_path <- paste0(
  "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/",
  "epi_caQTL_sceQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv"
)
fwrite(
  merged_data,
  file = output_path,
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

####################add COLOC top SNP rsid
all_data_0.5=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/epi_caQTL_sceQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv")
all_data_0.5[, with_bulk_coloc := ifelse(is.na(PPH4_eQTL), "No", "Yes")]
all_data_0.5=na.omit(all_data_0.5[,c(1:19,39)])
dim(all_data_0.5)
#[1] 149  20
##
EAS=fread("/data1/gy/public/1kg/plinkfile_EAS_rsid/chr_all_1kgv3_2015_EAS.bim")
EAS <- EAS[
  nchar(V5) == 1L & nchar(V6) == 1L &
  V5 %chin% c("A","C","G","T") & V6 %chin% c("A","C","G","T") &
  V5 != V6
]
#
EAS1 <- copy(EAS)[, .(
  id = paste0("chr", V1,":",V4,":",V5,":",V6),
  rsid  = as.character(V2)
)]
EAS2 <- copy(EAS)[, .(
  id = paste0("chr", V1,":",V4,":",V6,":",V5),
  rsid  = as.character(V2)
)]
EAS_both <- rbind(EAS1, EAS2)
setkey(EAS_both, id)
##
coloc <- copy(all_data_0.5)
coloc[EAS_both, rsid := i.rsid, on = .(coloc_topSNP_caqtlID_sceQTL = id)]
all_data_0.5_with_rsid <- coloc
write.csv(all_data_0.5_with_rsid, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/epi_caQTL_sceQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.withbulkanno.with_rsID.csv", row.names = F,quote=F) 
q()
