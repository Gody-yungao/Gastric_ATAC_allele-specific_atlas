#####################################################
############caQTL-bulk_eQTL colocalization###########
#####################################################
###########################################
####1.95sample SNP frequence calculation###
###########################################
cd /data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis
mkdir plink
###
plink \
--vcf 95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.vcf.gz \
--make-bed --keep-allele-order \
--out ./plink/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr \
--double-id 

###freq
plink \
-bfile ./plink/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr \
-freq \
-out ./plink/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr

#######################
####2.bulk eQTL file###
#############/3########
bcftools query -f '%ID\t%REF\t%ALT\n' /data1/RNAseq/262samples/262geno.all.chr.Rsq03.MAF01.G05.M05.Hwe6e.vcf.gz > /data1/gy/262_eQTL_final_result/262geno.all.chr.Rsq03.MAF01.G05.M05.Hwe6e.snp_info.txt

#########
/Public/gaoyun/software/R-4.2.0/bin/R
###################GTF
gtf <- rtracklayer::import('/data1/gy/public/gtf/gencode.v29lift37.annotation.gtf')
gtf <- as.data.frame(gtf)
gtf<- dplyr::select(gtf,c(gene_name,gene_id))#,gene_biotype
gtf<-unique(gtf)

###################eQTL
library(tidyfst) 
library(fst) 
data <- read_fst("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.all.chr.fst")  
library(data.table)
data=as.data.table(data)
###################symbol
library(dplyr)  
#
gtf <- gtf %>%  
  mutate(gene_id = sub("_.*$", "", gene_id))
#
data <- data %>%  
  left_join(gtf, by = c("phenotype_id" = "gene_id"))
###################SNP REF&ALT
snp_info=fread("/data1/gy/262_eQTL_final_result/262geno.all.chr.Rsq03.MAF01.G05.M05.Hwe6e.snp_info.txt",header=F) 
colnames(snp_info)=c("SNPid","REF","ALT")
data <- data %>%  
  left_join(snp_info, by = c("variant_id" = "SNPid")) 
data=data[,c(1,2,11,12,3:10)]
####################
fwrite(data,"/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.with_symbol.all.chr.txt",col.names=T,row.names=F,quote=F,sep="\t")
data_maf0.01=subset(data,data$af<=0.99 & data$af>=0.01)
data_maf0.01$FDR=p.adjust(data_maf0.01$pval_nominal,method="BH")
fwrite(data_maf0.01,"/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR.with_symbol.all.chr.txt",col.names=T,row.names=F,quote=F,sep="\t")
data_maf0.01_FDR0.1=subset(data_maf0.01,data_maf0.01$FDR<0.1)
dim(data_maf0.01_FDR0.1)
#[1] 1127771      13
fwrite(data_maf0.01_FDR0.1,"/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR0.1.with_symbol.all.chr.txt",col.names=T,row.names=F,quote=F,sep="\t")
data_maf0.01_FDR0.05=subset(data_maf0.01,data_maf0.01$FDR<0.05)
dim(data_maf0.01_FDR0.05)
#[1] 882465     13
length(unique(data_maf0.01_FDR0.05$phenotype_id)) ##eGene
#[1] 11232
fwrite(data_maf0.01_FDR0.05,"/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR0.05.with_symbol.all.chr.txt",col.names=T,row.names=F,quote=F,sep="\t")
q()

################################################################
######3.Extract SNPs shared by the eQTL and caQTL datasets######
################################################################
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/per_chr
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
#########################################eQTL
eqtl=fread("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR.with_symbol.all.chr.txt")
eqtl=eqtl[,c("variant_id","af")]
eqtl=unique(eqtl)
dim(eqtl)
#[1] 6543340       2
eqtl[, c("chr", "bp", "Allele1", "Allele2") := tstrsplit(variant_id, ":", fixed = TRUE)]  
eqtl$chrbp=paste0("chr",eqtl$chr,":",eqtl$bp)
colnames(eqtl)[1:2]=c("eqtlID","MAF_eqtl")
eqtl=subset(eqtl,nchar(Allele1)==1 & nchar(Allele2)==1)
dim(eqtl)
#[1] 6009223       7
eqtl = eqtl[,c('chrbp','Allele1','Allele2','MAF_eqtl','eqtlID')]

###########################################caQTL freq
freq = fread('/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/plink/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.frq')
freq$caqtlID = freq$SNP
freq[, c("chr", "bp", "a1", "a2") := tstrsplit(SNP, ":", fixed = TRUE)]  
colnames(freq)[5]="MAF_caqtl"
freq = freq[,c('chr','bp','a1','a2','MAF_caqtl','caqtlID')]
###用chrbp代表SNPid
freq$chrbp = paste(freq$chr,freq$bp,sep=':')
freq = freq[,c('chrbp','a1','a2','MAF_caqtl','caqtlID')]

###########################################merge
eqtl = merge(eqtl,freq,by='chrbp')
dim(eqtl)
#[1] 4377200       9
#######
eqtl1 = eqtl[which((eqtl$Allele1 == eqtl$a1 & eqtl$Allele2 == eqtl$a2) | (eqtl$Allele1 == eqtl$a2 & eqtl$Allele2 == eqtl$a1)),]
dim(eqtl1)
#[1] 4258633       9
#######
eqtl2 = eqtl[!which((eqtl$Allele1 == eqtl$a1 & eqtl$Allele2 == eqtl$a2) | (eqtl$Allele1 == eqtl$a2 & eqtl$Allele2 == eqtl$a1)),]
dim(eqtl2)
#[1] 118567      9
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
#[1] 118468      9

#######
eqtl = rbind(eqtl1,eqtl2_1)
dim(eqtl)
#[1] 4377101       9
fwrite(eqtl,"/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/caQTL_eQTL_intersected_SNPlist.txt",col.names=T,row.names=F,quote=F)
q()

##############################################################################
######4.Prepare the caQTL and eQTL coloc input files for the shared SNPs######
##############################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
snplist=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/caQTL_eQTL_intersected_SNPlist.txt")
##caqtl
caqtl=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
caqtl=subset(caqtl,caqtl$var_id %in% snplist$caqtlID)
snplist_sub=snplist[,c("MAF_caqtl","caqtlID")]
caqtl=merge(caqtl,snplist_sub,by.x="var_id",by.y="caqtlID")
fwrite(caqtl,"/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_eQTL_coloc.txt",col.names=T,row.names=F,quote=F)
##eqtl
eqtl=fread("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR.with_symbol.all.chr.txt")
eqtl=subset(eqtl,eqtl$variant_id %in% snplist$eqtlID)
fwrite(eqtl,"/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/262sample_ciseQTL_1MB_result.cis_qtl_pairs.chr1_22.with_symbol.for_caQTL_eQTL_coloc.txt",col.names=T,row.names=F,quote=F)
q()


################################################################################
###5.Prepare the caQTL and eQTL coloc input files per chr for the shared SNPs###
################################################################################
######caQTL
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/caQTL_input_per_chr
cd /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/caQTL_input_per_chr
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
caqtl=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_eQTL_coloc.txt")
#
for (i in 1:22) {  
  #
  subset_data <- subset(caqtl, pheno_chr == paste0("chr",i))  
  #
  file_name <- paste0("95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_eQTL_coloc.chr", i,".txt")
  #
  fwrite(subset_data, file_name,sep="\t",col.names=T,row.names=F,quote=F)  
}  
q()

######eQTL
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/eQTL_input_per_chr
cd /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/eQTL_input_per_chr
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
eqtl=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/262sample_ciseQTL_1MB_result.cis_qtl_pairs.chr1_22.with_symbol.for_caQTL_eQTL_coloc.txt")
#
eqtl[, c("chr", "bp", "Allele1", "Allele2") := tstrsplit(variant_id, ":", fixed = TRUE)]  
eqtl=eqtl[,-c("Allele1", "Allele2")]
#
for (i in 1:22) {  
  #
  subset_data <- subset(eqtl, chr == i)  
  #
  file_name <- paste0("262sample_ciseQTL_1MB_result.cis_qtl_pairs.chr1_22.with_symbol.for_caQTL_eQTL_coloc.chr", i,".txt")
  #
  fwrite(subset_data, file_name,sep="\t",col.names=T,row.names=F,quote=F)  
}  
q()

###########################################################################################################
###6.Generate OCR-gene pairs for colocalization analysis (caOCR boundaries ±200kb within eGene TSS ±1Mb)###
###########################################################################################################
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/peak_gene_pair_per_chr

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
eqtl=fread("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR.with_symbol.all.chr.txt")
eqtl=subset(eqtl,FDR<0.05)
eqtl=dplyr::select(eqtl,c(phenotype_id,gene_name))
eqtl=unique(eqtl)
dim(eqtl)
#[1] 11232     2
gtf <- rtracklayer::import('/data1/gy/public/gtf/gencode.v29lift37.annotation.gtf')
gtf <- as.data.frame(gtf)
gtf=subset(gtf,gtf$type=="gene")
gtf<- dplyr::select(gtf,c(seqnames,start,end,strand,gene_id))#,gene_biotype
gtf<-unique(gtf)
gtf$TSS=gtf$start
gtf$TSS[gtf$strand=="-"]=gtf$end[gtf$strand=="-"]
gtf<-unique(gtf)
gtf$gene_id <- sub("_.*", "", gtf$gene_id)
gtf<- dplyr::select(gtf,c(seqnames,TSS,gene_id))
eqtl <- dplyr::left_join(eqtl, gtf, by = c("phenotype_id"="gene_id")) 
dim(eqtl)
#[1] 11232     4

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
#[1] 196911      8
fwrite(filtered_eqtl_caqtl, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/peak_gene_pair.for_caQTL_eQTL_coloc.txt",sep="\t",col.names=T,row.names=F,quote=F)

#############per chr
dir.create("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/peak_gene_pair_per_chr")
for (i in 1:22) {  
  #
  subset_data <- subset(filtered_eqtl_caqtl, pheno_chr == paste0("chr",i))  
  #
  outdir="/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/peak_gene_pair_per_chr/"
  file_name <- paste0(outdir,"peak_gene_pair.for_caQTL_eQTL_coloc.chr", i,".txt")
  #
  fwrite(subset_data, file_name,sep="\t",col.names=T,row.names=F,quote=F)  
}  
q()

##########################################
###7.caQTL-eQTL colocalization analysis###
##########################################
mkdir -p /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/per_chr
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(coloc)
#############
for (chr_num in 1:22) {  
chr <- paste0("chr", chr_num) 

####
inputdir="/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/"
caqtl_file=paste0(inputdir,"caQTL_input_per_chr/","95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_eQTL_coloc.",chr,".txt")
eqtl_file=paste0(inputdir,"eQTL_input_per_chr/","262sample_ciseQTL_1MB_result.cis_qtl_pairs.chr1_22.with_symbol.for_caQTL_eQTL_coloc.",chr,".txt")
peak_gene_file=paste0(inputdir,"peak_gene_pair_per_chr/","peak_gene_pair.for_caQTL_eQTL_coloc.",chr,".txt")
##
caqtl=as.data.table(fread(caqtl_file))
caqtl=caqtl[,c("var_id","pheno_id","p_nominal")]
colnames(caqtl)=c("caqtlID","pheno_id","Pcaqtl")
##
eqtl=as.data.table(fread(eqtl_file))
eqtl=eqtl[,c("variant_id","phenotype_id","pval_nominal")]
colnames(eqtl)=c("eqtlID","phenotype_id","Peqtl")
##
peak_gene=as.data.table(fread(peak_gene_file))
####
snplist=data.table(fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/caQTL_eQTL_intersected_SNPlist.txt"))
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

##################
smy_list <- list()  
#
for (i in 1:nrow(peak_gene)) {  
  #
  current_peakid <- peak_gene$pheno_id[i]  
  current_geneid <- peak_gene$phenotype_id[i]  
  current_genename <- peak_gene$gene_name[i]  
  current_chr <- peak_gene$pheno_chr[i]  
  current_dis <- peak_gene$distance[i] 
  
  #
  caqtl_summary <- caqtl[pheno_id == current_peakid]  
  #
  eqtl_summary <- eqtl[phenotype_id == current_geneid]  
  #
  summary <- merge(caqtl_summary, eqtl_summary, by='eqtlID')
  #
  if (nrow(summary) > 0) {  
    result <- coloc.abf(  
      dataset1 = list(snp = summary$caqtlID, pvalues = summary$Peqtl, type = "quant", N = 262, MAF = summary$MAF_eqtl),  
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
outdir="/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/per_chr/"
write.csv(smy_all_df,file=paste0(outdir,'caQTL_eQTL_coloc_200kb_1MB_PPH4.',chr,'.csv'),quote=F,row.names=F)
}

#####################all chr
dir <- "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/per_chr"  
csv_files <- list.files(path = dir, pattern = "\\.csv$", full.names = TRUE)  
#
all_data <- do.call(rbind, lapply(csv_files, function(file) {  
  read.csv(file, stringsAsFactors = FALSE)  
}))  
##
dim(all_data)
#[1] 196740     13
all_data_0.5=subset(all_data,all_data$PPH4>=0.5)

##########colocalized caOCR-eGene pair
dim(unique(all_data_0.5))
#[1] 2053   13
##########colocalized caOCR num
length(unique(all_data_0.5$Peak))
#[1] 1350
##########colocalized eGene num
length(unique(all_data_0.5$Geneid))
#[1] 1205

##
write.csv(all_data, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4.chr1_22.csv", row.names = F,quote=F) 
write.csv(all_data_0.5, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv", row.names = F,quote=F) 

####################add COLOC top SNP rsid
all_data_0.5=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv")

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
coloc[EAS_both, rsid := i.rsid, on = .(coloc_topSNP_caqtlID = id)]
all_data_0.5_with_rsid <- coloc
write.csv(all_data_0.5_with_rsid, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.with_rsID.csv", row.names = F,quote=F) 
q()

