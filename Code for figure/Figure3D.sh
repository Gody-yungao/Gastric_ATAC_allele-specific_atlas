##############################
###########Figure3D###########
##############################
#################################################################
#######1.42GWAS from BBJ published in 2020 (PMID:32514122)#######
#################################################################
########################################
#######i.BBJ GWAS GARFIELD format#######
#########################################
cd /data1/gy/public/BBJ_2020_GWAS/raw
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
###
# 
gwas_files <- list.files(pattern = "^bbj.*\\.vcf.gz$") 
#
output_dir <- "/data1/gy/public/BBJ_2020_GWAS/garfield_format/"  
# 
results_summary <- data.frame(
  file = character(),
  significant_count = integer(),
  stringsAsFactors = FALSE
)

# 
for (i in seq_along(gwas_files)) {

  #
  name <- sub("\\.vcf\\.gz$", "", gwas_files[i])
  GWAS <- fread(gwas_files[i], skip = "#CHROM", select = c("#CHROM","POS",name))
  #
  setnames(GWAS, c("chromosome","base_pair_location","SAMPLE"))
  #
  parts <- tstrsplit(GWAS$SAMPLE, ":", fixed = TRUE)
  #
  GWAS[, LP := as.numeric(parts[[3]])]
  #
  GWAS[, p_value := 10^(-LP)]
  #
  GWAS <- GWAS[, .(chromosome, base_pair_location, p_value)]
  
  #
  sig_count <- sum(GWAS$p_value < 5e-8, na.rm = TRUE)

  #
  results_summary <- rbind(results_summary,
                           data.frame(file = name,
                                      significant_count = sig_count,
                                      stringsAsFactors = FALSE))


  #
  dir.create(paste0(output_dir, name), showWarnings = FALSE, recursive = TRUE)

  #
  for (m in 1:22) {
    GWAS_chr <- subset(GWAS, chromosome == m)
    fwrite(GWAS_chr[, c(2,3)],
           file = paste0(output_dir, name, "/chr", m),
           col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
  }

  cat("Processed:", name, " (Found", sig_count, "significant SNPs)\n")
}

#
fwrite(results_summary, file = paste0(output_dir, "GWAS_significant_counts.tsv"), sep = "\t")

##########################################
#######ii.anno file GARFILED format#######
##########################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
######EAS
EAS=fread("/data1/gy/public/LDSC_ref/1000G.EAS.QC.1_22.bim",header=F)
EAS=EAS[,c(1,4)]
######caSNP
caSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_caQTL_FDR0.1.SNP.vcf")
dim(caSNP)
#[1] 291761      5
caSNP=caSNP[,c(1,2)]
######asSNP
asSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_ASCA_FDR0.1.SNP.vcf")
dim(asSNP)
#[1] 10470     5
asSNP=asSNP[,c(1,2)]
##
for (i in 1:22){
  EAS_chr=subset(EAS,V1==i)
  caSNP_chr=subset(caSNP,caSNP$'#CHROM'==paste0("chr",i))
  asSNP_chr=subset(asSNP,asSNP$'#CHROM'==paste0("chr",i))
  # 在 EAS_chr 中添加 anno 列  
  EAS_chr <- EAS_chr %>%  
  mutate(anno1 = ifelse(V4 %in% caSNP_chr$POS, 1, 0)) %>%
  mutate(anno2 = ifelse(V4 %in% asSNP_chr$POS, 1, 0)) %>%
  mutate(anno3 = ifelse(V4 %in% caSNP_chr$POS | V4 %in% asSNP_chr$POS, 1, 0))
  EAS_chr$anno=paste0(EAS_chr$anno1,EAS_chr$anno2,EAS_chr$anno3)
  fwrite(EAS_chr[,c(2,6)], paste0("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/anno/chr",i), sep = " ", col.names = F, quote = F)
}

#########################################
##########iii.garfield-prep-chr##########
#########################################
cd /data1/gy/software/garfield-v2
#
GWAS_files=$(ls /data1/gy/public/BBJ_2020_GWAS/garfield_format)  

#
GWAS_array=($GWAS_files)

# 
for GWAS in "${GWAS_array[@]:0:42}"; do
   for i in {1..22}; do  
      ./garfield-prep-chr \
      -ptags /data1/gy/public/garfield-data_EAS/tags/r01/chr$i \
      -ctags /data1/gy/public/garfield-data_EAS/tags/r08/chr$i \
      -maftss /data1/gy/public/garfield-data_EAS/maftssd/chr$i \
      -pval /data1/gy/public/BBJ_2020_GWAS/garfield_format/$GWAS/chr$i \
      -ann  /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/anno/chr$i \
      -o    /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-prep-chr_BBJ2020/garfield.prep.$GWAS.out \
      -chr $i  
   done
   echo "Finished $GWAS"
done

########################################
##########iv.garfield-Meff-Padj##########
########################################
cd /data1/gy/software/garfield-v2
#
GWAS_files=$(ls /data1/gy/public/BBJ_2020_GWAS/garfield_format)  

#
GWAS_array=($GWAS_files)

#
for GWAS in "${GWAS_array[@]:0:42}"; do  
   #
   /Public/gaoyun/software/R-4.2.0/bin/Rscript garfield-Meff-Padj.R \
   -i /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-prep-chr_BBJ2020/garfield.prep.$GWAS.out \
   -o /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-Meff-Padj_BBJ2020/garfield.Meff.$GWAS.out
done  

###################################
##########v.garfield-test##########
###################################
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-test_BBJ2020
cd /data1/gy/software/garfield-v2
#
GWAS_files=$(ls /data1/gy/public/BBJ_2020_GWAS/garfield_format)  

#
GWAS_array=($GWAS_files)

#
for GWAS in "${GWAS_array[@]:0:42}"; do  
/Public/gaoyun/software/R-4.2.0/bin/Rscript garfield-test.R \
-i /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-prep-chr_BBJ2020/garfield.prep.$GWAS.out \
-o /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-test_BBJ2020/garfield.test.$GWAS.out \
-l /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/anno/link_file.txt \
-pt 5e-8 \
-b m5,n5,t5
done  

############################################
##########vi.garfield-test results##########
############################################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-test_BBJ2020")
library(data.table)
library(dplyr)
# 
enrichment_files <- list.files(pattern = "garfield.test.bbj.*$") 
#
data_rbind <- data.table()  
#
for (file in enrichment_files) {  
    data <- fread(file) 
    
    #
    trait <- sub(".*(bbj-a-[0-9]+).*", "\\1", basename(file))
    
    #
    data[, trait := trait]  
    
    #
    data_rbind <- rbind(data_rbind, data)  
}
#
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield_enrichment_result_BBJ2020")
BBJ=openxlsx::read.xlsx("/data1/gy/public/BBJ_2020_GWAS/raw/2020BBJ_GWAS.xlsx")
library(dplyr)
data_rbind=left_join(data_rbind,BBJ,by=c("trait"="BBJid"))
fwrite(data_rbind,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield_enrichment_result_BBJ2020/caQTL_ASCA_enrichment_in_BBJ_GWAS.GARFIELD.result",sep="\t",quote=F,col.names=T)
q()

#################################################################
#######2.4GWAS from BBJ published in 2023 (PMID:38036781)#######
#################################################################
#########################################
#######i.BBJ GWAS GARFIELD format#######
#########################################
cd /data1/gy/public/BBJ_2023_GWAS/raw
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
###
#
gwas_files <- list.files(pattern = "GCST.*\\.tsv.gz$") 
#
output_dir <- "/data1/gy/public/BBJ_2023_GWAS/garfield_format/"  
#
results_summary <- data.frame(
  file = character(),
  significant_count = integer(),
  stringsAsFactors = FALSE
)

#
for (i in seq_along(gwas_files)) {

  #
  GWAS <- fread(gwas_files[i], select = c("chromosome","base_pair_location","p_value"))

  #
  sig_count <- sum(GWAS$p_value < 5e-8, na.rm = TRUE)

  #
  results_summary <- rbind(results_summary,
                           data.frame(file = gwas_files[i],
                                      significant_count = sig_count,
                                      stringsAsFactors = FALSE))

  #
  if (sig_count > 0) {

    #
    base_name <- sub("\\.tsv.gz$", "", gwas_files[i]) 

    #
    dir.create(paste0(output_dir, base_name), showWarnings = FALSE, recursive = TRUE)

    #
    for (m in 1:22) {
      GWAS_chr <- subset(GWAS, chromosome == m)
      fwrite(GWAS_chr[, c(2,3)],
             file = paste0(output_dir, base_name, "/chr", m),
             col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
    }

    cat("Processed:", gwas_files[i], " (Found", sig_count, "significant SNPs)\n")

  } else {
    cat("Skipped:", gwas_files[i], " (no p<5e-8 SNPs)\n")
  }
}

#
fwrite(results_summary, file = paste0(output_dir, "GWAS_significant_counts.tsv"), sep = "\t")

##########################################
#######ii.anno file GARFILED format#######
##########################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
######EAS
EAS=fread("/data1/gy/public/LDSC_ref/1000G.EAS.QC.1_22.bim",header=F)
EAS=EAS[,c(1,4)]
######caSNP
caSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_caQTL_FDR0.1.SNP.vcf")
dim(caSNP)
#[1] 291761      5
caSNP=caSNP[,c(1,2)]
######asSNP
asSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_ASCA_FDR0.1.SNP.vcf")
dim(asSNP)
#[1] 10470     5
asSNP=asSNP[,c(1,2)]
##
for (i in 1:22){
  EAS_chr=subset(EAS,V1==i)
  caSNP_chr=subset(caSNP,caSNP$'#CHROM'==paste0("chr",i))
  asSNP_chr=subset(asSNP,asSNP$'#CHROM'==paste0("chr",i))
  #
  EAS_chr <- EAS_chr %>%  
  mutate(anno1 = ifelse(V4 %in% caSNP_chr$POS, 1, 0)) %>%
  mutate(anno2 = ifelse(V4 %in% asSNP_chr$POS, 1, 0)) %>%
  mutate(anno3 = ifelse(V4 %in% caSNP_chr$POS | V4 %in% asSNP_chr$POS, 1, 0))
  EAS_chr$anno=paste0(EAS_chr$anno1,EAS_chr$anno2,EAS_chr$anno3)
  fwrite(EAS_chr[,c(2,6)], paste0("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/anno/chr",i), sep = " ", col.names = F, quote = F)
}

#########################################
##########iii.garfield-prep-chr##########
#########################################
cd /data1/gy/software/garfield-v2
#
GWAS_files=$(ls /data1/gy/public/BBJ_2023_GWAS/garfield_format)  
#
GWAS_array=($GWAS_files)

# 
for GWAS in "${GWAS_array[@]:0:4}"; do
   for i in {1..22}; do  
      ./garfield-prep-chr \
      -ptags /data1/gy/public/garfield-data_EAS/tags/r01/chr$i \
      -ctags /data1/gy/public/garfield-data_EAS/tags/r08/chr$i \
      -maftss /data1/gy/public/garfield-data_EAS/maftssd/chr$i \
      -pval /data1/gy/public/BBJ_2023_GWAS/garfield_format/$GWAS/chr$i \
      -ann  /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/anno/chr$i \
      -o    /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-prep-chr_BBJ2023/garfield.prep.$GWAS.out \
      -chr $i  
   done
   echo "Finished $GWAS"
done

########################################
##########iv.garfield-Meff-Padj##########
########################################
cd /data1/gy/software/garfield-v2
#
GWAS_files=$(ls /data1/gy/public/BBJ_2023_GWAS/garfield_format)  

#
GWAS_array=($GWAS_files)

#
for GWAS in "${GWAS_array[@]:0:4}"; do  
   #
   /Public/gaoyun/software/R-4.2.0/bin/Rscript garfield-Meff-Padj.R \
   -i /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-prep-chr_BBJ2023/garfield.prep.$GWAS.out \
   -o /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-Meff-Padj_BBJ2023/garfield.Meff.$GWAS.out
done  

###################################
##########v.garfield-test##########
###################################
cd /data1/gy/software/garfield-v2
#
GWAS_files=$(ls /data1/gy/public/BBJ_2023_GWAS/garfield_format)  

#
GWAS_array=($GWAS_files)

#
for GWAS in "${GWAS_array[@]:0:4}"; do  
  /Public/gaoyun/software/R-4.2.0/bin/Rscript garfield-test.R \
  -i /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-prep-chr_BBJ2023/garfield.prep.$GWAS.out \
  -o /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-test_BBJ2023/garfield.test.$GWAS.out \
  -l /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/anno/link_file.txt \
  -pt 5e-8 \
  -b m5,n5,t5
done  

############################################
##########vi.garfield-test results##########
############################################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield-test_BBJ2023")
library(data.table)
library(dplyr)
#
enrichment_files <- list.files(pattern = "garfield.test.GCST.*$") 
#
data_rbind <- data.table() 
# 
for (file in enrichment_files) {  
    data <- fread(file)
    
    #
    if (data[21, "NThresh"] != 0) {  
        #
        trait <- sub(".*(GCST[0-9]+).*", "\\1", basename(file))
        
        #
        data[, trait := trait]  
        
        #
        data_rbind <- rbind(data_rbind, data)  
    }  
}  
#
dir.create("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield_enrichment_result_BBJ2023")
BBJ=openxlsx::read.xlsx("/data1/gy/public/BBJ_2023_GWAS/raw/2023BBJ_GWAS.xlsx")
library(dplyr)
data_rbind=left_join(data_rbind,BBJ,by=c("trait"="accessionId"))
fwrite(data_rbind,"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield_enrichment_result_BBJ2023/caQTL_ASCA_enrichment_in_BBJ_GWAS.GARFIELD.result",sep="\t",quote=F,col.names=T)
q()

#######################################################################
##########3.enrichment forestplot for significantly enriched ##########
#######################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATAC_for_review/Figure3D/output")
library(data.table)
library(dplyr)
data1=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield_enrichment_result_BBJ2020/caQTL_ASCA_enrichment_in_BBJ_GWAS.GARFIELD.result")
data1=data1[,-21]
data2=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/garfield_enrichment_result_BBJ2023/caQTL_ASCA_enrichment_in_BBJ_GWAS.GARFIELD.result")
data=rbind(data1,data2)
data=subset(data,Annotation=="caQTL|ASCA")
#
data[NAnnotThesh == 0, `:=`(  
    OR = 1,  
    Beta = 0,  
    SE = 0,  
    CI95_lower = 0,  
    CI95_upper = 0  
)]  


#
data$color <- ifelse(grepl("gastric", data$Disease, ignore.case = TRUE), "Gastric-related", "non-Gastric-related") 
#
data=subset(data,Pvalue < 0.05 & OR > 1)
data=subset(data,!(Disease %in% c("duodenal ulcer, gastric ulcer","peptic ulcer disease"))) ## These are removed because they are duplicate/overlapping phenotypes with gastric and duodenal ulcers.
data[, Disease := factor(Disease, levels = rev(c("Gastric cancer","gastric ulcer","Chronic hepatitis C","duodenal ulcer","Chronic hepatitis B","Hepatocellular carcinoma","Rheumatoid arthritis","Graves' disease",
                                                 "Coronary artery disease","Atopic dermatitis","Asthma","Type 2 diabetes")))]  

##forest plot
library(ggplot2)
plot_forest <- function(sub_data, PThresh) {  
    ggplot(sub_data, aes(x = log2(OR), y = Disease)) +  
        geom_point(aes(color = color), size = 3) +
        geom_errorbarh(aes(xmin = log2(exp(CI95_lower)), xmax = log2(exp(CI95_upper)), color = color), height = 0.2) +
        scale_color_manual(values = c("Gastric-related" = "#E31A1C", "non-Gastric-related" = "#225EA8"))+ 
        labs(title = "",  
             x = "OR",  
             y = "Disease") +  
        scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6), 
                           labels = c("1", "2", "4", "8", "16", "32", "64")) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
        theme_minimal() +  
        theme(legend.position = "none", 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"), 
              axis.ticks.length = unit(0.25, "cm"), 
              axis.ticks = element_line(colour = "black")) +  
        theme(axis.ticks.y = element_line(size = 0.3), 
              axis.ticks.x = element_line(size = 0.3))
}  

##
PThresh_sub=5e-08
sub_data =data[PThresh==PThresh_sub]
p <- plot_forest(sub_data, PThresh_sub)
ggsave(p,filename=paste0("caQTL_ASCA_enrichment_in_BBJ_GWAS.PThresh_",PThresh_sub,".garfield.forestplot.pdf"),height=7,width=6) ##Figure3D
