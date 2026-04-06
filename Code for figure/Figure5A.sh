######################
#######Figure5A#######
######################
##########################################
#######1.caSNP/asSNP GC GWAS signal#######
##########################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(qqman)
################
######caSNP
caSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_caQTL_FDR0.1.SNP.vcf")
dim(caSNP)
#[1] 291761      5
######asSNP
asSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/VEP_geno/input/95sample_ASCA_FDR0.1.SNP.vcf")
dim(asSNP)
#[1] 10470     5
######GWAS
GWAS=fread("/data1/gy/EAS_GWAS_meta/GWAS_QCandFilter.inEAS.new.metaResult")
GWAS$chr=paste0("chr",GWAS$chr)
GWAS$ID=paste(GWAS$chr,GWAS$bp,GWAS$Allele1,GWAS$Allele2,sep=":")
GWAS=GWAS[,c('P-value','ID')]
#################
######caSNP
caSNP1=caSNP
caSNP2=caSNP
caSNP1$ID=paste(caSNP1$'#CHROM',caSNP1$POS,caSNP1$REF,caSNP1$ALT,sep=":")
caSNP2$ID=paste(caSNP2$'#CHROM',caSNP2$POS,caSNP2$ALT,caSNP2$REF,sep=":")
caSNP=rbind(caSNP1,caSNP2)
caSNP=merge(caSNP,GWAS)
caSNP = caSNP[order(caSNP$'P-value',decreasing = F),]
caSNP$obs=-log10(caSNP$'P-value')
caSNP$exp=-log10(ppoints(length(caSNP$'P-value')))
caSNP=caSNP[,c('obs','exp')]
caSNP$group="caQTL"
######asSNP
asSNP1=asSNP
asSNP2=asSNP
asSNP1$ID=paste(asSNP1$'#CHROM',asSNP1$POS,asSNP1$REF,asSNP1$ALT,sep=":")
asSNP2$ID=paste(asSNP2$'#CHROM',asSNP2$POS,asSNP2$ALT,asSNP2$REF,sep=":")
asSNP=rbind(asSNP1,asSNP2)
asSNP=merge(asSNP,GWAS)
asSNP = asSNP[order(asSNP$'P-value',decreasing = F),]
asSNP$obs=-log10(asSNP$'P-value')
asSNP$exp=-log10(ppoints(length(asSNP$'P-value')))
asSNP=asSNP[,c('obs','exp')]
asSNP$group="ASCA"
#####GWAS
GWAS = GWAS[order(GWAS$'P-value',decreasing = F),]
GWAS$obs = -log10(GWAS$'P-value')
GWAS$exp = -log10(ppoints(length(GWAS$'P-value')))
GWAS=GWAS[,c('obs','exp')]
GWAS$group="Genome-wide"
#####
data=rbind(caSNP,asSNP,GWAS)
data$group=factor(data$group,level=c("caQTL","ASCA","Genome-wide"))

######################
#######2.QQplot#######
######################
library(ggplot2)
setwd("/data1/gy/ATAC_for_review/Figure5A/output")
#
p <- ggplot(data, aes(x = exp, y = obs)) +  
  geom_point(aes(color = group), alpha = 0.8, size = 3, shape = 16) +
  geom_abline() +  
  scale_x_continuous(expand = c(0.01, 0.01)) +  
  labs(x = expression(paste("Expected -log"[10], "(", italic('P'), ")")),  
       y = expression(paste("Observed -log"[10], "(", italic('P'), ")"))) +  
  theme_bw() +  
  theme(panel.grid = element_blank(),  
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12, color = 'black', face = 'plain'),  
        axis.title = element_text(size = 12, color = 'black', face = 'plain'),  
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "none") + 
        scale_color_manual(values = c("caQTL" = "#EF3B2C", "ASCA" = "#EC7014", "Genome-wide" = "black"),   
                     guide = "legend")
ggsave(p,file="caQTL_ASCA_Genome-wide.GWAS_signal.QQplot.tiff",dpi=600) ##Figure5A
