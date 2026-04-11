###########################################
###########FigureS10E&FigureS10G###########
###########################################
mkdir -p /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NIPAL1
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
############################################################
#############1.create caQTL sceQTL GWAS SNPlist#############
############################################################
#######caqtl
peakID="chr4:48072992-48073492"
caqtl=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
caqtl=subset(caqtl,pheno_id %in% peakID)
caqtl=caqtl[,c("pheno_id","var_id","p_nominal")]
caqtl$chrbp <- sub("([^:]+:)([0-9]+).*", "\\1\\2", caqtl$var_id)  
caqtl$REF <- sub(".*:([A-Z])(:[A-Z])$", "\\1", caqtl$var_id)  
caqtl$ALT <- sub(".*:([A-Z]):([A-Z])$", "\\2", caqtl$var_id)  
colnames(caqtl)[2:3]=c("caqtlID","caqtlP")

#######203sceqtl
symbol="NIPAL1"
sceqtl=fread("/data1/gy/sceQTL/eQTL_CHN100K/epi/epi/epi_sceQTL_result.allpairs.withsymbol.withFDR.txt")
sceqtl=subset(sceqtl,gene_name %in% symbol)
sceqtl=sceqtl[,c("gene_id","gene_name","variant_id","pval_nominal")]
sceqtl$chrbp <- sub("([^:]+:)([0-9]+).*", "\\1\\2", sceqtl$variant_id)  
sceqtl$ref <- sub(".*:([A-Z])(:[A-Z])$", "\\1", sceqtl$variant_id)  
sceqtl$alt <- sub(".*:([A-Z]):([A-Z])$", "\\2", sceqtl$variant_id)  
colnames(sceqtl)[3:4]=c("eqtlID","eqtlP")

#######GWAS
GWAS=fread("/data1/gy/EAS_GWAS_meta/GWAS_QCandFilter.new.metaResult")
GWAS$chrbp=paste0("chr",GWAS$chr,":",GWAS$bp) 
GWAS=GWAS[,c("RSID","chrbp","Allele1","Allele2","P-value")]
colnames(GWAS)[5]=c("GWASP")

##########################################
#############2.files for plot#############
##########################################
##LD
LDr2=fread("/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/LDr2/chr4:48072992-48073492.ld")
LDr2=LDr2[,c(3,6,7)]

########################################
##########i.caqtl_sceqtl_coloc##########
########################################
caqtl_sceqtl_coloc=merge(caqtl,sceqtl,by="chrbp")
dim(caqtl_sceqtl_coloc)
#[1] 845  12
caqtl_sceqtl_coloc1 = caqtl_sceqtl_coloc[which((caqtl_sceqtl_coloc$REF == caqtl_sceqtl_coloc$ref & caqtl_sceqtl_coloc$ALT == caqtl_sceqtl_coloc$alt) 
                                     | (caqtl_sceqtl_coloc$REF == caqtl_sceqtl_coloc$alt & caqtl_sceqtl_coloc$ALT == caqtl_sceqtl_coloc$ref)),]
dim(caqtl_sceqtl_coloc1)
#[1] 845  12
caqtl_sceqtl_coloc = caqtl_sceqtl_coloc1
dim(caqtl_sceqtl_coloc)
#[1] 845  12
##
caqtl_sceqtl_coloc=merge(caqtl_sceqtl_coloc,LDr2,by.x="caqtlID",by.y="SNP_B")
dim(caqtl_sceqtl_coloc)
#[1] 842  14
fwrite(caqtl_sceqtl_coloc,"/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NIPAL1/caqtl_sceqtl_coloc.txt",col.names=T,row.names=F,quote=F,sep="\t")

##################################################
##########ii.caqtl_sceqtl_GWAS_hyprcoloc##########
##################################################
caqtl_sceqtl_coloc=merge(caqtl,sceqtl,by="chrbp")
hyprcoloc=merge(caqtl_sceqtl_coloc,GWAS,by="chrbp")
dim(hyprcoloc)
#[1] 843  16
hyprcoloc1 = hyprcoloc[which((hyprcoloc$Allele1 == hyprcoloc$REF & hyprcoloc$Allele2 == hyprcoloc$ALT) 
                          | (hyprcoloc$Allele1 == hyprcoloc$ALT & hyprcoloc$Allele2 == hyprcoloc$REF)),]
dim(hyprcoloc1)
#[1] 842  16
#######
hyprcoloc2 = hyprcoloc[!which((hyprcoloc$Allele1 == hyprcoloc$REF & hyprcoloc$Allele2 == hyprcoloc$ALT) 
                          | (hyprcoloc$Allele1 == hyprcoloc$ALT & hyprcoloc$Allele2 == hyprcoloc$REF)),]
dim(hyprcoloc2)
#[1]  1 16
hyprcoloc = hyprcoloc1
dim(hyprcoloc)
#[1] 842  16
##
hyprcoloc=merge(hyprcoloc,LDr2,by.x="caqtlID",by.y="SNP_B")
dim(hyprcoloc)
#[1] 842  18
fwrite(hyprcoloc,"/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NIPAL1/caqtl_sceqtl_GWAS_hyprcoloc.txt",col.names=T,row.names=F,quote=F,sep="\t")

################################
#############3.plot#############
################################
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
######################################
##########i.caqtl_sceqtl_coloc########
######################################
caqtl_sceqtl_coloc=fread("/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NIPAL1/caqtl_sceqtl_coloc.txt")
###
ref1 <- "chr4:48090040"  
ref2 <- "chr4:48073300"  
shape = ifelse(caqtl_sceqtl_coloc$chrbp == ref1 |caqtl_sceqtl_coloc$chrbp == ref2, 23, 21)
names(shape) = caqtl_sceqtl_coloc$chrbp
###
size = ifelse(caqtl_sceqtl_coloc$chrbp == ref1 |caqtl_sceqtl_coloc$chrbp == ref2, 3, 2)
names(size) = caqtl_sceqtl_coloc$chrbp
###
caqtl_sceqtl_coloc$label = ifelse(caqtl_sceqtl_coloc$chrbp == ref1 |caqtl_sceqtl_coloc$chrbp == ref2, caqtl_sceqtl_coloc$caqtlID, '')
###
color = as.character(cut(caqtl_sceqtl_coloc$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"), include.lowest=TRUE))
names(color) = caqtl_sceqtl_coloc$chrbp
###
title = "COLOC.PP4 = 0.93"
p_caqtl_sceqtl_coloc = ggplot(caqtl_sceqtl_coloc,aes(x=-log10(eqtlP),y=-log10(caqtlP)))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=caqtl_sceqtl_coloc[caqtl_sceqtl_coloc$label!='',],aes(x=-log10(eqtlP),y=-log10(caqtlP),fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape,guide='none')+
          scale_size_manual(values=size,guide='none')+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab("Epithelium sc-eQTL -log10(P)")+
          ylab("caQTL -log10(P)")+
          labs(title = title)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1)))
          )
#
ref_sceqtlP <- -log10(caqtl_sceqtl_coloc[caqtl_sceqtl_coloc$chrbp == ref1, ]$eqtlP)  
ref_caqtlP <- -log10(caqtl_sceqtl_coloc[caqtl_sceqtl_coloc$chrbp == ref1, ]$caqtlP)  
#
p_caqtl_sceqtl_coloc <- p_caqtl_sceqtl_coloc +   
  geom_point(aes(x = ref_sceqtlP, y = ref_caqtlP),   
             shape = 23, size = 3, fill = "purple")  
  
###
setwd("/data1/gy/ATAC_for_review/FigureS10E&FigureS10G/output")
legend_box = data.frame(x = 0.25, y = seq(0.85, 0.69, -0.04))
p_caqtl_sceqtl_coloc_withlegend = ggdraw(p_caqtl_sceqtl_coloc)+
              geom_rect(data = legend_box,
              aes(xmin = x, xmax = x + 0.04, ymin = y, ymax = y + 0.04),
              color = "black",
              fill = rev(c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"))) +
              draw_label("1", x = legend_box$x[1] + 0.04, y = legend_box$y[1]+0.04, hjust = -0.7, size = 7) +
              draw_label("0.8", x = legend_box$x[1] + 0.04, y = legend_box$y[1], hjust = -0.3, size = 7) +
              draw_label("0.6", x = legend_box$x[2] + 0.04, y = legend_box$y[2], hjust = -0.3, size = 7) +
              draw_label("0.4", x = legend_box$x[3] + 0.04, y = legend_box$y[3], hjust = -0.3, size = 7) +
              draw_label("0.2", x = legend_box$x[4] + 0.04, y = legend_box$y[4], hjust = -0.3, size = 7) +
              draw_label("0", x = legend_box$x[5] + 0.04, y = legend_box$y[5], hjust = -0.7, size = 7) +
              draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.02, y = legend_box$y[1], vjust = -2, size = 7) 
ggsave(p_caqtl_sceqtl_coloc_withlegend,
       filename="chr4:48072992-48073492~epi_NIPAL1.caqtl_sceqtl_coloc.scatterplot.pdf",
       width=4,height=4)  ##FigureS10E

################################
##########ii.hyprcoloc##########(3D scatterplot)
################################
library("scatterplot3d")
hyprcoloc=fread("/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NIPAL1/caqtl_sceqtl_GWAS_hyprcoloc.txt")
###
ref1 <- "chr4:48090040"  ##hyprcoloc top SNP = caQTL-GWAS coloc top SNP
ref2 <- "chr4:48073300"  
shape = ifelse(hyprcoloc$chrbp == ref1 | hyprcoloc$chrbp == ref2 , 18, 20)
names(shape) = hyprcoloc$chrbp
###
size = ifelse(hyprcoloc$chrbp == ref1 | hyprcoloc$chrbp == ref2 , 3, 2)
names(size) = hyprcoloc$chrbp
###
hyprcoloc$label = ifelse(hyprcoloc$chrbp == ref1 | hyprcoloc$chrbp == ref2 , hyprcoloc$caqtlID, '')
###
color = as.character(cut(hyprcoloc$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"), include.lowest=TRUE))
names(color) = hyprcoloc$chrbp
color[names(color) == ref1] <- "#6F4597"  
###
setwd("/data1/gy/ATAC_for_review/FigureS10E&FigureS10G/output")
title = "Hyprcoloc.PP = 0.69"
pdf("chr4:48072992-48073492~sceQTL~GWAS.hyprcoloc.scatterplot.pdf",
    width=7,height=7)  ##FigureS10G
p_hyprcoloc = scatterplot3d(
              x = -log10(hyprcoloc$eqtlP),
              y = -log10(hyprcoloc$GWASP),
              z = -log10(hyprcoloc$caqtlP),
              pch = shape,
              angle = 40, 
              color = color,
              cex.lab = 1,cex.axis = 1,cex.symbols = 2,
              main = title,
              type = 'p',
              xlab = "Epithelium sc-eQTL -log10(P)",
              ylab = "GWAS -log10(P)",
              zlab = "caQTL -log10(P)")
dev.off()
q()
