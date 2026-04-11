######################################################
###########FigureS10A&FigureS10D&FigureS10H###########
######################################################
mkdir -p /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
##########################################################
#############1.create caQTL eQTL GWAS SNPlist#############
##########################################################
#######caqtl
peakID="chr4:48072992-48073492"
caqtl=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
caqtl=subset(caqtl,pheno_id %in% peakID)
caqtl=caqtl[,c("pheno_id","var_id","p_nominal")]
caqtl$chrbp <- sub("([^:]+:)([0-9]+).*", "\\1\\2", caqtl$var_id)  
caqtl$REF <- sub(".*:([A-Z])(:[A-Z])$", "\\1", caqtl$var_id)  
caqtl$ALT <- sub(".*:([A-Z]):([A-Z])$", "\\2", caqtl$var_id)  
colnames(caqtl)[2:3]=c("caqtlID","caqtlP")

#######262eqtl
symbol="NFXL1"
eqtl=fread("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.with_symbol.all.chr.txt")
eqtl=subset(eqtl,gene_name %in% symbol)
eqtl=eqtl[,c("phenotype_id","gene_name","variant_id","pval_nominal","REF","ALT")]
eqtl$chrbp <- sub("([^:]+:)([0-9]+).*", "\\1\\2", eqtl$variant_id)  
eqtl$chrbp <- paste0("chr",eqtl$chrbp)
eqtl=eqtl[,c(1:4,7,5:6)]
colnames(eqtl)[3:4]=c("eqtlID","eqtlP")
colnames(eqtl)[6:7]=c("ref","alt")

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
######################################
##########i.caqtl_GWAS_coloc##########
######################################
caqtl_GWAS_coloc=merge(caqtl,GWAS,by="chrbp")
dim(caqtl_GWAS_coloc)
#[1] 879  10
caqtl_GWAS_coloc1 = caqtl_GWAS_coloc[which((caqtl_GWAS_coloc$Allele1 == caqtl_GWAS_coloc$REF & caqtl_GWAS_coloc$Allele2 == caqtl_GWAS_coloc$ALT) 
                                     | (caqtl_GWAS_coloc$Allele1 == caqtl_GWAS_coloc$ALT & caqtl_GWAS_coloc$Allele2 == caqtl_GWAS_coloc$REF)),]
dim(caqtl_GWAS_coloc1)
#[1] 878  10
#######
caqtl_GWAS_coloc2 = caqtl_GWAS_coloc[!which((caqtl_GWAS_coloc$Allele1 == caqtl_GWAS_coloc$REF & caqtl_GWAS_coloc$Allele2 == caqtl_GWAS_coloc$ALT) 
                                     | (caqtl_GWAS_coloc$Allele1 == caqtl_GWAS_coloc$ALT & caqtl_GWAS_coloc$Allele2 == caqtl_GWAS_coloc$REF)),]
dim(caqtl_GWAS_coloc2)
#[1]  1 10
caqtl_GWAS_coloc = caqtl_GWAS_coloc1
dim(caqtl_GWAS_coloc)
#[1] 878  10
##
caqtl_GWAS_coloc=merge(caqtl_GWAS_coloc,LDr2,by.x="caqtlID",by.y="SNP_B")
dim(caqtl_GWAS_coloc)
#[1] 878  12
fwrite(caqtl_GWAS_coloc,"/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1/caqtl_GWAS_coloc.txt",col.names=T,row.names=F,quote=F,sep="\t")

#########################################
##########ii.caqtl_262eqtl_coloc##########
#########################################
caqtl_eqtl_coloc=merge(caqtl,eqtl,by="chrbp")
dim(caqtl_eqtl_coloc)
#[1] 845  12
caqtl_eqtl_coloc1 = caqtl_eqtl_coloc[which((caqtl_eqtl_coloc$REF == caqtl_eqtl_coloc$ref & caqtl_eqtl_coloc$ALT == caqtl_eqtl_coloc$alt) 
                                     | (caqtl_eqtl_coloc$REF == caqtl_eqtl_coloc$alt & caqtl_eqtl_coloc$ALT == caqtl_eqtl_coloc$ref)),]
dim(caqtl_eqtl_coloc1)
#[1] 833  12
caqtl_eqtl_coloc = caqtl_eqtl_coloc1
dim(caqtl_eqtl_coloc)
#[1] 833  12
##
caqtl_eqtl_coloc=merge(caqtl_eqtl_coloc,LDr2,by.x="caqtlID",by.y="SNP_B")
dim(caqtl_eqtl_coloc)
#[1] 828  14
fwrite(caqtl_eqtl_coloc,"/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1/caqtl_eqtl_coloc.txt",col.names=T,row.names=F,quote=F,sep="\t")

#################################################
##########iii.caqtl_eqtl_GWAS_hyprcoloc##########
#################################################
caqtl_eqtl_coloc=merge(caqtl,eqtl,by="chrbp")
hyprcoloc=merge(caqtl_eqtl_coloc,GWAS,by="chrbp")
dim(hyprcoloc)
#[1] 841  16
hyprcoloc1 = hyprcoloc[which((hyprcoloc$Allele1 == hyprcoloc$REF & hyprcoloc$Allele2 == hyprcoloc$ALT) 
                          | (hyprcoloc$Allele1 == hyprcoloc$ALT & hyprcoloc$Allele2 == hyprcoloc$REF)),]
dim(hyprcoloc1)
#[1] 840  16
#######筛选caQTL和GWAS不相同的SNP
hyprcoloc2 = hyprcoloc[!which((hyprcoloc$Allele1 == hyprcoloc$REF & hyprcoloc$Allele2 == hyprcoloc$ALT) 
                          | (hyprcoloc$Allele1 == hyprcoloc$ALT & hyprcoloc$Allele2 == hyprcoloc$REF)),]
dim(hyprcoloc2)
#[1]  1 16
hyprcoloc = hyprcoloc1
dim(hyprcoloc)
#[1] 840  16
##添加LD
hyprcoloc=merge(hyprcoloc,LDr2,by.x="caqtlID",by.y="SNP_B")
dim(hyprcoloc)
#[1] 840  18
fwrite(hyprcoloc,"/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1/caqtl_eqtl_GWAS_hyprcoloc.txt",col.names=T,row.names=F,quote=F,sep="\t")

################################
#############3.plot#############
################################
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
######################################
##########i.caqtl_GWAS_coloc##########
######################################
caqtl_GWAS_coloc=fread("/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1/caqtl_GWAS_coloc.txt")
###
ref1 <- "chr4:48090040"  
ref2 <- "chr4:48073300"  
shape = ifelse(caqtl_GWAS_coloc$chrbp == ref1 |caqtl_GWAS_coloc$chrbp == ref2, 23, 21)
names(shape) = caqtl_GWAS_coloc$chrbp
###
size = ifelse(caqtl_GWAS_coloc$chrbp == ref1 |caqtl_GWAS_coloc$chrbp == ref2, 3, 2)
names(size) = caqtl_GWAS_coloc$chrbp
###
caqtl_GWAS_coloc$label = ifelse(caqtl_GWAS_coloc$chrbp == ref1 |caqtl_GWAS_coloc$chrbp == ref2, caqtl_GWAS_coloc$caqtlID, '')
###
color = as.character(cut(caqtl_GWAS_coloc$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"), include.lowest=TRUE))
names(color) = caqtl_GWAS_coloc$chrbp
###
title = "COLOC.PP4 = 0.80"
p_caqtl_GWAS_coloc = ggplot(caqtl_GWAS_coloc,aes(x=-log10(GWASP),y=-log10(caqtlP)))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=caqtl_GWAS_coloc[caqtl_GWAS_coloc$label!='',],aes(x=-log10(GWASP),y=-log10(caqtlP),fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape,guide='none')+
          scale_size_manual(values=size,guide='none')+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab("GWAS -log10(P)")+
          ylab("caQTL -log10(P)")+
          labs(title = title)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1)))
          )
#
ref_GWASP <- -log10(caqtl_GWAS_coloc[caqtl_GWAS_coloc$chrbp == ref1, ]$GWASP)  
ref_caqtlP <- -log10(caqtl_GWAS_coloc[caqtl_GWAS_coloc$chrbp == ref1, ]$caqtlP)  
#
p_caqtl_GWAS_coloc <- p_caqtl_GWAS_coloc +   
  geom_point(aes(x = ref_GWASP, y = ref_caqtlP),   
             shape = 23, size = 3, fill = "purple")  
  
###
legend_box = data.frame(x = 0.25, y = seq(0.85, 0.69, -0.04))
p_caqtl_GWAS_coloc_withlegend = ggdraw(p_caqtl_GWAS_coloc)+
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

##
setwd("/data1/gy/ATAC_for_review/FigureS10A&FigureS10D&FigureS10H/output")
ggsave(p_caqtl_GWAS_coloc_withlegend,
       filename="chr4:48072992-48073492~GWAS.caqtl_GWAS_coloc.scatterplot.pdf",
       width=4,height=4)  ##FigureS10A

######################################
##########ii.caqtl_eqtl_coloc##########
######################################
caqtl_eqtl_coloc=fread("/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1/caqtl_eqtl_coloc.txt")
###
ref1 <- "chr4:48090040"  
ref2 <- "chr4:48073300"  
shape = ifelse(caqtl_eqtl_coloc$chrbp == ref1 |caqtl_eqtl_coloc$chrbp == ref2, 23, 21)
names(shape) = caqtl_eqtl_coloc$chrbp
###
size = ifelse(caqtl_eqtl_coloc$chrbp == ref1 |caqtl_eqtl_coloc$chrbp == ref2, 3, 2)
names(size) = caqtl_eqtl_coloc$chrbp
###
caqtl_eqtl_coloc$label = ifelse(caqtl_eqtl_coloc$chrbp == ref1 |caqtl_eqtl_coloc$chrbp == ref2, caqtl_eqtl_coloc$caqtlID, '')
###
color = as.character(cut(caqtl_eqtl_coloc$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"), include.lowest=TRUE))
names(color) = caqtl_eqtl_coloc$chrbp
###
title = "COLOC.PP4 = 0.72"
p_caqtl_eqtl_coloc = ggplot(caqtl_eqtl_coloc,aes(x=-log10(eqtlP),y=-log10(caqtlP)))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=caqtl_eqtl_coloc[caqtl_eqtl_coloc$label!='',],aes(x=-log10(eqtlP),y=-log10(caqtlP),fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape,guide='none')+
          scale_size_manual(values=size,guide='none')+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab("Bulk eQTL -log10(P)")+
          ylab("caQTL -log10(P)")+
          labs(title = title)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(title, "\n") + 1)))
          )
#
ref_eqtlP <- -log10(caqtl_eqtl_coloc[caqtl_eqtl_coloc$chrbp == ref1, ]$eqtlP)  
ref_caqtlP <- -log10(caqtl_eqtl_coloc[caqtl_eqtl_coloc$chrbp == ref1, ]$caqtlP)  
#
p_caqtl_eqtl_coloc <- p_caqtl_eqtl_coloc +   
  geom_point(aes(x = ref_eqtlP, y = ref_caqtlP),   
             shape = 23, size = 3, fill = "purple")  
  
###
legend_box = data.frame(x = 0.25, y = seq(0.85, 0.69, -0.04))
p_caqtl_eqtl_coloc_withlegend = ggdraw(p_caqtl_eqtl_coloc)+
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

##
setwd("/data1/gy/ATAC_for_review/FigureS10A&FigureS10D&FigureS10H/output")
ggsave(p_caqtl_eqtl_coloc_withlegend,
       filename="chr4:48072992-48073492~NFXL1.caqtl_eqtl_coloc.scatterplot.pdf",
       width=4,height=4) ##FigureS10D

#################################
##########iii.hyprcoloc##########(3D scatterplot)
#################################
###
system("mkdir -p /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1/LDr2")
##LDr2 ref:hyprcoloc top SNP
system("plink --bfile /data1/gy/ATACseq_RWAS/RWAS_STITCH/ldref/1KG_EAS_no_indel_shared --r2 --ld-snp chr4:48076486:G:A --ld-window-kb 400 --ld-window 99999 --ld-window-r2 0 --out /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1/LDr2/chr4:48072992-48073492")
#install.packages("scatterplot3d")

library("scatterplot3d")
hyprcoloc=fread("/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/colocalization_plot/NFXL1/caqtl_eqtl_GWAS_hyprcoloc.txt")
###
ref1 <- "chr4:48076486"  ##hyprcoloc top SNP
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
setwd("/data1/gy/ATAC_for_review/FigureS10A&FigureS10D&FigureS10H/output")
title = "Hyprcoloc.PP = 0.47"
pdf("chr4:48072992-48073492~eQTL~GWAS.hyprcoloc.scatterplot.pdf",
    width=7,height=7) ##FigureS10H
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
              xlab = "Bulk eQTL -log10(P)",
              ylab = "GWAS -log10(P)",
              zlab = "caQTL -log10(P)")
dev.off()
q()
