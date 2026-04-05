#####################################
#############FigureS8G-H#############
#####################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
###################################################
#############1.caQTL eQTL GWAS SNPlist#############
###################################################
#######caqtl
peakID="chr14:34896054-34896554"
caqtl=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")
caqtl=subset(caqtl,pheno_id %in% peakID)
caqtl=caqtl[,c("pheno_id","var_id","p_nominal")]
caqtl$chrbp <- sub("([^:]+:)([0-9]+).*", "\\1\\2", caqtl$var_id)  
caqtl$REF <- sub(".*:([A-Z])(:[A-Z])$", "\\1", caqtl$var_id)  
caqtl$ALT <- sub(".*:([A-Z]):([A-Z])$", "\\2", caqtl$var_id)  
colnames(caqtl)[2:3]=c("caqtlID","caqtlP")

#######262eqtl
symbol="SPTSSA"
eqtl=fread("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.with_symbol.all.chr.txt")
eqtl=subset(eqtl,gene_name %in% symbol)
eqtl=eqtl[,c("phenotype_id","gene_name","variant_id","pval_nominal","REF","ALT")]
eqtl$chrbp <- sub("([^:]+:)([0-9]+).*", "\\1\\2", eqtl$variant_id)  
eqtl$chrbp <- paste0("chr",eqtl$chrbp)
eqtl=eqtl[,c(1:4,7,5:6)]
colnames(eqtl)[3:4]=c("eqtlID","eqtlP")
colnames(eqtl)[6:7]=c("ref","alt")

#######sceqtl
symbol="SPTSSA"
sceqtl=fread("/data1/gy/sceQTL/eQTL_CHN100K/epi/epi/epi_sceQTL_result.allpairs.withsymbol.withFDR.txt")
sceqtl=subset(sceqtl,gene_name %in% symbol)
sceqtl=sceqtl[,c("gene_id","gene_name","variant_id","pval_nominal")]
sceqtl$chrbp <- sub("([^:]+:)([0-9]+).*", "\\1\\2", sceqtl$variant_id)  
sceqtl$ref <- sub(".*:([A-Z])(:[A-Z])$", "\\1", sceqtl$variant_id)  
sceqtl$alt <- sub(".*:([A-Z]):([A-Z])$", "\\2", sceqtl$variant_id)  
colnames(sceqtl)[3:4]=c("eqtlID","eqtlP")

############################################
#############2.file for plot#############
############################################
##LDr2 file
LDr2=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/chr14:34896054-34896554/LDr2/chr14:34896054-34896554.ld")
LDr2=LDr2[,c(3,6,7)]

#########################################
##########i.caqtl_262eqtl_coloc##########
#########################################
caqtl_eqtl_coloc=merge(caqtl,eqtl,by="chrbp")
dim(caqtl_eqtl_coloc)
#[1] 798  12
caqtl_eqtl_coloc1 = caqtl_eqtl_coloc[which((caqtl_eqtl_coloc$REF == caqtl_eqtl_coloc$ref & caqtl_eqtl_coloc$ALT == caqtl_eqtl_coloc$alt) 
                                     | (caqtl_eqtl_coloc$REF == caqtl_eqtl_coloc$alt & caqtl_eqtl_coloc$ALT == caqtl_eqtl_coloc$ref)),]
dim(caqtl_eqtl_coloc1)
#[1] 762  12
caqtl_eqtl_coloc = caqtl_eqtl_coloc1
dim(caqtl_eqtl_coloc)
#[1] 762  12
##add LDr2
caqtl_eqtl_coloc=merge(caqtl_eqtl_coloc,LDr2,by.x="caqtlID",by.y="SNP_B")
dim(caqtl_eqtl_coloc)
#[1] 762  14
fwrite(caqtl_eqtl_coloc,"/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/chr14:34896054-34896554/colocalization_plot/caqtl_eqtl_coloc.SPTSSA.txt",col.names=T,row.names=F,quote=F,sep="\t")

#########################################
##########ii.caqtl_sceqtl_coloc##########
#########################################
caqtl_sceqtl_coloc=merge(caqtl,sceqtl,by="chrbp")
dim(caqtl_sceqtl_coloc)
#[1] 823  12
caqtl_sceqtl_coloc1 = caqtl_sceqtl_coloc[which((caqtl_sceqtl_coloc$REF == caqtl_sceqtl_coloc$ref & caqtl_sceqtl_coloc$ALT == caqtl_sceqtl_coloc$alt) 
                                     | (caqtl_sceqtl_coloc$REF == caqtl_sceqtl_coloc$alt & caqtl_sceqtl_coloc$ALT == caqtl_sceqtl_coloc$ref)),]
dim(caqtl_sceqtl_coloc1)
#[1] 823  12
caqtl_sceqtl_coloc = caqtl_sceqtl_coloc1
dim(caqtl_sceqtl_coloc)
#[1] 823  12
##add LDr2
caqtl_sceqtl_coloc=merge(caqtl_sceqtl_coloc,LDr2,by.x="caqtlID",by.y="SNP_B")
dim(caqtl_sceqtl_coloc)
#[1] 823  14
fwrite(caqtl_sceqtl_coloc,"/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/chr14:34896054-34896554/colocalization_plot/caqtl_sceqtl_coloc.SPTSSA.txt",col.names=T,row.names=F,quote=F,sep="\t")

################################
#############3.plot#############
################################
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
######################################
##########i.caqtl_eqtl_coloc##########
######################################
caqtl_eqtl_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/chr14:34896054-34896554/colocalization_plot/caqtl_eqtl_coloc.SPTSSA.txt")
###
ref1 <- "chr14:34895977"  
shape = ifelse(caqtl_eqtl_coloc$chrbp == ref1, 23, 21)
names(shape) = caqtl_eqtl_coloc$chrbp
###
size = ifelse(caqtl_eqtl_coloc$chrbp == ref1, 3, 2)
names(size) = caqtl_eqtl_coloc$chrbp
###
caqtl_eqtl_coloc$label = ifelse(caqtl_eqtl_coloc$chrbp == ref1, caqtl_eqtl_coloc$caqtlID, '')
###
color = as.character(cut(caqtl_eqtl_coloc$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"), include.lowest=TRUE))
names(color) = caqtl_eqtl_coloc$chrbp
###
title = "COLOC.PP4 = 0.19"
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

ref_eqtlP <- -log10(caqtl_eqtl_coloc[caqtl_eqtl_coloc$chrbp == ref1, ]$eqtlP)  
ref_caqtlP <- -log10(caqtl_eqtl_coloc[caqtl_eqtl_coloc$chrbp == ref1, ]$caqtlP)  
#
p_caqtl_eqtl_coloc <- p_caqtl_eqtl_coloc +   
  geom_point(aes(x = ref_eqtlP, y = ref_caqtlP),   
             shape = 23, size = 3, fill = "purple")  
  
###add LD r2 legend
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
ggsave(p_caqtl_eqtl_coloc_withlegend,
       filename="caqtl_eqtl_coloc.SPTSSA.scatterplot.pdf",
       width=4,height=4)  ##FigureS8H

#########################################
##########ii.caqtl_sceqtl_coloc##########
#########################################
caqtl_sceqtl_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/chr14:34896054-34896554/colocalization_plot/caqtl_sceqtl_coloc.SPTSSA.txt")
###
ref1 <- "chr14:34895977"  
shape = ifelse(caqtl_sceqtl_coloc$chrbp == ref1, 23, 21)
names(shape) = caqtl_sceqtl_coloc$chrbp
###
size = ifelse(caqtl_sceqtl_coloc$chrbp == ref1, 3, 2)
names(size) = caqtl_sceqtl_coloc$chrbp
###
caqtl_sceqtl_coloc$label = ifelse(caqtl_sceqtl_coloc$chrbp == ref1, caqtl_sceqtl_coloc$caqtlID, '')
###
color = as.character(cut(caqtl_sceqtl_coloc$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"), include.lowest=TRUE))
names(color) = caqtl_sceqtl_coloc$chrbp
###
title = "COLOC.PP4 = 0.96"
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
  
###add LD r2 legend
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
       filename="caqtl_sceqtl_coloc.SPTSSA.scatterplot.pdf",
       width=4,height=4) ##FigureS8G
