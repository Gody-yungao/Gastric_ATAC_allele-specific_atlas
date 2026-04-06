####################################
########Figure6F&FigureS11G#########
####################################
###########################################################
######1.cell-type-specific expression 0f GATA family#######
###########################################################
conda activate scATACseq_Signac1.9.0_final
R
###########
library(tidyverse)
library(Seurat)
library(RColorBrewer) 
library(patchwork)
library(ggplot2)
library(dplyr)
library(data.table)
##
load('/data1/GC/scRNAseq/normal_203/MIXALL_final.Rdata')
GATA_TF=c("GATA1","GATA2","GATA3","GATA4","GATA5","GATA6","TRPS1")

################################################plot
setwd("/data1/gy/ATAC_for_review/Figure6F&FigureS11G/output")
#
Idents(MIXALL_final) <- MIXALL_final@meta.data$celltype  
#
custom_order <- c("Epithelium", "Tcell", "Bcell", "myeloid", "Mast", "Fibroblasts", "Endotheliocyte")  
#
Idents(MIXALL_final) <- factor(Idents(MIXALL_final), levels = rev(custom_order))  

#
p1=DotPlot(object = MIXALL_final, features=GATA_TF, scale = T)& 
  theme_bw()&
  geom_point(shape=21, aes(size=pct.exp),stroke=1)& 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color = 'black',angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = 'black', face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1, 1, 1, 1),'cm'),
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),
        legend.frame = element_rect(colour = "black"),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.title = element_text(color = 'black', face = "bold", size=9))&
  scale_color_gradientn(colours = colorRampPalette(c("navy","white","firebrick3"))(100))&
  labs(tag = "")&
  theme(plot.tag.position = c(0.3, 1.05),
        plot.tag = element_text(size = 12,face = "bold"))&
  guides(size=guide_legend(title="Proportion of\nexpressing cells"),
         colour=guide_colorbar(title="Average\nexpression"))
#
ggsave(p1,filename="GATA_TF_across_cells.major_type.dotplot.pdf",width=6.5,height=6) ##Figure6F
q()

#################################################################################################################
######2.individual-level rs875179 genotype stratified correlation between GATA4/6 and NIPAL1 in epithelium#######
#################################################################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(data.table)
##
gtf <- rtracklayer::import('/data1/gy/public/gtf/gencode.v29lift37.annotation.gtf')
gtf <- as.data.frame(gtf)
gtf<- dplyr::select(gtf,c(gene_name,gene_id))
gtf$gene_id=substr(gtf$gene_id,1,15)
gtf<-unique(gtf)
gtf=subset(gtf,gene_name %in% c("GATA4","GATA5","GATA6","NIPAL1"))

##epi exp
type="epi"
exp_data <- fread(paste0("/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig/", type, "/", type, "/exp_", type, "_bed.bed.gz")) 
#
exp_data <- exp_data[, -c(1:3)]  
#
exp_data <- left_join(gtf, exp_data, by = c("gene_id" = "ID"))  
rownames(exp_data)=exp_data$gene_name
exp_data=exp_data[,-c(1:2)]

##
geno=fread(paste0("/data1/gy/sceQTL/genotype/", type, "/", type, "_CHN100K_filter05.dosage.txt"))
geno=subset(geno,ID=="4:48073300")
geno$ID[1]="rs875179"
#
sample_cols <- names(geno)[10:ncol(geno)]  
#
geno[, (sample_cols) := lapply(.SD, function(x) {  
  fcase(  
    x == 0, paste0(REF, REF),
    x == 1, paste0(REF, ALT),
    x == 2, paste0(ALT, ALT)
  )  
}), .SDcols = sample_cols] 
geno=as.data.frame(geno)
rownames(geno)=geno$ID
geno=geno[,-c(1:9)]

##
data=rbind(exp_data,geno)
data=as.data.frame(t(data))
data$sample=rownames(data)

#
data[, 1:4] <- lapply(data[, 1:4], as.numeric)
data$rs875179 = factor(data$rs875179, levels = c("CC", "CG", "GG"))
# 
data_long <- data %>%  
  select(NIPAL1, GATA4, GATA6, rs875179) %>%  
  pivot_longer(  
    cols = c("GATA4", "GATA6"),  
    names_to = "X_gene",  
    values_to = "X_value"  
  )  

###############plot
setwd("/data1/gy/ATAC_for_review/Figure6F&FigureS11G/output")
p_facet <- ggplot(data_long, aes(x = X_value, y = NIPAL1, color = rs875179)) +  
  geom_point() +  
  geom_smooth(method = "lm", se = TRUE) +  
  scale_color_manual(values = c("CC" = "#769FCD", "CG" = "#BCBCBC", "GG" = "#C86B85")) +  

  #
  stat_cor(
    method = "pearson",
    r.digits = 2, 
    p.digits = 3,
    label.x.npc = "left",
    label.y.npc = "top"
  ) +  

  facet_wrap(
    ~X_gene,
    scales = "free_x",
    strip.position = "bottom"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(
    y = "NIPAL1 Expression",
    color = "rs875179 genotype"
  )
#
ggsave(p_facet,filename="GATA4_GATA6_NIPAL1.scatterplot.pdf",width=9,height=5) ##FigureS11G
q()
