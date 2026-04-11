#######################
#######FigureS9A#######
#######################
###############GARFILED anno files were created in Figure3D.sh
############################################
#######1.GC GWAS-meta GARFIELD format#######
############################################
mkdir -p /data1/gy/EAS_GWAS_meta/garfield_format
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
GWAS <- fread("/data1/gy/EAS_GWAS_meta/GWAS_QCandFilter.inEAS.new.metaResult")  
GWAS=GWAS[,c("chr","bp","P-value")]   
for (m in 1:22) {
  GWAS_chr=subset(GWAS,chr==m)
  #
  GWAS_chr <- GWAS_chr %>%  
  arrange(bp)  
  fwrite(GWAS_chr[,2:3], paste0("/data1/gy/EAS_GWAS_meta/garfield_format/chr",m),col.names=F,row.names=F,sep=" ",quote=F)  
}
q()

#################################
#######2.garfield-prep-chr#######
#################################
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-prep-chr
cd /data1/gy/software/garfield-v2
GWAS=GC-meta
#
for i in {1..22}; do  
   #
   ./garfield-prep-chr \
   -ptags /data1/gy/public/garfield-data_EAS/tags/r01/chr$i \
   -ctags /data1/gy/public/garfield-data_EAS/tags/r08/chr$i \
   -maftss /data1/gy/public/garfield-data_EAS/maftssd/chr$i \
   -pval /data1/gy/EAS_GWAS_meta/garfield_format/chr$i \
   -ann /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/anno/chr$i \
   -o /data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-prep-chr/garfield.prep.$GWAS.out \
   -chr $i  
done

##################################
#######3.garfield-Meff-Padj#######
##################################
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-Meff-Padj
cd /data1/gy/software/garfield-v2
GWAS=GC-meta
/Public/gaoyun/software/R-4.2.0/bin/Rscript garfield-Meff-Padj.R \
-i /data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-prep-chr/garfield.prep.$GWAS.out \
-o /data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-Meff-Padj/garfield.Meff.$GWAS.out

##############################
#######4.garfield-test#######
##############################
mkdir -p /data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-test
/Public/gaoyun/software/R-4.2.0/bin/Rscript garfield-test.R \
-i /data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-prep-chr/garfield.prep.$GWAS.out \
-o /data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-test/garfield.test.$GWAS.out \
-l /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/garfield/anno/link_file.txt \
-pt 5e-2,5e-4,5e-6,5e-8 \
-b m5,n5,t5

###########################################
##########5.garfield-test results##########
###########################################
/Public/gaoyun/software/R-4.2.0/bin/R
dir.create("/data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield_enrichment_result")
setwd("/data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield_enrichment_result")
library(data.table)
library(dplyr)
data=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield-test/garfield.test.GC-meta.out")
data <- data %>%  
  arrange(desc(PThresh), ID) 
data$CI95_lower=exp(data$CI95_lower)
data$CI95_upper=exp(data$CI95_upper)
data <- data %>%
  mutate(OR = round(OR, 2)) %>%
  mutate(CI95_lower = round(CI95_lower, 2)) %>%
  mutate(CI95_upper = round(CI95_upper, 2))
data=data[,c(2,14,3,7,8)]
fwrite(data,"garfield_enrichment_result.for_forestplot.txt",sep="\t",quote=F,col.names=T)
##Manually compile and format it as an xlsx input file

############################################
##########6.enrichment forestplot##########
############################################
conda activate R_base
R
setwd("/data1/gy/ATACseq_RWAS/RWAS_STITCH/garfield/garfield_enrichment_result")
library(openxlsx)
library(dplyr)
data=read.xlsx("garfield_enrichment_result.for_forestplot.xlsx")
##
data= data %>%
  mutate(id = factor(letters[1:nrow(data)]))
#
data <- data %>%  
  mutate(color = case_when(  
    Subgroups == "caQTL" ~ "red",  
    Subgroups == "ASCA" ~ "orange",  
    Subgroups == "combined" ~ "pink",  
    TRUE ~ NA_character_
  ))  

##
setwd("/data1/gy/ATAC_for_review/FigureS9A/output")
library(ggplot2)
p1=ggplot(data, aes(x = log2(OR), y = rev(id))) +  
        geom_point(aes(color = color), size = 3) +
        geom_errorbarh(aes(xmin = log2(CI95_lower), xmax = log2(CI95_upper), color = color), height = 0.2) +
        scale_color_manual(values = c("red" = "#EF3B2C", "orange" = "#EC7014", "pink" = "#DD3497"))+
        labs(title = "",  
             x = "OR",  
             y = "") +  
        scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5),
                           labels = c("1", "2", "4", "8", "16", "32"),
                           limits = c(0, NA),
                           expand = c(0, 0)) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) + 
        theme_minimal() +  
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.ticks.length = unit(0.25, "cm"),
              axis.ticks.x = element_line(size = 0.3),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.line.y = element_blank())
##
p2 <- ggplot(data,aes(x=1,y=rev(id)))+
  geom_text(aes(label=Subgroups),hjust=0.5,size=3)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(rep(0,4),"cm")
        )+
  labs(x="Subgroups")+
  scale_x_discrete(position = "top")+
  theme(axis.title.x = element_text(hjust = 0.5))
##
p3 <- ggplot(data,aes(x=1,y=rev(id)))+
  geom_text(aes(label=`OR(95%CI)`),hjust=0.5,size=3)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(rep(0,4),"cm")
        )+
  labs(x="OR(95%CI)")+
  scale_x_discrete(position = "top")+
  theme(axis.title.x = element_text(hjust = 0.5))

##
library(patchwork)
p=p2+p3+p1+plot_layout(widths = c(0.1,0.1,0.3))
ggsave(p,filename="garfield_enrichment_result.forestplot.pdf") ##FigureS9A
