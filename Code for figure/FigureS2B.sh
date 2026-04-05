################################
############FigureS2B###########
################################
cd /data1/gy/software
wget https://ernstlab.biolchem.ucla.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM.zip
conda activate /Public/gaoyun/miniconda3/envs/ATACseq

##############################
########1.BinarizeBam#########
##############################
mkdir -p /data1/gy/chromHMM_stomach/binarydir
cd /data1/gy/chromHMM_stomach/binarydir
ChromHMM=/data1/gy/software/ChromHMM/ChromHMM.jar
chrlengthfile=/data1/gy/software/ChromHMM/CHROMSIZES/hg19.txt
inputdir_ENCODE=/data1/gy/public/stomach_ENCODE_region/bam
inputdir_mylab=/data1/gy/fresh_tissue_epigeno
cellmarkfiletable_ENCODE=/data1/gy/chromHMM_stomach/stomach_37Y-ENCODE_metadata_with_input.txt
cellmarkfiletable_mylab=/data1/gy/chromHMM_stomach/1443364-N_metadata_with_IgG.txt
outputdir=/data1/gy/chromHMM_stomach/binarydir
##ENCODE
java -mx4000M -jar ${ChromHMM} BinarizeBam \
-gzip \
-f 0 \
-g 0 \
-b 200 \
${chrlengthfile} ${inputdir_ENCODE} ${cellmarkfiletable_ENCODE} ${outputdir}

##NJ-GaEpi
java -mx4000M -jar ${ChromHMM} BinarizeBam \
-gzip \
-f 0 \
-g 0 \
-b 200 \
-paired \
${chrlengthfile} ${inputdir_mylab} ${cellmarkfiletable_mylab} ${outputdir}

############################
########2.LearnModel########
############################
ChromHMM=/data1/gy/software/ChromHMM/ChromHMM.jar
inputdir=/data1/gy/chromHMM_stomach/binarydir
outputdir=/data1/gy/chromHMM_stomach/segmentation
chrlengthfile=/data1/gy/software/ChromHMM/CHROMSIZES/hg19.txt
# Command
seq 8 15 | while read id
do
  java -mx4000M -jar ${ChromHMM} LearnModel \
  -b 200 -color 255,0,0 \
  -i stomach \
  -l ${chrlengthfile} \
  -s 123 \
  -noautoopen \
  -p 5 \
  ${inputdir} \
  ${outputdir}/state${id} \
  ${id} \
  hg19
done

####################################################
########3.plot chromHMM emission_11 heatmap#########
####################################################
/Public/gaoyun/software/R-4.2.0/bin/R
#
setwd("/data1/gy/ATAC_for_review/FigureS2B/output")
library(patchwork)
library(ComplexHeatmap)
###################################emissions_11
#
df<-read.delim('/data1/gy/chromHMM_stomach/segmentation/state11/emissions_11_stomach.txt',header = T)
df <- df[,c('State..Emission.order.', 'H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K27me3', 'H3K9me3')]
df <- df[c(2, 1, 7, 4, 6, 3, 8, 9, 11, 5, 10),]
df$State..Emission.order. <- c('Active TSS', 'Flanking active TSS', 'Poised promoter',
                               'Active enhancer', 'Weak enhancer', 'Poised enhancer',
                               'Repressed polycomb1', 'Repressed polycomb2',
                               'Heterochromatin',
                               'Quiescent1', 'Quiescent2')
df$class <- c(rep("Promoter",3),rep("Enhancer",3),rep("Repressed",2),rep("Heterochromatin",1),rep("Quiescent",2))
#
df$State..Emission.order.<-factor(df$State..Emission.order.,levels = rev(df$State..Emission.order.))
df$class<-factor(df$class,levels = rev(df$class))
#
rownames(df)=df[,1]
df=df[,-1]
df1<-df[,1:5,drop=F]
df2<-df[,6,drop=F]

#####
class_colors <- c("Promoter" = "#D73027", "Enhancer" = "#FDAE61", "Repressed" = "#3288BD", "Heterochromatin" = "#5E4FA2", "Quiescent" = "#878787")  
left_anno <- HeatmapAnnotation(  
    df = df2[, "class", drop = FALSE], 
    col = list(  
        class = class_colors  
    ),  
    which = "row",  
    annotation_label= " ",
    annotation_legend_param = list(   
        class = list(title = "class",
                     at = c("Promoter", "Enhancer", "Repressed", "Heterochromatin", "Quiescent")) 
    ),  
    annotation_name_side = "bottom",  
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),  
    border = TRUE,  
    show_annotation_name = TRUE  
)  
##heatmap
pdf("emissions_11_stomach.heatmap.pdf", width = 6.2, height = 7) ##FigureS2B
p = Heatmap(as.matrix(df1),
            cluster_rows = F ,
            cluster_columns = F ,
            show_column_dend = F,
            show_row_dend = F,
            show_column_names = T,
            show_row_names = T,
            row_title = NULL,
            column_title = NULL,
            row_gap = unit(0, "mm"), 
            left_annotation = left_anno,
            heatmap_legend_param = list(title = 'State'), 
            col = colorRampPalette(c("white", "#071ac8"))(100),
            column_names_gp = gpar(fontsize = 8),
            border = "black",
            border_gp = gpar(col = "black", lwd = 2)) 
p
dev.off()
