###############################
#########FigureS3A-C###########
###############################
###############################################################
#######1.5sample scATACseq stomach major cell umap plot########FigureS3A
###############################################################
conda activate scATACseq_Signac1.9.0_final
#pip2 install macs2
#
R
setwd("/data1/gy/ATAC_for_review/FigureS3A-C/output")
set.seed(123456)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg19) ##hg19
library(ggplot2)
library(dplyr)
library(stringr)
library(rlang)
library(tidyr)
library(purrr)
library(tidyverse)
library(GenomicRanges)
library(IRanges)
library(ggpubr)
library(patchwork)
library(data.table)
library(biovizBase)
##
load('/data/blj/singlecell/ARC/int_signac_v2/integrated1.Rdata') #########published PMID: 40112817
##major cell type
use_colors=c(Epithelium = "#3182BDFF","Tcell" = "#FDAE6BFF","Bcell" = "#31A354FF",Myelocyte = "#756BB1FF",
  Mast = "#6BAED6FF",Fibroblasts = "#FD8D3CFF",ECs = "#74C476FF")
integrated1$type02 <- factor(integrated1$type02, levels=c('Epithelium','Tcell','Bcell','Myelocyte','Mast','Fibroblasts','ECs'), ordered=TRUE)
p <- DimPlot(integrated1, group.by = "type02", label = F)+
     scale_color_manual(values = use_colors)
ggsave('integrated1.major_cell_type.umapplot.pdf', p, width = 6.5, height = 5.5) ##FigureS3A


##########################################################
#######2.scATACseq call peak for 7 major cell type########
##########################################################
setwd("/data1/gy/scATACseq_multi/scATACseq_5Nsample_integrated/Signac_pipeline")
metadata <- data.frame(integrated1@meta.data)
metadata = metadata[,c("orig.ident","celltype02","type02")]
dir.create("./metadata")
write.csv(cbind("cell"=rownames(metadata),metadata),"./metadata/metadata.csv",row.names=F,quote=F)

##
umap_matrix = Embeddings(integrated1, reduction = "umap") 
dir.create("./umap_matrix")
write.csv(cbind("cell"=rownames(umap_matrix),umap_matrix),"./umap_matrix/umap_matrix.csv",row.names=F,quote=F)

#########call peak for 7 major cell type
#default: MACS2 -q0.05
DefaultAssay(integrated1) <- 'ATAC'
peaks_type <- CallPeaks(
  object = integrated1,
  group.by = "type02" 
)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions 
peaks_type <- keepStandardChromosomes(peaks_type, pruning.mode = "coarse")
#
old.seqnames <- seqlevels(blacklist_hg19)
new.seqnames <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
names(new.seqnames) <- old.seqnames
blacklist_hg19 <- renameSeqlevels(blacklist_hg19, new.seqnames)
peaks_type <- subsetByOverlaps(x = peaks_type, ranges = blacklist_hg19, invert = TRUE)
dir.create("./MACS2_call_peaks")
save(peaks_type,file = "./MACS2_call_peaks/MACS2_call_peaks_type.Rdata")
##
peaks_type_bed_df <- data.frame(  
                     chr = paste0("chr",as.character(seqnames(peaks_type))),  
                     start = start(peaks_type) - 1,  # 0-based
                     end = end(peaks_type),  
                     name = paste0(as.character(seqnames(peaks_type)),"-",as.character(start(peaks_type)),"-",as.character(end(peaks_type))),
                     type = peaks_type$peak_called_in
)  
write.table(  
  peaks_type_bed_df,   
  file = "./MACS2_call_peaks/MACS2_call_peaks_type.bed",  
  sep = "\t",  
  quote = FALSE,  
  col.names = FALSE,
  row.names = FALSE  
)  
##keep chr1-22
peaks_type_bed_df <- subset(peaks_type_bed_df,chr %in% paste0("chr",1:22)) 
write.table(  
  peaks_type_bed_df,   
  file = "./MACS2_call_peaks/MACS2_call_peaks_type.chr1_22.bed",  
  sep = "\t",  
  quote = FALSE, 
  col.names = FALSE,
  row.names = FALSE  
)  

#######################################################covariate for logistic regression (nCount_peaks)
DefaultAssay(integrated1) <- "ATAC" 
Idents(integrated1) <- "type02"

peaks_new <- read.table("./MACS2_call_peaks/MACS2_call_peaks_type.chr1_22.bed")
peaks_new$V1=substr(peaks_new$V1,4,nchar(peaks_new$V1))
peaks_new$V2=peaks_new$V2+1
peaks_new <- makeGRangesFromDataFrame(
  peaks_new,
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3",
  keep.extra.columns = F
)

frags_use <- Fragments(integrated1)[c(1, 7, 13, 19, 25)]
new_counts <- FeatureMatrix(
  fragments = frags_use,
  features = peaks_new,
  cells = colnames(integrated1),
  process_n = 2000  
)

chrom_assay_new <- CreateChromatinAssay(
  counts = new_counts,
  fragments = frags_use
)

integrated1[["ATAC_new"]] <- chrom_assay_new
DefaultAssay(integrated1) <- "ATAC_new"
integrated1$nCount_peaks_new <- colSums(integrated1[["ATAC_new"]]@counts)

save(integrated1, file = "/data1/gy/scATACseq_multi/scATACseq_5Nsample_integrated/Signac_pipeline/MACS2_call_peaks/integrated1_updated.RData")

##################################################################
#######3.5sample scATACseq stomach major celltype diffpeak########
##################################################################
setwd("/data1/gy/scATACseq_multi/scATACseq_5Nsample_integrated/Signac_pipeline")
################################################1-vs-rest[need ~4days]
types_todo <- c("Fibroblasts","Epithelium","Bcell","Tcell","Myelocyte","Mast","ECs")
markers_list_rest <- list()

for(tp in types_todo){
  message("Running 1-vs-rest for: ", tp)
  markers_list_rest[[tp]] <- FindMarkers(
    object = integrated1,
    ident.1 = tp,
    min.pct = 0.01,
    logfc.threshold = 0,
    test.use = "LR",
    latent.vars = "nCount_peaks_new"
  )
  #
  markers_list_rest[[tp]]$cluster <- tp
  fwrite(cbind("Peak"=rownames(markers_list_rest[[tp]]),markers_list_rest[[tp]]),
         paste0("./DARs_across_type/diff_peak.Findallmarker.pct0.01.",tp,".tsv"),sep="\t",col.names=T,row.names=F,quote=F)
}

######################################
#######4.scATACseq OCR barplot########FigureS3B
######################################
setwd("/data1/gy/scATACseq_multi/scATACseq_5Nsample_integrated/Signac_pipeline")
peaks_type_bed_df=fread("./MACS2_call_peaks/MACS2_call_peaks_type.chr1_22.bed")
colnames(peaks_type_bed_df)=c("chr","start","end","OCR","type")
##add celltype diff anno
for(celltype in c("Epithelium", "Tcell", "Bcell", "Myelocyte", "Mast", "Fibroblasts", "ECs")){
  celltype_diffpeak=fread(paste0("./DARs_across_type/diff_peak.Findallmarker.pct0.01.",celltype,".tsv"))
  celltype_diffpeak_sig=subset(celltype_diffpeak,p_val_adj<0.05 & avg_log2FC>0)
  peaks_type_bed_df[, paste0(celltype, "_diff") := OCR %in% celltype_diffpeak_sig$Peak]
}

##specific-peak for single cell type
diff_cols <- grep("_diff$", names(peaks_type_bed_df), value = TRUE)
celltypes <- sub("_diff$", "", diff_cols)

#
peaks_type_bed_df[, n_diff := rowSums(.SD), .SDcols = diff_cols]

#
for (ct in celltypes) {
  peaks_type_bed_df[, paste0(ct, "_specific") :=
    get(paste0(ct, "_diff")) == TRUE &
    n_diff == 1 &
    type == ct
  ]
}

#
peaks_type_bed_df[, n_diff := NULL]
fwrite(peaks_type_bed_df,"./MACS2_call_peaks/MACS2_call_peaks_type.with_diffanno_pct0.01.chr1_22.bed",sep="\t",col.names=T,row.names=F,quote=F)

#
library(stringr)
types <- unique(unlist(strsplit(peaks_type_bed_df$type, ",")))  
peaks_type_bed_df[, (types) := lapply(types, function(ct) as.integer(str_detect(type, ct)))]  
# 
library(ggplot2)  
library(scales) 
#
cell_order <- c('Epithelium','Tcell','Bcell','Myelocyte','Mast','Fibroblasts','ECs')  

###
cell_counts <- peaks_type_bed_df[  
  !is.na(type),  
  .(cell_type = factor(unlist(strsplit(gsub(" ", "", type), ",")), levels = cell_order)),  
  by = .(OCR)  
][, .(count = .N), by = cell_type][  
  order(match(cell_type, cell_order))] 

###
specific_cols <- grep("_specific$", names(peaks_type_bed_df), value = TRUE)
celltypes <- sub("_specific$", "", specific_cols)
cell_specific_counts <- rbindlist(lapply(seq_along(celltypes), function(i) {
  ct <- celltypes[i]
  col <- specific_cols[i]
  data.table(
    cell_type = ct,
    specific_count = sum(peaks_type_bed_df[[col]], na.rm = TRUE)
  )
}))

##
data=merge(cell_counts,cell_specific_counts)
data$cell_type=factor(data$cell_type,levels=cell_order)
#
cols = c(Epithelium = "#3182BDFF",Tcell = "#FDAE6BFF",Bcell = "#31A354FF",Myelocyte = "#756BB1FF",
         Mast = "#6BAED6FF",Fibroblasts = "#FD8D3CFF",ECs = "#74C476FF")

##barplot
setwd("/data1/gy/ATAC_for_review/FigureS3A-C/output")
##
p = ggplot(data, aes(x = cell_type)) +  
  #
  geom_col(
    aes(y = count, fill = cell_type),
    width = 0.7,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  #
  geom_text(
    aes(y = count, label = count),
    vjust = -0.6,
    size = 2.5,
    color = "black",
    lineheight = 0.8
  ) +
  # 
  geom_line(
    aes(y = specific_count, group = 1, color = "Celltype-specific OCR"),
    linewidth = 0.5
  ) +
  geom_point(
    aes(y = specific_count, color = "Celltype-specific OCR"),
    size = 3,
    shape = 19,
    show.legend = TRUE
  ) +
  #
  geom_text(
    aes(y = specific_count, label = specific_count),
    vjust = -0.8,  
    color = "black", 
    size = 2.5
  ) +
  #
  scale_y_continuous(
    name = "Number of sc-ATACseq OCRs of major cell types",
    expand = expansion(mult = c(0, 0.15))
  ) +
  scale_fill_manual(
    values = cols, 
    guide = "none"  
  ) +
  scale_color_manual(
    name = NULL,  
    values = "#EF6548",
    labels = "Celltype-specific OCR"
  ) +
  labs(x = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(1, 1),
    legend.justification = c(1, 1)
  )
ggsave(p,filename="scATACseq_5Nsample_integrated.OCR.type_anno.pertype.with_specific.barplot.pdf",width=6,height=4.5) ##FigureS3B

###########################################################################
#######5.5sample scATACseq stomach major cell specific peak heatmap########FigureS3C
###########################################################################
load("/data1/gy/scATACseq_multi/scATACseq_5Nsample_integrated/Signac_pipeline/MACS2_call_peaks/integrated1_updated.RData")
peaks=fread("/data1/gy/scATACseq_multi/scATACseq_5Nsample_integrated/Signac_pipeline/MACS2_call_peaks/MACS2_call_peaks_type.with_diffanno_pct0.01.chr1_22.bed")
Epithelium_OCRs <- peaks[Epithelium_specific == TRUE, OCR]
Tcell_OCRs <- peaks[Tcell_specific == TRUE, OCR]
Bcell_OCRs <- peaks[Bcell_specific == TRUE, OCR]
Myelocyte_OCRs <- peaks[Myelocyte_specific == TRUE, OCR]
Mast_OCRs <- peaks[Mast_specific == TRUE, OCR]
Fibroblasts_OCRs <- peaks[Fibroblasts_specific == TRUE, OCR]
ECs_OCRs <- peaks[ECs_specific == TRUE, OCR]
target_peaks=c(Epithelium_OCRs,Tcell_OCRs,Bcell_OCRs,Myelocyte_OCRs,Mast_OCRs,Fibroblasts_OCRs,ECs_OCRs)

##
DefaultAssay(integrated1) <- "ATAC_new"
#
Idents(integrated1) <- "type02"

#
avg_accessibility <- AverageExpression(
  integrated1,
  assays = "ATAC_new",
  features = target_peaks
)$ATAC_new

#
cell_order <- c("Epithelium", "Tcell", "Bcell","Myelocyte", "Mast", "Fibroblasts", "ECs")
avg_accessibility <- avg_accessibility[, cell_order, drop = FALSE]

#######################################heatmap
library(ComplexHeatmap)
library(circlize)

#
peak_labels <- rep(NA, length(target_peaks))

#
peak_labels[target_peaks %in% Epithelium_OCRs]   <- "Epithelium"
peak_labels[target_peaks %in% Tcell_OCRs]        <- "Tcell"
peak_labels[target_peaks %in% Bcell_OCRs]        <- "Bcell"
peak_labels[target_peaks %in% Myelocyte_OCRs]    <- "Myelocyte"
peak_labels[target_peaks %in% Mast_OCRs]         <- "Mast"
peak_labels[target_peaks %in% Fibroblasts_OCRs]  <- "Fibroblasts"
peak_labels[target_peaks %in% ECs_OCRs]          <- "ECs"

peak_info <- data.frame(
  OCR = target_peaks,
  specific_celltype = peak_labels,
  stringsAsFactors = FALSE
)
rownames(peak_info)=peak_info$OCR

#
mat <- t(avg_accessibility)
#
mat_z <- scale(mat, center = TRUE, scale = TRUE)

#
color_mapping <- colorRamp2(c(-3, 0, 3), c("#377EB8", "white", "#E41A1C")) 

#
type_colors <- c(Epithelium = "#3182BDFF",Tcell = "#FDAE6BFF",Bcell = "#31A354FF",Myelocyte = "#756BB1FF",
         Mast = "#6BAED6FF",Fibroblasts = "#FD8D3CFF",ECs = "#74C476FF")

#
top_anno <- HeatmapAnnotation(
  specific_celltype = peak_info$specific_celltype,
  col = list(specific_celltype = type_colors),
  annotation_label = "Cell type specificity",
  annotation_legend_param = list(
    specific_celltype = list(
      title = "Cell type specificity",
      at = c("Epithelium", "Tcell", "Bcell", "Myelocyte", "Mast", 
             "Fibroblasts", "ECs"),
      labels = c("Epithelium", "T cell", "B cell", "Myelocyte",
                 "Mast cell", "Fibroblast", "Endothelial cell")
    )
   ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  border = TRUE
)

# heatmap
p <- Heatmap(
  mat_z,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  col = color_mapping,
  top_annotation = top_anno,
  use_raster = TRUE, 
  border = "black",
  border_gp = gpar(col = "black", lwd = 2),
  raster_quality = 3,
  heatmap_legend_param = list(
    title = 'ATAC Zscore',
    at = c(-3, -1.5, 0, 1.5, 3),
    labels = c("-3", "-1.5", "0", "1.5", "3")
  )
)

##
pdf("scATACseq_5Nsample_integrated.7major_celltypes.specific_OCR.heatmap.pdf", width = 12, height = 4)  ##FigureS3C
draw(p)  
dev.off() 

