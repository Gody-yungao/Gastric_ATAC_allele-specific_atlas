############################
#########Figure1D###########
############################
R
#############################################
#######1.bulk OCR cell-type annotation#######
#############################################
##############################add specific anno
library(data.table)
library(dplyr)
library(GenomicRanges) 
macs2_type=fread("/data1/gy/scATACseq_multi/scATACseq_5Nsample_integrated/Signac_pipeline/MACS2_call_peaks/MACS2_call_peaks_type.with_diffanno_pct0.01.chr1_22.bed")
# find all *_specific
specific_cols <- grep("_specific$", names(macs2_type), value = TRUE)
#
cell_types <- sub("_specific$", "", specific_cols)
#
macs2_type[, specific_anno := "FALSE"]
#
for (i in seq_along(specific_cols)) {
  macs2_type[get(specific_cols[i]) == TRUE, specific_anno := cell_types[i]]
}
macs2_type=macs2_type[,c(1:5,20)]
table(macs2_type$specific_anno)
#      Bcell         ECs  Epithelium       FALSE Fibroblasts        Mast 
#       4209         908       22763      141967        8923          56 
#  Myelocyte       Tcell 
#        617        4865
colnames(macs2_type)=paste0("V",c(1:6))

##
OCR=fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed")
OCR$V4=paste0(OCR$V1,":",OCR$V2+1,"-",OCR$V3)
OCR=OCR[,1:4]

# 
OCR_gr <- OCR %>%   
  mutate(V2 = V2 + 1) %>%  # 0-based → 1-based 
  makeGRangesFromDataFrame(  
    seqnames.field = "V1",  
    start.field = "V2",  
    end.field = "V3",  
    keep.extra.columns = TRUE  
  )  

macs2_type_gr <- macs2_type %>%   
  mutate(V2 = V2 + 1) %>%
  makeGRangesFromDataFrame(  
    seqnames.field = "V1",  
    start.field = "V2",
    end.field = "V3",  
    keep.extra.columns = TRUE  
  )  
#
overlaps <- findOverlaps(OCR_gr, macs2_type_gr, maxgap = 0L)  
overlap_pairs <- data.table(  
  OCR_idx = queryHits(overlaps),  
  macs2_idx = subjectHits(overlaps)  
)  

#
calculate_overlap <- function(ocr_start, ocr_end, macs_start, macs_end) {  
  pmax(0, pmin(ocr_end, macs_end) - pmax(ocr_start, macs_start) + 1)  
}  

#
overlap_pairs[, overlap_length := calculate_overlap(  
  OCR[OCR_idx, V2+1], OCR[OCR_idx, V3],  
  macs2_type[macs2_idx, V2+1], macs2_type[macs2_idx, V3]  
)]  

#
annotations <- merge(  
  overlap_pairs,  
  macs2_type[, .(macs2_idx = .I, V1, V2, V3, V4, type = V5, type_specific = V6)],  
  by = "macs2_idx",  
  all.x = TRUE  
)  

#
best_annotations <- annotations[  
  order(-overlap_length),   
  .SD[1],   
  by = OCR_idx  
]  

#
OCR_annotated <- merge(  
  OCR[, .(V1, V2, V3, V4, OCR_idx = .I)],  
  best_annotations[, .(OCR_idx, type, type_specific, overlap_length, macs2_peakID = V4)],  
  by = "OCR_idx",  
  all.x = TRUE  
)[, OCR_idx := NULL]

dim(OCR_annotated)
#[1] 109187      7
colnames(OCR_annotated)[1:4]=c("chr","start","end","OCR")
dir.create("/data1/gy/ATACseq_RWAS/ATACseq/OCR_celltype_anno/macs2_type_peak_anno")
fwrite(OCR_annotated,"/data1/gy/ATACseq_RWAS/ATACseq/OCR_celltype_anno/macs2_type_peak_anno/95sample_IterativeOverlapPeakSet.type_anno.txt",col.names=T,row.names=F,quote=F,sep="\t")

#########################################
#####2.plot OCR celltype annobarplot#####
#########################################
library(ggplot2)  
library(scales)
library(ggpattern)  
setwd("/data1/gy/ATAC_for_review/Figure1D/output")
OCR_annotated=fread("/data1/gy/ATACseq_RWAS/ATACseq/OCR_celltype_anno/macs2_type_peak_anno/95sample_IterativeOverlapPeakSet.type_anno.txt")
#
cell_order <- c('Epithelium','Tcell','Bcell','Myelocyte','Mast','Fibroblasts','ECs')  

#
cell_counts <- OCR_annotated[
  !is.na(type),
  .(
    #
    cell_type = factor(unlist(strsplit(gsub(" ", "", type), ",")), levels = cell_order),
    
    # 
    is_independent = type_specific != "FALSE"
  ),
  by = .(OCR)
][,
  .(
    count = .N,                               
    independent_count = sum(is_independent)   
  ),
  by = cell_type
][,
  percent := count / nrow(OCR_annotated) 
][order(match(cell_type, cell_order))]

#
cell_counts[, independent_percent := independent_count / nrow(OCR_annotated)]

# 
max_count <- max(cell_counts$count)  
max_percent <- max(cell_counts$percent)  
y2_ratio <- max_count / max_percent 

#
cols = c(Epithelium = "#3182BDFF",Tcell = "#FDAE6BFF",Bcell = "#31A354FF",Myelocyte = "#756BB1FF",
         Mast = "#6BAED6FF",Fibroblasts = "#FD8D3CFF",ECs = "#74C476FF")
		 
##barplot
library(ggplot2)
library(scales)
#
p <- ggplot(na.omit(cell_counts), aes(x = cell_type)) +  
  geom_col(  
    aes(y = count, fill = cell_type),  
    width = 0.7,  
    alpha = 0.8,  
    show.legend = FALSE  
  ) +  
  geom_text(  
    aes(
      y = count,  
      label = sprintf("%d\n(%.1f%%)", count, percent*100)  
    ),
    vjust = -0.6, size = 2.8, color = "black", lineheight = 0.8   
  ) +
  geom_text(  
    aes(
      y = independent_count,  
      label = sprintf("%d\n(%.1f%%)", independent_count, independent_percent*100)  
    ),
    vjust = -0.6, size = 2.8, color = "black", lineheight = 0.8   
  ) +
  geom_line(
    aes(
      y = independent_percent * y2_ratio,
      group = 1,
      color = "Celltype-specific OCR"
    ),
    linewidth = 0.8,
    show.legend = TRUE
  ) +
  geom_point(
    aes(
      y = independent_percent * y2_ratio,
      color = "Celltype-specific OCR"
    ),
    size = 2,
    shape = 19,
    show.legend = TRUE
  ) +
  scale_y_continuous(
    name = "Number of ATAC-seq OCRs overlapped with\nOCRs of major cell type from scATAC-seq",  
    expand = expansion(mult = c(0, 0.15)),  
    sec.axis = sec_axis(
      ~ ./y2_ratio,  
      name = "Percentage of ATAC-seq OCRs overlapped with\nOCRs of major cell type from scATAC-seq",  
      labels = scales::label_percent(accuracy = 1),
      breaks = scales::pretty_breaks(n = 5)
    )
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
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 10),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    plot.margin = margin(10, 15, 10, 15),
    legend.position = c(1, 1),
    legend.justification = c(1, 1)
  )
ggsave(p,filename="95sample_IterativeOverlapPeakSet.type_anno.pertype.barplot.pdf",width=6,height=4.5) ##Figure1D
