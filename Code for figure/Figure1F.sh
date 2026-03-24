###########################################
################Figure1F###################
###########################################
R
library(dplyr)
library(data.table)
##
stomach=fread("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/cluster_binary_file/pantissue.pantissue_OCR.IterativeOverlapPeakSet.TMM.below4_tissue_specific.stomach.limma_result_FDR0.001.binary_matrix.txt",header=T)
##
peak_anno=read.delim("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/iterative_peak_filtering/10.final_IterativeOverlapPeakSet_pantissue/pantissue_IterativeOverlapPeakSet.UCSC_hg19_known_gene.anno_by_chipseeker.with_symbol.txt",header=T)
peak_anno=peak_anno[,c(6,8,14,15,16)]
##
stomach_peak_anno=subset(peak_anno,V4 %in% colnames(stomach))
dim(stomach_peak_anno)
#[1] 9598    5
##
stomach_gene=unique(stomach_peak_anno$geneId)

#####################pathway enrichment
library(msigdbr)
library(GSEABase)
library(clusterProfiler)
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("tidyverse")
library(org.Hs.eg.db)
##GO:BP
erich.go.BP_stomach= enrichGO(gene =stomach_gene,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
erich.go.BP_stomach_result <- erich.go.BP_stomach@result
write.csv(erich.go.BP_stomach_result,"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/TMM_distal_heatmap/stomach_specific_pathway/enrich.go.BP.stomach_specific.csv")
#
library(enrichplot)
dotplot(erich.go.BP_stomach,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10) #top10
#
erich.go.BP_stomach_result <- erich.go.BP_stomach_result %>%
  mutate(
    GeneRatio_numerator = as.numeric(sub("/.*", "", GeneRatio)),
    GeneRatio_denominator = as.numeric(sub(".*?/(.*)", "\\1", GeneRatio)),
    BgRatio_numerator = as.numeric(sub("/.*", "", BgRatio)),
    BgRatio_denominator = as.numeric(sub(".*?/(.*)", "\\1", BgRatio)))
erich.go.BP_stomach_result$RichFactor=erich.go.BP_stomach_result$GeneRatio_numerator/erich.go.BP_stomach_result$GeneRatio_denominator
go_select <- erich.go.BP_stomach_result %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  head(10)  %>%
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

##
setwd("/data1/gy/ATAC_for_review/Figure1F/output")
library(ggplot2)
go_select$pval_log <- -log10(go_select$pvalue)
scale_factor <- max(go_select$pval_log) / max(go_select$RichFactor)
p <- ggplot(go_select, aes(y = reorder(Description, pval_log))) +
  geom_bar(aes(x = pval_log), stat = "identity", fill = "#7f4ea860") +
  geom_path(aes(x = RichFactor * scale_factor, y = reorder(Description, pval_log), group = 1),
            color = "black", size = 0.75) +
  geom_point(aes(x = RichFactor * scale_factor),
            shape = 21,
            fill = "#7f4ea8",
            color = "black",
            size = 3) +
  scale_x_continuous(
    name = "-log10(pvalue)",
    sec.axis = sec_axis(~ . / scale_factor, name = "GeneRatio")
  ) +
  labs(y = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(colour = "black", linewidth = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
ggsave("enrich.go.BP.stomach_specific.top10.barplot.pdf",p,width=7,height=6) ##Figure1F
