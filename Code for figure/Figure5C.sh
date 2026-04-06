##########################
#########Figure5C#########
##########################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
RWAS=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter/sig/RWAS.fdr0.1.sig.result.txt")
TWAS=fread("/data1/gy/ATACseq_RWAS/TWAS_262/result/TWAS.withFDR.FDR0.1.result.txt")
TWAS$TSS=TWAS$P1-1e6
epiTWAS=fread("/data1/gy/ATACseq_RWAS/epiTWAS_203/result/epiTWAS.withFDR.FDR0.1.result.txt")
epiTWAS$TSS=epiTWAS$P1-1e6

####################################################
########1.bulk TWAS add RWAS OCR anno(±1MB)#########
####################################################
library(data.table)
library(dplyr)
library(GenomicRanges)
library(stringr)
library(purrr)

# bulkTWAS: columns CHR, TSS, TWAS.Z, ID/gene_name
# RWAS   : columns CHR, P0, P1, TWAS.Z, ID
bulk <- as.data.table(TWAS)
rwas <- as.data.table(RWAS)

#
norm_chr <- function(x) {
  x <- as.character(x)
  x <- ifelse(str_detect(x, "^chr"), x, paste0("chr", x))
  x
}
bulk[, CHR := norm_chr(CHR)]
rwas[, CHR := norm_chr(CHR)]

#
stopifnot(all(c("CHR","TSS","TWAS.Z") %in% names(bulk)))
stopifnot(all(c("CHR","P0","P1","TWAS.Z") %in% names(rwas)))

##bulk TWAS gene （TSS ±1000 kb）
pad <- 1000000L
bulk_gr <- GRanges(
  seqnames = bulk$CHR,
  ranges   = IRanges(start = pmax(bulk$TSS - pad, 1L), end = bulk$TSS + pad),
  idx      = seq_len(nrow(bulk)),
  gene_id  = if ("ID" %in% names(bulk)) bulk$ID else if ("gene_name" %in% names(bulk)) bulk$gene_name else NA_character_,
  twas_z   = bulk$TWAS.Z
)

##RWAS OCR
rwas_gr <- GRanges(
  seqnames = rwas$CHR,
  ranges   = IRanges(start = rwas$P0, end = rwas$P1),
  ocr_id   = rwas$ID,
  rwas_z   = rwas$TWAS.Z,
  p0       = rwas$P0,
  p1       = rwas$P1
)

#overlap
hits <- findOverlaps(bulk_gr, rwas_gr, ignore.strand = TRUE)

# 
annot <- data.table(
  bulk_idx = mcols(bulk_gr)$idx[queryHits(hits)],
  gene_id = mcols(bulk_gr)$gene_id[queryHits(hits)],
  bulk_TWAS_Z = mcols(bulk_gr)$twas_z[queryHits(hits)],
  ocr_id = mcols(rwas_gr)$ocr_id[subjectHits(hits)],
  rwas_TWAS_Z = mcols(rwas_gr)$rwas_z[subjectHits(hits)],
  ocr_chr = as.character(seqnames(rwas_gr))[subjectHits(hits)],
  ocr_p0 = mcols(rwas_gr)$p0[subjectHits(hits)],
  ocr_p1 = mcols(rwas_gr)$p1[subjectHits(hits)]
)

#Directional consistency filtering: bulkTWAS and RWAS TWAS.Z scores have the same sign
annot_dir <- annot %>%
  mutate(direction_consistent = (bulk_TWAS_Z * rwas_TWAS_Z) > 0) %>%
  filter(direction_consistent)

#
if (nrow(annot_dir) == 0L) {
  bulk$OCR_within_1MB_dirCons <- NA_character_
  bulk$OCR_within_1MB_dirCons_n <- 0L
} else {
  #
  ocr_label <- annot_dir %>%
    mutate(ocr_region = paste0(ocr_chr, ":", ocr_p0, "-", ocr_p1,
                               " (", ocr_id, ", Z=", signif(rwas_TWAS_Z, 3), ")")) %>%
    arrange(bulk_idx, ocr_chr, ocr_p0) %>%
    group_by(bulk_idx) %>%
    summarise(
      OCR_within_1MB_dirCons = paste(ocr_region, collapse = "; "),
      OCR_within_1MB_dirCons_n = n(),
      .groups = "drop"
    )
  
  #
  bulk <- bulk %>%
    mutate(rowid = row_number()) %>%
    left_join(ocr_label, by = c("rowid" = "bulk_idx")) %>%
    select(-rowid)
}

#
bulk <- bulk %>%
  mutate(has_OCR_dirCons_1MB = OCR_within_1MB_dirCons_n > 0)

#
bulk_result <- bulk
data.table::setDT(bulk_result)
bulk_result[]
fwrite(bulk_result,"/data1/gy/ATACseq_RWAS/TWAS_262/result/TWAS.withFDR.FDR0.1.result.withRWASanno.txt",col.names=T,row.names=F,quote=F,sep="\t")
##
table(bulk_result$has_OCR_dirCons_1MB)
#TRUE NA
#  24 21

#################################################
########2.epiTWAS add RWAS OCR anno(±1MB)#########
#################################################
# epiTWAS: columns CHR, TSS, TWAS.Z, ID/gene_name
# RWAS   : columns CHR, P0, P1, TWAS.Z, ID

epi <- as.data.table(epiTWAS)
rwas <- as.data.table(RWAS)

#
norm_chr <- function(x) {
  x <- as.character(x)
  x <- ifelse(str_detect(x, "^chr"), x, paste0("chr", x))
  x
}
epi[, CHR := norm_chr(CHR)]
rwas[, CHR := norm_chr(CHR)]

#
stopifnot(all(c("CHR","TSS","TWAS.Z") %in% names(epi)))
stopifnot(all(c("CHR","P0","P1","TWAS.Z") %in% names(rwas)))

##epi TWAS gene （TSS ±1000 kb）
pad <- 1000000L
epi_gr <- GRanges(
  seqnames = epi$CHR,
  ranges   = IRanges(start = pmax(epi$TSS - pad, 1L), end = epi$TSS + pad),
  idx      = seq_len(nrow(epi)),
  gene_id  = if ("ID" %in% names(epi)) epi$ID else if ("gene_name" %in% names(epi)) epi$gene_name else NA_character_,
  twas_z   = epi$TWAS.Z
)

##RWAS OCR
rwas_gr <- GRanges(
  seqnames = rwas$CHR,
  ranges   = IRanges(start = rwas$P0, end = rwas$P1),
  ocr_id   = rwas$ID,
  rwas_z   = rwas$TWAS.Z,
  p0       = rwas$P0,
  p1       = rwas$P1
)

#overlap
hits <- findOverlaps(epi_gr, rwas_gr, ignore.strand = TRUE)

#
annot <- data.table(
  epi_idx = mcols(epi_gr)$idx[queryHits(hits)],
  gene_id = mcols(epi_gr)$gene_id[queryHits(hits)],
  epi_TWAS_Z = mcols(epi_gr)$twas_z[queryHits(hits)],
  ocr_id = mcols(rwas_gr)$ocr_id[subjectHits(hits)],
  rwas_TWAS_Z = mcols(rwas_gr)$rwas_z[subjectHits(hits)],
  ocr_chr = as.character(seqnames(rwas_gr))[subjectHits(hits)],
  ocr_p0 = mcols(rwas_gr)$p0[subjectHits(hits)],
  ocr_p1 = mcols(rwas_gr)$p1[subjectHits(hits)]
)

#Directional consistency filtering: epiTWAS and RWAS TWAS.Z scores have the same sign
annot_dir <- annot %>%
  mutate(direction_consistent = (epi_TWAS_Z * rwas_TWAS_Z) > 0) %>%
  filter(direction_consistent)

#
if (nrow(annot_dir) == 0L) {
  epi$OCR_within_1MB_dirCons <- NA_character_
  epi$OCR_within_1MB_dirCons_n <- 0L
} else {
  #
  ocr_label <- annot_dir %>%
    mutate(ocr_region = paste0(ocr_chr, ":", ocr_p0, "-", ocr_p1,
                               " (", ocr_id, ", Z=", signif(rwas_TWAS_Z, 3), ")")) %>%
    arrange(epi_idx, ocr_chr, ocr_p0) %>%
    group_by(epi_idx) %>%
    summarise(
      OCR_within_1MB_dirCons = paste(ocr_region, collapse = "; "),
      OCR_within_1MB_dirCons_n = n(),
      .groups = "drop"
    )
  
  #
  epi <- epi %>%
    mutate(rowid = row_number()) %>%
    left_join(ocr_label, by = c("rowid" = "epi_idx")) %>%
    select(-rowid)
}

#
epi <- epi %>%
  mutate(has_OCR_dirCons_1MB = OCR_within_1MB_dirCons_n > 0)

#
epi_result <- epi
data.table::setDT(epi_result)
epi_result[]
fwrite(epi_result,"/data1/gy/ATACseq_RWAS/epiTWAS_203/result/epiTWAS.withFDR.FDR0.1.result.withRWASanno.txt",col.names=T,row.names=F,quote=F,sep="\t")
##
table(epi_result$has_OCR_dirCons_1MB)
#TRUE NA
#   9  1

######################################################
########3.combine bulk&epi TWAS RWAS OCR anno#########
######################################################
bulk_sorted <- bulk_result %>%
  #
  mutate(CHR = as.character(CHR)) %>%
  #
  mutate(
    chr_num = case_when(
      str_detect(CHR, "^chr") ~ str_remove(CHR, "^chr"),
      TRUE ~ CHR
    ),
    chr_num = case_when(
      chr_num %in% as.character(1:22) ~ as.integer(chr_num),
      chr_num %in% c("X","x") ~ 23L,
      chr_num %in% c("Y","y") ~ 24L,
      chr_num %in% c("M","MT","m","mt") ~ 25L,
      TRUE ~ 99L
    )
  ) %>%
  #
  arrange(chr_num, TSS) %>%
  select(gene_name, chr_num, TSS, has_OCR_dirCons_1MB)
colnames(bulk_sorted)[4]="bulkTWAS_has_OCR_dirCons_1MB"
bulk_sorted$bulkTWAS_gene=TRUE
##
epi_sorted <- epi_result %>%
  #
  mutate(CHR = as.character(CHR)) %>%
  #
  mutate(
    chr_num = case_when(
      str_detect(CHR, "^chr") ~ str_remove(CHR, "^chr"),
      TRUE ~ CHR
    ),
    chr_num = case_when(
      chr_num %in% as.character(1:22) ~ as.integer(chr_num),
      chr_num %in% c("X","x") ~ 23L,
      chr_num %in% c("Y","y") ~ 24L,
      chr_num %in% c("M","MT","m","mt") ~ 25L,
      TRUE ~ 99L
    )
  ) %>%
  #
  arrange(chr_num, TSS) %>%
  select(gene_name, chr_num, TSS, has_OCR_dirCons_1MB)
colnames(epi_sorted)[4]="epiTWAS_has_OCR_dirCons_1MB"
epi_sorted$epiTWAS_gene=TRUE
##
res=full_join(bulk_sorted,epi_sorted,by="gene_name")
##
res_merged <- res %>%
  mutate(
    chr_num = coalesce(chr_num.x, chr_num.y),
    TSS     = coalesce(TSS.x, TSS.y)
  ) %>%
  select(
    #
    gene_name,
    chr_num, TSS,
    bulkTWAS_gene, epiTWAS_gene,
    bulkTWAS_has_OCR_dirCons_1MB,
    epiTWAS_has_OCR_dirCons_1MB,
  ) %>%
  arrange(chr_num, TSS)
##
res_merged2 <- res_merged %>%
  mutate(
    has_OCR_dirCons_1MB = coalesce(bulkTWAS_has_OCR_dirCons_1MB, FALSE) |
                          coalesce(epiTWAS_has_OCR_dirCons_1MB,  FALSE)
  )%>%
  select(
    #
    gene_name,
    chr_num, TSS,
    bulkTWAS_gene, epiTWAS_gene,
    has_OCR_dirCons_1MB
  )
##
dir.create("/data1/gy/ATACseq_RWAS/RWAS_STITCH/TWAS_epiTWAS_withRWAS")
fwrite(res_merged2,"/data1/gy/ATACseq_RWAS/RWAS_STITCH/TWAS_epiTWAS_withRWAS/TWAS_epiTWAS.withRWASanno.dirCons_1MB.txt",col.names=T,row.names=F,quote=F,sep="\t")
##
dim(res_merged2)
#[1] 52  6

###########################
########3.heatmap##########
###########################
res_merged2=fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/TWAS_epiTWAS_withRWAS/TWAS_epiTWAS.withRWASanno.dirCons_1MB.txt")
# 
res_merged2[, bulk2 := fcoalesce(bulkTWAS_gene, FALSE)]
res_merged2[, epi2  := fcoalesce(epiTWAS_gene,  FALSE)]

#
res_merged2[, grp := fifelse(bulk2 & !epi2, 1L,
                       fifelse(bulk2 &  epi2, 2L, 3L))]

#
setorder(res_merged2, grp, chr_num, TSS)

#
res_merged2[, c("bulk2","epi2","grp") := NULL]

#
res_num <- res_merged2 %>%
  mutate(
    across(
      c(bulkTWAS_gene, epiTWAS_gene),
      ~ ifelse(is.na(.), 0L, ifelse(. == TRUE, 1L, 0L))
    )
  ) %>%
  mutate(
    across(
      c(has_OCR_dirCons_1MB),
      ~ ifelse(is.na(.), 0L, ifelse(. == TRUE, 2L, 0L))
    )
  )
res_num=as.data.frame(res_num)
rownames(res_num)=res_num$gene_name
res_num=res_num[,4:6]
plot_df=t(res_num)

###heatmap
setwd("/data1/gy/ATAC_for_review/Figure5C/output")
library(ComplexHeatmap)  
plot <- ComplexHeatmap::Heatmap(
    as.matrix(plot_df),
    col = c("white", "#145390", "#F08E59"),
    show_heatmap_legend = FALSE,
    rect_gp = grid::gpar(col = "white", lwd = 0.5),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_names_gp = grid::gpar(fontsize = 6, angle=45, fontface="italic"),
    column_names_side = "top",               
    row_names_gp = grid::gpar(fontsize = 8),
    row_names_side = "left",
    column_names_rot = 90,
    border = "black", 
    border_gp = gpar(col = "black", lwd = 0.5))

pdf("TWAS_epiTWAS.withRWASanno.dirCons_1MB.heatmap.pdf", width = 6, height = 2) ##Figure5C
plot
dev.off()
