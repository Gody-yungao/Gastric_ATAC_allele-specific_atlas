#########################
########Figure4D#########
#########################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
###################################
##########1.E-P link input#########
###################################
###TAD
TAD=fread("/data1/gy/epistasis_enhancer/HiC/TopDom_TAD/TAD_result/25000/chr1_22_stomach.25kb_res.250kb_win.TAD_result.domain.TSS_anno.ABCref_gencodev29.bed")
dim(TAD)
#[1] 6064    6
TAD$start=TAD$start+1

###loop
loop_sig=fread("/data1/gy/epistasis_enhancer/HiC/FitHiC_loop/fithicoutput_low_bound_20kb_res_10kb/FitHiC.spline_pass2.res10000.significances.qvalue0.05.TSS_anno.ABCref_gencodev29.txt")
dim(loop_sig)
#[1] 1236163      16
###
loop_sig[loop_sig == ""] <- NA
loop_sig <- loop_sig[!(is.na(loop_sig$anno1_gene_id) & is.na(loop_sig$anno2_gene_id)), ]
dim(loop_sig)
#[1] 160453     16
loop_sig$start1=loop_sig$start1+1
loop_sig$start2=loop_sig$start2+1

###ABC
ABC=fread("/data1/gy/epistasis_enhancer/ABC_model_final/stomach_output/ABC_6sample_stomachHiC_mergedResult/6sample_merged_result/final_merged_result_allgene/6sample_merged_ABC_prediction_result_for_2repSample_withhead_final.txt")
ABC=subset(ABC,class != "promoter")

####################################################################################
##########2.add HiC TAD & loop & ABC model anno to colocalized distal pairs#########
####################################################################################
###caQTL-bulk eQTL coloc
coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv")
coloc$Geneid=substr(coloc$Geneid,1,15)
#distal pair
coloc=subset(coloc,abs(coloc$Distance) > 1000)
#
coloc[, c("chr", "peak_start", "peak_end") := tstrsplit(Peak, "[:-]")]
coloc[, peak_start := as.integer(peak_start)]
coloc[, peak_end := as.integer(peak_end)]

#########################add loop anno
#
loop1 <- loop_sig[, .(loop_id = .I, chr = chr1, start = start1, end = end1, 
                      other_chr = chr2, other_start = start2, other_end = end2, 
                      loop_geneid = anno2_gene_id)]
loop2 <- loop_sig[, .(loop_id = .I, chr = chr2, start = start2, end = end2, 
                      other_chr = chr1, other_start = start1, other_end = end1, 
                      loop_geneid = anno1_gene_id)]
loop_long <- rbind(loop1, loop2)

#
setkey(loop_long, chr, start, end)
setkey(coloc, chr, peak_start, peak_end)

#
ovlp_loop <- foverlaps(
  coloc[, .(chr, peak_start, peak_end, Peak, Geneid)], 
  loop_long[, .(chr, start, end, loop_id, loop_geneid)],
  by.x = c("chr", "peak_start", "peak_end"),
  by.y = c("chr", "start", "end"),
  nomatch = 0L
)

##
ovlp_loop_filtered <- ovlp_loop[
  !is.na(loop_geneid) & sapply(
    seq_len(.N), 
    function(i) Geneid[i] %in% trimws(unlist(strsplit(loop_geneid[i], ",")))
  )
]
dim(ovlp_loop_filtered)
#[1] 237   9

# 
key_pairs <- unique(ovlp_loop_filtered[, .(Peak, Geneid)])
##
coloc[, loop_anno := fifelse(
  paste(Peak, Geneid) %in% paste(key_pairs$Peak, key_pairs$Geneid),
  TRUE, FALSE
)]

table(coloc$loop_anno)
#FALSE  TRUE 
# 1517   235

#########################add TAD anno
#
setkey(TAD, chr, start, end)
setkey(coloc, chr, peak_start, peak_end)

#
ovlp_TAD <- foverlaps(
  coloc[, .(chr, peak_start, peak_end, Peak, Geneid)], 
  TAD[, .(chr, start, end, anno_gene_id)],
  by.x = c("chr", "peak_start", "peak_end"),
  by.y = c("chr", "start", "end"),
  nomatch = 0L
)

##
ovlp_TAD_filtered <- ovlp_TAD[
  !is.na(anno_gene_id) & sapply(
    seq_len(.N), 
    function(i) Geneid[i] %in% trimws(unlist(strsplit(anno_gene_id[i], ",")))
  )
]
dim(ovlp_TAD_filtered)
#[1] 987   8

#
key_pairs <- unique(ovlp_TAD_filtered[, .(Peak, Geneid)])
##
coloc[, TAD_anno := fifelse(
  paste(Peak, Geneid) %in% paste(key_pairs$Peak, key_pairs$Geneid),
  TRUE, FALSE
)]
table(coloc$TAD_anno)
#FALSE  TRUE 
#  765   987

#########################add ABC model anno
#
setkey(ABC, chr, start, end)
setkey(coloc, chr, peak_start, peak_end)

#
ovlp_ABC <- foverlaps(
  coloc[, .(chr, peak_start, peak_end, Peak, Symbol)], 
  ABC[, .(chr, start, end, TargetGene)],
  by.x = c("chr", "peak_start", "peak_end"),
  by.y = c("chr", "start", "end"),
  nomatch = 0L
)

##
ovlp_ABC_filtered <- ovlp_ABC[
  TargetGene==Symbol
]
dim(ovlp_ABC_filtered)
#[1] 464   8

# 
key_pairs <- unique(ovlp_ABC_filtered[, .(Peak, Symbol)])
##
coloc[, ABC_anno := fifelse(
  paste(Peak, Symbol) %in% paste(key_pairs$Peak, key_pairs$Symbol),
  TRUE, FALSE
)]
table(coloc$ABC_anno)
#FALSE  TRUE 
# 1288   464

#########
coloc=coloc[,-c(14:16)]
dir.create("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno")
write.csv(coloc, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.all_pair.distal.TAD_loop_ABC_anno.csv", row.names = F,quote=F) 

########################################################################################
##########3.add HiC TAD & loop & ABC model anno to non-colocalized distal pairs#########
########################################################################################
###caQTL-bulk eQTL distal no coloc
coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4.chr1_22.csv")
coloc$Geneid=substr(coloc$Geneid,1,15)
coloc=subset(coloc,PPH4 < 0.5)
#
coloc=subset(coloc,abs(coloc$Distance) > 1000)
#
coloc[, c("chr", "peak_start", "peak_end") := tstrsplit(Peak, "[:-]")]
coloc[, peak_start := as.integer(peak_start)]
coloc[, peak_end := as.integer(peak_end)]
dim(coloc)
#[1] 192506     16

#########################add loop anno
#
loop1 <- loop_sig[, .(loop_id = .I, chr = chr1, start = start1, end = end1, 
                      other_chr = chr2, other_start = start2, other_end = end2, 
                      loop_geneid = anno2_gene_id)]
loop2 <- loop_sig[, .(loop_id = .I, chr = chr2, start = start2, end = end2, 
                      other_chr = chr1, other_start = start1, other_end = end1, 
                      loop_geneid = anno1_gene_id)]
loop_long <- rbind(loop1, loop2)

#
setkey(loop_long, chr, start, end)
setkey(coloc, chr, peak_start, peak_end)

#
ovlp_loop <- foverlaps(
  coloc[, .(chr, peak_start, peak_end, Peak, Geneid)], 
  loop_long[, .(chr, start, end, loop_id, loop_geneid)],
  by.x = c("chr", "peak_start", "peak_end"),
  by.y = c("chr", "start", "end"),
  nomatch = 0L
)

##
ovlp_loop_filtered <- ovlp_loop[
  !is.na(loop_geneid) & sapply(
    seq_len(.N), 
    function(i) Geneid[i] %in% trimws(unlist(strsplit(loop_geneid[i], ",")))
  )
]
dim(ovlp_loop_filtered)
#[1] 11916     9

#
key_pairs <- unique(ovlp_loop_filtered[, .(Peak, Geneid)])
##
coloc[, loop_anno := fifelse(
  paste(Peak, Geneid) %in% paste(key_pairs$Peak, key_pairs$Geneid),
  TRUE, FALSE
)]

table(coloc$loop_anno)
# FALSE   TRUE 
#180844  11662

#########################add TAD anno
#
setkey(TAD, chr, start, end)
setkey(coloc, chr, peak_start, peak_end)

#
ovlp_TAD <- foverlaps(
  coloc[, .(chr, peak_start, peak_end, Peak, Geneid)], 
  TAD[, .(chr, start, end, anno_gene_id)],
  by.x = c("chr", "peak_start", "peak_end"),
  by.y = c("chr", "start", "end"),
  nomatch = 0L
)

##
ovlp_TAD_filtered <- ovlp_TAD[
  !is.na(anno_gene_id) & sapply(
    seq_len(.N), 
    function(i) Geneid[i] %in% trimws(unlist(strsplit(anno_gene_id[i], ",")))
  )
]
dim(ovlp_TAD_filtered)
#[1] 37766     8

#
key_pairs <- unique(ovlp_TAD_filtered[, .(Peak, Geneid)])
##
coloc[, TAD_anno := fifelse(
  paste(Peak, Geneid) %in% paste(key_pairs$Peak, key_pairs$Geneid),
  TRUE, FALSE
)]
table(coloc$TAD_anno)
# FALSE   TRUE 
#154740  37766

#########################add ABC anno
#
setkey(ABC, chr, start, end)
setkey(coloc, chr, peak_start, peak_end)

#
ovlp_ABC <- foverlaps(
  coloc[, .(chr, peak_start, peak_end, Peak, Symbol)], 
  ABC[, .(chr, start, end, TargetGene)],
  by.x = c("chr", "peak_start", "peak_end"),
  by.y = c("chr", "start", "end"),
  nomatch = 0L
)

##
ovlp_ABC_filtered <- ovlp_ABC[
  TargetGene==Symbol
]
dim(ovlp_ABC_filtered)
#[1] 9344    8

#
key_pairs <- unique(ovlp_ABC_filtered[, .(Peak, Symbol)])
##
coloc[, ABC_anno := fifelse(
  paste(Peak, Symbol) %in% paste(key_pairs$Peak, key_pairs$Symbol),
  TRUE, FALSE
)]
table(coloc$ABC_anno)
# FALSE   TRUE 
#183178   9328

#########
coloc=coloc[,-c(14:16)]
write.csv(coloc, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno/caQTL_eQTL_coloc_200kb_1MB_PPH4_below0.5.all_pair.distal.TAD_loop_ABC_anno.csv", row.names = F,quote=F) 


########################################
##########4.enrichment analysis#########
########################################
distal_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.all_pair.distal.TAD_loop_ABC_anno.csv")
distal_coloc_control=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno/caQTL_eQTL_coloc_200kb_1MB_PPH4_below0.5.all_pair.distal.TAD_loop_ABC_anno.csv")
##
distal_coloc$group <- "colocalized"
distal_coloc_control$group <- "non-colocalized"
all_dat <- rbindlist(list(distal_coloc, distal_coloc_control), use.names = TRUE, fill = TRUE)
results=data.frame()

#Fisher's exact test
for (feature in c("loop_anno", "TAD_anno", "ABC_anno")) {
  # 
  tab <- table(factor(all_dat[[feature]], levels = c(TRUE, FALSE)), all_dat$group)
  fisher_test <- fisher.test(tab)
  results <- rbind(results, data.frame(FEATURE = feature,  
                                        Odds_Ratio = fisher_test$estimate,  
                                        P_Value = fisher_test$p.value,  
                                        CI_Lower = fisher_test$conf.int[1],  
                                        CI_Upper = fisher_test$conf.int[2]))  
}
#
dir.create("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno/enrichment")
write.csv(results, "/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno/enrichment/enrichment_result.csv", row.names = F,quote=F) 
#
head(results)
#     FEATURE Odds_Ratio CI_Lower CI_Upper       P_Value True_Colocalized
#      <fctr>      <num>    <num>    <num>         <num>            <int>
#1: loop_anno   2.402368 2.082222 2.761521  6.439799e-29              235
#2:  TAD_anno   5.286314 4.801919 5.821660 7.556484e-249              987
#3:  ABC_anno   7.075145 6.334200 7.888763 8.006161e-197              464
#   True_NonColocalized False_Colocalized False_NonColocalized
#                 <int>             <int>                <int>
#1:               11662              1517               180844
#2:               37766               765               154740
#3:                9328              1288               183178

##############################################
##########5.enrichment result barplot#########
##############################################
results=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_distal_result_anno/enrichment/enrichment_result.csv")

#
results[, FEATURE := factor(FEATURE, levels = FEATURE)]

#
results[, sig := fifelse(P_Value < 0.001, "***",
                  fifelse(P_Value < 0.01,  "**",
                  fifelse(P_Value < 0.05,  "*",  "ns")))]

#
rng <- range(results$CI_Upper, na.rm = TRUE)
y_margin <- diff(rng) * 0.05
results[, y_pos := CI_Upper + y_margin]

#
ann_dt <- results 

colors <- c(loop_anno="#FDD0A2", TAD_anno="#FDAE6B", ABC_anno="#FD8D3C")

setwd("/data1/gy/ATAC_for_review/Figure4D/output")
p <- ggplot(results, aes(x = FEATURE, y = Odds_Ratio, fill = FEATURE)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.15, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  #
  geom_text(data = ann_dt,
            aes(x = FEATURE, y = y_pos, label = sig),
            inherit.aes = FALSE, vjust = 0, size = 4) +
  #
  expand_limits(y = max(ann_dt$y_pos, na.rm = TRUE) * 1.05) +
  scale_fill_manual(values = colors, guide = "none") +
  labs(x = "Methods", y = "Odds Ratio (OR)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave(p,file="Colocalized_pairs.E-P_link.enrichment_result.barplot.pdf",width=4,height=5.5) ##Figure4D
