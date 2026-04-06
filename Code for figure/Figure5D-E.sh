#################################
###########Figure5D-E############
#################################
conda activate R_base
######################
########1.LD r2#######
######################
targetPeak="chr4:48072992-48073492"
ref="chr4:48090040:T:C"
#
plink --bfile /data1/gy/ATACseq_RWAS/RWAS_STITCH/ldref/1KG_EAS_no_indel_shared \
      --r2 \
      --ld-snp $ref \
      --ld-window-kb 400 \
      --ld-window 99999 \
      --ld-window-r2 0 \
      --out /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/LDr2/${targetPeak}

###############################
########2.locuszoom plot#######
###############################
conda activate R_base 
R
suppressMessages(library(data.table))  
suppressMessages(library(dplyr))  
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(data.table))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(ggbio))
suppressMessages(library(GenomicRanges))
suppressMessages(library(patchwork))
suppressMessages(library(stringr))  
suppressMessages(library(scales))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))

##
targetPeak <- "chr4:48072992-48073492"
ref1 <- "chr4:48090040:T:C"  ##rs11942098
ref2 <- "chr4:48073300:G:C"  ##rs875179
targetGene1 <- "NFXL1"
targetGeneid1 <- "ENSG00000170448"
targetGene2 <- "NIPAL1"
targetGeneid2 <- "ENSG00000163293"
outdir <- "/data1/gy/ATAC_for_review/Figure5D-E/output"

#
chr <- sub("^(chr\\d+):.*", "\\1", ref1)
ref1 <- sub("chr(\\d+):(\\d+):.*", "\\1:\\2", ref1)  
ref1_pos <- as.numeric(sub(".*:", "", ref1))
ref2 <- sub("chr(\\d+):(\\d+):.*", "\\1:\\2", ref2)  
ref2_pos <- as.numeric(sub(".*:", "", ref2))
#
start_end <- sub(".*:(\\d+)-(\\d+)", "\\1 \\2", targetPeak)  
start <- as.numeric(unlist(strsplit(start_end, " "))[1])  
end <- as.numeric(unlist(strsplit(start_end, " "))[2])  

##caQTL
caQTL <- fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt")  
caQTL = caQTL %>% filter(pheno_id == targetPeak)
caQTL <- caQTL %>% mutate(bp = str_extract(var_id, "(?<=:)[^:]+"))  
caQTL$pheno_chr <- substr(caQTL$pheno_chr, 4, nchar(caQTL$pheno_chr))  
caQTL$bp <- as.numeric(caQTL$bp)  
caQTL$chrbp = paste0(caQTL$pheno_chr,":",caQTL$bp)
caQTL$logP = -log10(caQTL$p_nominal)

##bulkeQTL
eQTL <- fread("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR.with_symbol.all.chr.txt")  
##NFXL1
eQTL1 = eQTL %>% filter(gene_name == targetGene1)
eQTL1[, c("chr", "bp", "Allele1", "Allele2") := tstrsplit(variant_id, ":", fixed = TRUE)]  
eQTL1$chrbp=paste0(eQTL1$chr,":",eQTL1$bp)
eQTL1$logP = -log10(eQTL1$pval_nominal)

##sceQTL epi
type="epi"
sceQTL <- fread(paste0("/data1/gy/sceQTL/eQTL_CHN100K/",type,"/",type,"/",type,"_sceQTL_result.allpairs.withsymbol.withFDR.txt"))
##NIPAL1
sceQTL2  = sceQTL %>% filter(gene_name == targetGene2)
sceQTL2 <- sceQTL2 %>% mutate(bp = str_extract(variant_id, "(?<=:)[^:]+"))   
sceQTL2 <- sceQTL2 %>% mutate(chr = str_extract(variant_id, "^[^:]+"))
sceQTL2$chrbp = paste0(sceQTL2$chr,":",sceQTL2$bp)
sceQTL2$chrbp = substr(sceQTL2$chrbp,4,nchar(sceQTL2$chrbp))
sceQTL2$bp =as.numeric(sceQTL2$bp) 
sceQTL2$logP = -log10(sceQTL2$pval_nominal)

##GWAS meta
GWAS <- fread("/data1/gy/EAS_GWAS_meta/GWAS_QCandFilter.new.metaResult")  
GWAS = subset(GWAS,chr==chr) 
colnames(GWAS)[10] <- "P_value"  
GWAS$chrbp=paste0(GWAS$chr,":",GWAS$bp)
GWAS$logP = -log10(GWAS$P_value)

##
Peak <- fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed", header = FALSE)  
Peak$name <- paste0(Peak$V1, ":", Peak$V2 + 1, "-", Peak$V3)  

##LDr2
LDr2 <- fread(paste0("/data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/", targetPeak, "/LDr2/",targetPeak,".ld"))  
LDr2$chrbp=paste0(LDr2$CHR_B,":",LDr2$BP_B)
LDr2 <- LDr2[, c(6,7,8)]  
color = as.character(cut(LDr2$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"), include.lowest=TRUE))
names(color) = LDr2$chrbp

##
caQTL <- subset(caQTL, chrbp %in% LDr2$chrbp)
eQTL1 <- subset(eQTL1, chrbp %in% LDr2$chrbp)  
eQTL2 <- subset(eQTL2, chrbp %in% LDr2$chrbp)  
sceQTL1 <- subset(sceQTL1, chrbp %in% LDr2$chrbp)  
sceQTL2 <- subset(sceQTL2, chrbp %in% LDr2$chrbp)  
GWAS <- subset(GWAS, chrbp %in% LDr2$chrbp)  

###
region_start <- start-200000
region_end <- end+200000

txdb <- AnnotationDbi::loadDb("/data1/gy/public/genome/hg19/txdb_v19_hg19.sqlite")
gr = GenomicRanges::GRanges(seqnames = chr, ranges = IRanges(region_start, region_end))
gr_Peak=with(Peak, GRanges(V1, IRanges(start = V2 + 1, end = V3, name = name)))  

##############Gene plot
p_gene <- ggplot() + theme_classic() +
  geom_alignment(
    txdb,
    which = gr,
    cds.rect.h = 0.1,
    color = "black",
    fill = "black",
    label.size = 3,
    arrow.rate = 0.015,
    length = unit(0.2, "cm"),
    gap.geom = 'arrow'
  ) +
  ylim(c(0.75,1.25)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x = element_text(size=14), axis.text.x = element_text(size=12)
  ) +ylab("") + xlab(paste0(chr,' (Mb)'))+
  scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(region_start,region_end)) +
  geom_vline(xintercept = ref1_pos, linetype=2) +
  geom_vline(xintercept = ref2_pos, linetype=2) +
  geom_point(aes(x=ref1_pos, y=1), shape=23, size=3, fill="purple")+
  geom_point(aes(x=ref2_pos, y=1), shape=23, size=3, fill="#D7261C")

##
gr_Peak_sub = gr_Peak[gr_Peak %within% gr]
#
df_loc = as.data.frame(gr_Peak_sub) 
df_loc$status = 0
df_loc$status[which(df_loc$start==start)] = 1
##color
col = c()
if( any(df_loc$status==0)){
    col = append(col, "black")
}
if( any(df_loc$status==1)){
    col = append(col, "red")
}
#
df_loc$status <- as.factor(df_loc$status)  

#############OCR plot
p_Peak = ggplot(df_loc) +   
  geom_segment(aes(x = start, y = 0, xend = end, yend = 0, color = status), size = 8) +   
  scale_x_continuous(expand = c(0, 0), limits = c(region_start, region_end)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +   
  scale_color_manual(values = col) +   
  theme_bw() +   
  theme(  
    legend.position = "none",   
    axis.title = element_blank(),   
    axis.ticks.x = element_blank(), 
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),   
    axis.line.y = element_blank(),   
    panel.grid = element_blank(),   
    axis.text.y = element_blank(),  
    panel.background = element_blank(),   
    strip.background = element_blank(),   
    rect = element_rect(fill = "white", linetype = 0)  
  ) 

###ABC anno
ABC=fread("/data1/gy/epistasis_enhancer/ABC_model_final/stomach_output/ABC_6sample_stomachHiC_mergedResult/6sample_merged_result/final_merged_result_allgene/6sample_merged_ABC_prediction_result_for_2repSample_withhead_final.txt")
ABC_subset <- ABC[chr == chr & start <= 48073492 & end >= 48072992]
ABC_subset
#      chr    start      end    class TargetGene isSelfPromoter
#   <char>    <int>    <int>   <char>     <char>         <lgcl>
#1:   chr4 48072965 48074066 enhancer     NIPAL1          FALSE
subset(ABC,TargetGene=="NIPAL1" & class=="promoter")
#      chr    start      end    class TargetGene isSelfPromoter
#   <char>    <int>    <int>   <char>     <char>         <lgcl>
#1:   chr4 48017691 48019785 promoter     NIPAL1           TRUE
anchor <- data.frame(
  region = c("P", "E"),
  xmin   = c(48017691, 48072965),
  xmax   = c(48019785, 48074066)
)

baseline   <- -0.05 
peak_space <- 0.55 

underline.df <- transform(anchor, y = baseline)

link.df <- data.frame(
  x    = rowMeans(anchor[anchor$region == "P", c("xmin", "xmax")]),
  xend = rowMeans(anchor[anchor$region == "E", c("xmin", "xmax")]),
  y    = baseline,
  yend = baseline
)

###############ABC anno plot
p_Peak_withABC <- p_Peak +
  #
  geom_segment(data = underline.df,
               aes(x = xmin, xend = xmax, y = y, yend = y),
               size   = 1,
               colour = "grey30",
               lineend = "round") +
  #
  geom_curve(data = link.df,
             aes(x = x, y = y, xend = xend, yend = yend),
             curvature = 0.35,
             angle     = 90,
             size      = 1.2,
             colour    = "#D1495B") +
  #
  scale_y_continuous(limits = c(baseline - 0.2, peak_space),
                     expand = c(0, 0))

############### caQTL plot
###
shape_caQTL = ifelse(caQTL$chrbp == ref1|caQTL$chrbp == ref2, 23, 21)
names(shape_caQTL) = caQTL$chrbp
###
size_caQTL = ifelse(caQTL$chrbp == ref1|caQTL$chrbp == ref2, 3, 2)
names(size_caQTL) = caQTL$chrbp
###
caQTL$label = ifelse(caQTL$chrbp == ref1|caQTL$chrbp == ref2, caQTL$var_id, '')
###
caQTL_title = paste0(targetPeak," caQTL")
p_caQTL = ggplot(caQTL,aes(x=bp,y=logP))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=caQTL[caQTL$label!='',],aes(x=bp,logP,fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape_caQTL,guide='none')+
          scale_size_manual(values=size_caQTL,guide='none')+
          scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(region_start,region_end))+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab(paste0(chr,' (Mb)'))+
          ylab(bquote(-log[10]*'(P)'))+
          labs(title = caQTL_title)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(caQTL_title, "\n") + 1)))
          ) + 
         scale_y_continuous(limits = c(-0.1, max(caQTL$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0, 0.05)))
p_caQTL =p_caQTL + geom_vline(xintercept = ref1_pos, linetype=2) +
          geom_point(aes(x=ref1_pos, y=caQTL[caQTL$bp ==ref1_pos,]$logP), shape=23, size=3, fill="purple")  
p_caQTL =p_caQTL + geom_vline(xintercept = ref2_pos, linetype=2) +
          geom_point(aes(x=ref2_pos, y=caQTL[caQTL$bp ==ref2_pos,]$logP), shape=23, size=3, fill="#D7261C")  

###############Bulk eQTL plot
##NFXL1
###
shape_eQTL1 = ifelse(eQTL1$chrbp == ref1|eQTL1$chrbp == ref2, 23, 21)
names(shape_eQTL1) = eQTL1$chrbp
###
size_eQTL1 = ifelse(eQTL1$chrbp == ref1|eQTL1$chrbp == ref2, 3, 2)
names(size_eQTL1) = eQTL1$chrbp
###
eQTL1$label = ifelse(eQTL1$chrbp == ref1|eQTL1$chrbp == ref2, eQTL1$variant_id, '')
###
eQTL1$bp=as.numeric(eQTL1$bp)
eQTL_title1 = paste0(targetGene1," bulk eQTL")
p_eQTL1 = ggplot(eQTL1,aes(x=bp,y=logP))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=eQTL1[eQTL1$label!='',],aes(x=bp,logP,fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape_eQTL1,guide='none')+
          scale_size_manual(values=size_eQTL1,guide='none')+
          scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(region_start,region_end))+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab(paste0(chr,' (Mb)'))+
          ylab(bquote(-log[10]*'(P)'))+
          labs(title = eQTL_title1)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(eQTL_title1, "\n") + 1)))
          ) + 
          scale_y_continuous(limits = c(-0.1, max(eQTL1$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0, 0.05)))
p_eQTL1 =p_eQTL1 + geom_vline(xintercept = ref1_pos, linetype=2) +
          geom_point(aes(x=ref1_pos, y=eQTL1[eQTL1$bp ==ref1_pos,]$logP), shape=23, size=3, fill="purple")
p_eQTL1 =p_eQTL1 + geom_vline(xintercept = ref2_pos, linetype=2) +
          geom_point(aes(x=ref2_pos, y=eQTL1[eQTL1$bp ==ref2_pos,]$logP), shape=23, size=3, fill="#D7261C")

############### sceQTL plot
##NIPAL1
###
shape_sceQTL2 = ifelse(sceQTL2$chrbp == ref1|sceQTL2$chrbp == ref2, 23, 21)
names(shape_sceQTL2) = sceQTL2$chrbp
###
size_sceQTL2 = ifelse(sceQTL2$chrbp == ref1|sceQTL2$chrbp == ref2, 3, 2)
names(size_sceQTL2) = sceQTL2$chrbp
###
sceQTL2$label = ifelse(sceQTL2$chrbp == ref1|sceQTL2$chrbp == ref2, sceQTL2$variant_id, '')
###
sceQTL_title2 = paste0(targetGene2," Epithelium sc-eQTL")
p_sceQTL2 = ggplot(sceQTL2,aes(x=bp,y=logP))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=sceQTL2[sceQTL2$label!='',],aes(x=bp,logP,fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape_sceQTL2,guide='none')+
          scale_size_manual(values=size_sceQTL2,guide='none')+
          scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(region_start,region_end))+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab(paste0(chr,' (Mb)'))+
          ylab(bquote(-log[10]*'(P)'))+
          labs(title = sceQTL_title2)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(sceQTL_title2, "\n") + 1)))
          ) + 
          scale_y_continuous(limits = c(-0.1, max(sceQTL2$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0, 0.05)))
p_sceQTL2 =p_sceQTL2 + geom_vline(xintercept = ref1_pos, linetype=2) +
          geom_point(aes(x=ref1_pos, y=sceQTL2[sceQTL2$bp ==ref1_pos,]$logP), shape=23, size=3, fill="purple")
p_sceQTL2 =p_sceQTL2 + geom_vline(xintercept = ref2_pos, linetype=2) +
          geom_point(aes(x=ref2_pos, y=sceQTL2[sceQTL2$bp ==ref2_pos,]$logP), shape=23, size=3, fill="#D7261C")   

###############GWAS plot
###
shape_GWAS = ifelse(GWAS$chrbp == ref1|GWAS$chrbp == ref2, 23, 21)
names(shape_GWAS) = GWAS$chrbp
###
size_GWAS = ifelse(GWAS$chrbp == ref1|GWAS$chrbp == ref2, 3, 2)
names(size_GWAS) = GWAS$chrbp
###
GWAS$label = ifelse(GWAS$chrbp == ref1|GWAS$chrbp == ref2, GWAS$chrbp, '')
###
GWAS_title = "GC GWAS"
p_GWAS = ggplot(GWAS,aes(x=bp,y=logP))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=GWAS[GWAS$label!='',],aes(x=bp,logP,fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape_GWAS,guide='none')+
          scale_size_manual(values=size_GWAS,guide='none')+
          scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(region_start,region_end))+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab(paste0(chr,' (Mb)'))+
          ylab(bquote(-log[10]*'(P)'))+
          labs(title = GWAS_title)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(GWAS_title, "\n") + 1)))
          ) + 
          scale_y_continuous(limits = c(-0.1, max(GWAS$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0, 0.05)))
p_GWAS =p_GWAS + geom_vline(xintercept = ref1_pos, linetype=2) +
          geom_point(aes(x=ref1_pos, y=GWAS[GWAS$bp ==ref1_pos,]$logP), shape=23, size=3, fill="purple") 
p_GWAS =p_GWAS + geom_vline(xintercept = ref2_pos, linetype=2) +
          geom_point(aes(x=ref2_pos, y=GWAS[GWAS$bp ==ref2_pos,]$logP), shape=23, size=3, fill="#D7261C") 

###############Plot togerher
p_together <- p_caQTL + scale_y_continuous(limits = c(-0.1, max(caQTL$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0.05, 0.15))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
              p_GWAS + scale_y_continuous(limits = c(-0.1, max(GWAS$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0.05, 0.15))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
              p_eQTL1 + scale_y_continuous(limits = c(-0.1, max(eQTL1$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0.05, 0.15))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
              p_sceQTL2 + scale_y_continuous(limits = c(-0.1, max(sceQTL2$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0.05, 0.15))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
              p_Peak_withABC +
              p_gene +
              plot_layout(nrow = 7, heights = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.4))
###add LDr2 legend
legend_box = data.frame(x = 0.9, y = seq(0.95, 0.87, -0.02))
p_together_withlegend = ggdraw(p_together)+
              geom_rect(data = legend_box,
              aes(xmin = x, xmax = x + 0.02, ymin = y, ymax = y + 0.02),
              color = "black",
              fill = rev(c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"))) +
              draw_label("1", x = legend_box$x[1] + 0.02, y = legend_box$y[1]+0.02, hjust = -0.7, size = 7) +
              draw_label("0.8", x = legend_box$x[1] + 0.02, y = legend_box$y[1], hjust = -0.3, size = 7) +
              draw_label("0.6", x = legend_box$x[2] + 0.02, y = legend_box$y[2], hjust = -0.3, size = 7) +
              draw_label("0.4", x = legend_box$x[3] + 0.02, y = legend_box$y[3], hjust = -0.3, size = 7) +
              draw_label("0.2", x = legend_box$x[4] + 0.02, y = legend_box$y[4], hjust = -0.3, size = 7) +
              draw_label("0", x = legend_box$x[5] + 0.02, y = legend_box$y[5], hjust = -0.7, size = 7) +
              draw_label(parse(text = "r^2"), x = legend_box$x[1]+0.01, y = legend_box$y[1], vjust = -2, size = 7)
ggplot2::ggsave(p_together_withlegend,filename=paste0(outdir,"/",targetPeak,"~",targetGene1,"-",targetGene2,".caQTL_GWAS_bulkeQTL_sceQTL.locuszoomplot.pdf"),width=220, height=280, units="mm") ##Figure5D-E
