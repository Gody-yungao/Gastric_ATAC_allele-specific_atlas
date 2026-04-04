##################################
#############Figure4F#############
##################################
######################
########1.LD r2#######
######################
conda activate R_base
targetPeak="chr14:34896054-34896554"
ref="chr14:34895977:G:A" ##rs10132233
targetGene="SPTSSA"
# 
plink --bfile /data1/gy/ATACseq_RWAS/RWAS_STITCH/ldref/1KG_EAS_no_indel \
      --r2 \
      --ld-snp $ref \
      --ld-window-kb 400 \
      --ld-window 99999 \
      --ld-window-r2 0 \
      --out /data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/${targetPeak}/LDr2/${targetPeak}

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

##
targetPeak <- "chr14:34896054-34896554"
ref <- "chr14:34895977:G:A" ##rs10132233
targetGene <- "SPTSSA"
targetGeneid <- "ENSG00000165389"
outdir <- "/data1/gy/ATAC_for_review/Figure4F/output"

#
chr <- sub("^(chr\\d+):.*", "\\1", ref)
ref <- sub("chr(\\d+):(\\d+):.*", "\\1:\\2", ref)  
ref_pos <- as.numeric(sub(".*:", "", ref))

#
start_end <- sub(".*:(\\d+)-(\\d+)", "\\1 \\2", targetPeak)  
start <- as.numeric(unlist(strsplit(start_end, " "))[1])  
end <- as.numeric(unlist(strsplit(start_end, " "))[2])  

#
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
eQTL = eQTL %>% filter(gene_name == targetGene)
eQTL[, c("chr", "bp", "Allele1", "Allele2") := tstrsplit(variant_id, ":", fixed = TRUE)]  
eQTL$chrbp=paste0(eQTL$chr,":",eQTL$bp)
eQTL$logP = -log10(eQTL$pval_nominal)

##sceQTL epi
type="epi"
sceQTL <- fread(paste0("/data1/gy/sceQTL/eQTL_CHN100K/",type,"/",type,"/",type,"_sceQTL_result.allpairs.withsymbol.withFDR.txt"))
sceQTL  = sceQTL %>% filter(gene_name == targetGene)
sceQTL <- sceQTL %>% mutate(bp = str_extract(variant_id, "(?<=:)[^:]+"))   
sceQTL <- sceQTL %>% mutate(chr = str_extract(variant_id, "^[^:]+"))
sceQTL$chrbp = paste0(sceQTL$chr,":",sceQTL$bp)
sceQTL$chrbp = substr(sceQTL$chrbp,4,nchar(sceQTL$chrbp))
sceQTL$bp =as.numeric(sceQTL$bp) 
sceQTL$logP = -log10(sceQTL$pval_nominal)

##
Peak <- fread("/data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed", header = FALSE)  
Peak$name <- paste0(Peak$V1, ":", Peak$V2 + 1, "-", Peak$V3)  

#
LDr2 <- fread(paste0("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_sceQTL_coloc_celltype_specific_OCR/", targetPeak, "/LDr2/",targetPeak,".ld"))  
LDr2$chrbp=paste0(LDr2$CHR_B,":",LDr2$BP_B)
LDr2 <- LDr2[, c(6,7,8)]  
color = as.character(cut(LDr2$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c("#1A2B58", "#82BBE1", "#67AD32", "#E89920", "#D7261C"), include.lowest=TRUE))
names(color) = LDr2$chrbp

# 
caQTL <- subset(caQTL, chrbp %in% LDr2$chrbp)
eQTL <- subset(eQTL, chrbp %in% LDr2$chrbp)  
sceQTL <- subset(sceQTL, chrbp %in% LDr2$chrbp)  

###
region_start <- start-200000
region_end <- end+200000

txdb <- AnnotationDbi::loadDb("/data1/gy/public/genome/hg19/txdb_v19_hg19.sqlite")
gr = GenomicRanges::GRanges(seqnames = chr, ranges = IRanges(region_start, region_end))
gr_Peak=with(Peak, GRanges(V1, IRanges(start = V2 + 1, end = V3, name = name)))  

########Gene plot
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
  geom_vline(xintercept = ref_pos, linetype=2) +
  geom_point(aes(x=ref_pos, y=1), shape=23, size=3, fill="purple")

############Peak plot
gr_Peak_sub = gr_Peak[gr_Peak %within% gr]
# convert to data.frame
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
##
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

###############caQTL plot
###
shape_caQTL = ifelse(caQTL$chrbp == ref, 23, 21)
names(shape_caQTL) = caQTL$chrbp
###
size_caQTL = ifelse(caQTL$chrbp == ref, 3, 2)
names(size_caQTL) = caQTL$chrbp
###
caQTL$label = ifelse(caQTL$chrbp == ref, caQTL$var_id, '')
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
p_caQTL =p_caQTL + geom_vline(xintercept = ref_pos, linetype=2) +
          geom_point(aes(x=ref_pos, y=caQTL[caQTL$bp ==ref_pos,]$logP), shape=23, size=3, fill="purple")  
		  
############### eQTL plot
###
shape_eQTL = ifelse(eQTL$chrbp == ref, 23, 21)
names(shape_eQTL) = eQTL$chrbp
###
size_eQTL = ifelse(eQTL$chrbp == ref, 3, 2)
names(size_eQTL) = eQTL$chrbp
###
eQTL$label = ifelse(eQTL$chrbp == ref, eQTL$variant_id, '')
###
eQTL_title = paste0(targetGene," bulk eQTL")
eQTL$bp=as.numeric(eQTL$bp)
p_eQTL = ggplot(eQTL,aes(x=bp,y=logP))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=eQTL[eQTL$label!='',],aes(x=bp,logP,fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape_eQTL,guide='none')+
          scale_size_manual(values=size_eQTL,guide='none')+
          scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(region_start,region_end))+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab(paste0(chr,' (Mb)'))+
          ylab(bquote(-log[10]*'(P)'))+
          labs(title = eQTL_title)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(eQTL_title, "\n") + 1)))
          ) + 
          scale_y_continuous(limits = c(-0.1, max(eQTL$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0, 0.05)))
p_eQTL =p_eQTL + geom_vline(xintercept = ref_pos, linetype=2) +
          geom_point(aes(x=ref_pos, y=eQTL[eQTL$bp ==ref_pos,]$logP), shape=23, size=3, fill="purple")
 
###############epithelium sceQTL plot
###
shape_sceQTL = ifelse(sceQTL$chrbp == ref, 23, 21)
names(shape_sceQTL) = sceQTL$chrbp
###
size_sceQTL = ifelse(sceQTL$chrbp == ref, 3, 2)
names(size_sceQTL) = sceQTL$chrbp
###
sceQTL$label = ifelse(sceQTL$chrbp == ref, sceQTL$variant_id, '')
###
sceQTL_title = paste0(targetGene," Epithelium sc-eQTL")
p_sceQTL = ggplot(sceQTL,aes(x=bp,y=logP))+
          geom_point(aes(fill=chrbp,size=chrbp,shape=chrbp),alpha=0.8)+
          geom_point(data=sceQTL[sceQTL$label!='',],aes(x=bp,logP,fill=chrbp,size=chrbp,shape=chrbp))+
          scale_fill_manual(values=color,guide='none')+
          scale_shape_manual(values=shape_sceQTL,guide='none')+
          scale_size_manual(values=size_sceQTL,guide='none')+
          scale_x_continuous(labels=function(x){sprintf('%.1f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(region_start,region_end))+
          ggrepel::geom_text_repel(aes(label=label), max.overlaps = Inf)+
          xlab(paste0(chr,' (Mb)'))+
          ylab(bquote(-log[10]*'(P)'))+
          labs(title = sceQTL_title)+
          theme_classic()+
          theme(plot.margin=unit(c(0.5, 1, 0.5, 0.5), "lines"),
          plot.title = element_text(
              hjust = 1e-2, 
              margin = margin(b = -12 * (stringr::str_count(sceQTL_title, "\n") + 1)))
          ) + 
          scale_y_continuous(limits = c(-0.1, max(sceQTL$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0, 0.05))) 
p_sceQTL =p_sceQTL + geom_vline(xintercept = ref_pos, linetype=2) +
          geom_point(aes(x=ref_pos, y=sceQTL[sceQTL$bp ==ref_pos,]$logP), shape=23, size=3, fill="purple")

###############together
p_together <- p_caQTL + scale_y_continuous(limits = c(-0.1, max(caQTL$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0.05, 0.15))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
              p_sceQTL + scale_y_continuous(limits = c(-0.1, max(sceQTL$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0.05, 0.15))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
              p_eQTL + scale_y_continuous(limits = c(-0.1, max(eQTL$logP, na.rm = TRUE) * 1.05), expand = expansion(mult = c(0.05, 0.15))) + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
              p_Peak +
              p_gene +
              plot_layout(nrow = 5, heights = c(0.7, 0.7, 0.7, 0.2, 0.4))
###LD r2 legend
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
ggplot2::ggsave(p_together_withlegend,filename=paste0(outdir,"/",targetPeak,"~",targetGene,".caQTL_bulkeQTL_sceQTL_GWAS.locuszoomplot.pdf"),width=150, height=200, units="mm") ##Figure4F

