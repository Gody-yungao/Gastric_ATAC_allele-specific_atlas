#################################
###########Figure5F-I############
#################################
###################################
###########1.caQTL plot############
###################################
/Public/gaoyun/software/R-4.2.0/bin/R
rm(list = ls())
library(data.table)
##
chr = 4
bp = 48073300
peak = "chr4:48072992-48073492"
##
caQTL=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.txt",header=T)
caQTL=subset(caQTL,(pheno_id == peak) & (var_start == bp))
##
exp=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/exp/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.bed.gz")
exp=subset(exp,Geneid %in% peak)
##
dosage=fread("/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.phased.geno.dosage")
dosage=subset(dosage,CHROM == paste0("chr",chr) & POS ==bp)

##
colnames(exp)[4]="ID"
input=as.data.frame(rbind(dosage[,c(3,6:100)],exp[,c(4,7:101)]))
rownames(input)=input$ID
input=input[,-1]
input=as.data.frame(t(input))

############################################Plot
library(ggplot2)  
##
SNP = colnames(input)[1]  
OCR = colnames(input)[2]  
outdir = "/data1/gy/ATAC_for_review/Figure5F-I/output/"

##########
REF <- dosage$REF[1]  
ALT <- dosage$ALT[1]   
beta = caQTL$beta[1]  
P_value = caQTL$p_nominal[1]  
#
input[, SNP] <- ifelse(input[, SNP] == 0, paste0(REF, REF),   
                            ifelse(input[, SNP] == 1, paste0(REF, ALT),   
                                   ifelse(input[, SNP] == 2, paste0(ALT, ALT), NA)))  
#
count_00 <- sum(input[, SNP] == paste0(REF, REF))  
count_01 <- sum(input[, SNP] == paste0(REF, ALT)) 
count_11 <- sum(input[, SNP] == paste0(ALT, ALT))  
count_00
#[1] 8
count_01
#[1] 37
count_11
#[1] 50
#
if (count_00 < count_11) {  
  input[, SNP] <- ifelse(input[, SNP] == paste0(REF, ALT), paste0(ALT,REF), input[, SNP])  
}  
#
if (count_00 < count_11) { 
  input[, SNP] <- factor(input[, SNP], levels = c(paste0(ALT, ALT), paste0(ALT, REF), paste0(REF, REF)))  
}else{
  input[, SNP] <- factor(input[, SNP], levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
}

#
if (count_00 < count_11) { 
  beta_new = -beta  
}else{
  beta_new = beta
}

##
p <- ggplot(data = input, aes(x = input[, SNP], y = input[, OCR])) +
    geom_boxplot(
        aes(color = input[, SNP]), fill = "white",
        width = 0.6, alpha = 1,
        outlier.shape = NA
    ) +
    geom_jitter(aes(color = input[, SNP]), width = 0.2, size = 1.3) +
    scale_color_manual(values = c("#2B8CBE", "#525252", "#CB181D"), name = "") +
    labs(  
    x = "rs875179 genotype",  
    y = "Normalized accessibility",  
    title = OCR,
    subtitle = bquote(  
      atop(
        "caQTL",
        italic(beta) ~ "= " ~ .(beta_new) * ", " ~ italic(P) ~ "= " ~ .(P_value)
       )  
     )  
    )
p <- p +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "black", hjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        legend.position = "none"
    )
# 
ggsave(filename = paste0(outdir, "rs875179", "~", OCR, ".caQTL.scatterplot.pdf"),
         plot = p, width = 3, height = 3.5)  ##Figure5F
q()

###################################
###########2.ASCA plot############
###################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)  
library(tidyr)  
library(ggplot2) 
ASCA=fread("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR")

##
chr = 4
bp = 48073300
peak = "chr4:48072992-48073492"
data=subset(ASCA,CHR==paste0("chr",chr) & POS==bp)

###################################REF
data_long_REF <- data %>%  
  select(CHR, POS, RSID, IND.C0.COUNT.REF, ALL.AF, ALL.BBINOM.P, C0.BBINOM.FDR) %>%  
  mutate(REF = sub(".*:(.*):(.*)", "\\1", RSID),
         ALT = sub(".*:(.*):(.*)", "\\2", RSID)) %>%
  mutate(Type = "REF",Count=IND.C0.COUNT.REF) %>%  
  separate_rows(Count, sep = ",") %>%
  mutate(Count = as.numeric(Count)) %>%  
  #
  mutate(BASE = REF) %>%  
  #
  select(-REF, -ALT,-IND.C0.COUNT.REF)%>%  
  #
  group_by(RSID) %>%  
  mutate(Individual = paste0("IND", row_number()))
###################################ALT
data_long_ALT <- data %>%  
  select(CHR, POS, RSID, IND.C0.COUNT.ALT, ALL.AF, ALL.BBINOM.P, C0.BBINOM.FDR) %>%  
  mutate(REF = sub(".*:(.*):(.*)", "\\1", RSID),
         ALT = sub(".*:(.*):(.*)", "\\2", RSID)) %>%
  mutate(Type = "ALT",Count=IND.C0.COUNT.ALT) %>%  
  separate_rows(Count, sep = ",") %>%
  mutate(Count = as.numeric(Count)) %>%  
  #
  mutate(BASE = ALT) %>%  
  #
  select(-REF, -ALT,-IND.C0.COUNT.ALT)%>%  
  #
  group_by(RSID) %>%  
  mutate(Individual = paste0("IND", row_number()))
###combine
data_long=rbind(data_long_REF,data_long_ALT)

#################
library(dplyr)
library(tidyr)
##
data_long <- data %>%
  select(CHR, POS, RSID, IND.C0, IND.C0.COUNT.REF, IND.C0.COUNT.ALT, ALL.AF, ALL.BBINOM.P, C0.BBINOM.FDR) %>%  
  mutate(REF = sub(".*:(.*):(.*)", "\\1", RSID),
         ALT = sub(".*:(.*):(.*)", "\\2", RSID)) %>%
  mutate(
    IND.C0 = strsplit(IND.C0, ","),
    IND.C0.COUNT.REF = strsplit(IND.C0.COUNT.REF, ","),
    IND.C0.COUNT.ALT = strsplit(IND.C0.COUNT.ALT, ",")
  ) %>%
  unnest(cols = c(IND.C0, IND.C0.COUNT.REF, IND.C0.COUNT.ALT)) %>%
  mutate(
    IND.C0 = as.numeric(IND.C0),
    IND.C0.COUNT.REF = as.numeric(IND.C0.COUNT.REF),
    IND.C0.COUNT.ALT = as.numeric(IND.C0.COUNT.ALT)
  )

##haplotype count
res <- data_long %>%
  group_by(IND.C0) %>%
  summarise(
    CHR  = first(CHR),
    POS  = first(POS),
    RSID = first(RSID),
    IND.C0.COUNT.REF = sum(IND.C0.COUNT.REF, na.rm = TRUE),
    IND.C0.COUNT.ALT = sum(IND.C0.COUNT.ALT, na.rm = TRUE),
    ALL.AF = first(ALL.AF),
    ALL.BBINOM.P = first(ALL.BBINOM.P),
    C0.BBINOM.FDR = first(C0.BBINOM.FDR),
    REF = first(REF),
    ALT = first(ALT),
    .groups = "drop"
  )
dim(res)
#[1] 37  11

###
res_ref=res[,c(1:5,7:10)]
res_alt=res[,c(1:4,6:9,11)]
colnames(res_ref)[5]="IND.C0.COUNT"
colnames(res_alt)[5]="IND.C0.COUNT"
colnames(res_ref)[9]="BASE"
colnames(res_alt)[9]="BASE"
res=rbind(res_ref,res_alt)

#
af_value <- unique(res$ALL.AF)[1]   
alt_ref_ratio_log2 = log2(af_value/(1-af_value))
p_value <- unique(res$ALL.BBINOM.P)[1]   
ref <- res$BASE[1]  
alt <- res$BASE[nrow(res)]  

#
res$BASE <- factor(res$BASE, levels = c(alt, ref))  

#########################Plot
outdir = "/data1/gy/ATAC_for_review/Figure5F-I/output/"
p <- ggplot(res, aes(x = BASE, y = IND.C0.COUNT, group = IND.C0)) +  
  geom_point(aes(color = BASE), size = 1) +
  geom_line(color = "black", linetype = "solid") +
  scale_color_manual(values = c("#2B8CBE", "#CB181D"), name = "") +
  labs(  
    x = "rs875179 allele",  
    y = "Haplotype ATAC read counts",  
    title = peak,
    subtitle = bquote(  
      atop(
        "ASCA",
        "log2(G/C) = " ~ .(-alt_ref_ratio_log2) * ", " ~ italic(P) ~ "= " ~ .(p_value)
       )  
     )  
    )+
  scale_y_continuous(
    limits = c(0, NA)
  ) 
p <- p +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "black", hjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        legend.position = "none"
    )

#
ggsave(p, file = paste0(outdir, "rs875179", "~", peak, ".ASCA.scatterplot.pdf"), width = 3, height = 3.5)  ##Figure5G
q()

###########################################
###########3.262 bulk eQTL plot############
###########################################
/Public/gaoyun/software/R-4.2.0/bin/R
rm(list = ls())
library(data.table)
##
chr = 4
pos = 48073300
symbol1 = "NFXL1"
geneid1 = "ENSG00000170448"

##
eQTL=fread("/data1/gy/262_eQTL_final_result/262sample_eqtl.cis_qtl_pairs.maf0.01.with_FDR.with_symbol.all.chr.txt",header=T)
library(dplyr)
eQTL <- eQTL %>% 
  mutate(bp = as.integer(sub("^[^:]+:([0-9]+):.*$", "\\1", variant_id)))
eQTL=subset(eQTL,(eQTL$gene_name == symbol1 |eQTL$gene_name == symbol2) & eQTL$bp == pos)
##
exp=fread("/data1/gy/262_eQTL_final_result/262RNAseq_FID_exp_qnorm.txt")
exp$FID=substr(exp$FID,1,15)
exp=subset(exp,FID %in% geneid1)
##dosage
dosage=fread("/data1/gy/262_eQTL_final_result/geno/262geno.all.chr.Rsq03.MAF01.G05.M05.Hwe6e.geno.dosage")
dosage=subset(dosage,dosage$CHROM == chr & dosage$POS ==pos)

##
colnames(exp)[1]="ID"
input=as.data.frame(rbind(dosage[,c(3,6:ncol(dosage)), with = FALSE],exp))
rownames(input)=input$ID
input=as.data.frame(t(input[,-1]))

############################################Plot
library(ggplot2)  
##
SNP = colnames(input)[1]  
outdir = "/data1/gy/ATAC_for_review/Figure5F-I/output/"

##########
REF <- dosage$REF[1]  
ALT <- dosage$ALT[1]   
beta1 = eQTL$slope[2]  
P_value1 = eQTL$pval_nominal[2]  
beta2 = eQTL$slope[1]  
P_value2 = eQTL$pval_nominal[1]  
#
input[, SNP] <- ifelse(input[, SNP] == 0, paste0(REF, REF),   
                            ifelse(input[, SNP] == 1, paste0(REF, ALT),   
                                   ifelse(input[, SNP] == 2, paste0(ALT, ALT), NA)))  
#
count_00 <- sum(input[, SNP] == paste0(REF, REF))  
count_01 <- sum(input[, SNP] == paste0(REF, ALT)) 
count_11 <- sum(input[, SNP] == paste0(ALT, ALT))
count_00
#[1] 122
count_01
#[1] 116
count_11
#[1] 24
#
if (count_00 < count_11) {  
  input[, SNP] <- ifelse(input[, SNP] == paste0(REF, ALT), paste0(ALT,REF), input[, SNP])  
}  
#
if (count_00 < count_11) { 
  input[, SNP] <- factor(input[, SNP], levels = c(paste0(ALT, ALT), paste0(ALT, REF), paste0(REF, REF)))  
}else{
  input[, SNP] <- factor(input[, SNP], levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
}

#
if (count_00 < count_11) { 
  beta_new1 = -beta1  
}else{
  beta_new1 = beta1
}

##
##############NFXL1
p1 <- ggplot(data = input, aes(x = input[, SNP], y = input[, geneid1])) +
    geom_boxplot(
        aes(color = input[, SNP]), fill = "white",
        width = 0.6, alpha = 1,
        outlier.shape = NA
    ) +
    geom_jitter(aes(color = input[, SNP]), width = 0.2, size = 1.3) +
    scale_color_manual(values = c("#4EB3D3", "#737373", "#EF3B2C"), name = "") +
    labs(  
    x = "rs875179 genotype",  
    y = "Normalized expression",  
    title = symbol1,
    subtitle = bquote(  
      atop( 
        "bulk eQTL",
        italic(beta) ~ "= " ~ .(beta_new1) * ", " ~ italic(P) ~ "= " ~ .(P_value1)
       )  
     )  
    )
p1 <- p1 +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "black", hjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        legend.position = "none"
    )

#
ggsave(filename = paste0(outdir, "rs875179", "~", symbol1, ".262bulk_eQTL.scatterplot.pdf"),   
         plot = p1, width = 3, height = 3.5)  ##Figure5H
q()

####################################################
###########4.203 epithelium sc-eQTL plot############
####################################################
/Public/gaoyun/software/R-4.2.0/bin/R
rm(list = ls())
library(data.table)
library(dplyr)
##
chr = 4
bp = 48073300
symbol = "NIPAL1"
geneid = "ENSG00000163293"
##
sceQTL=fread("/data1/gy/sceQTL/eQTL_CHN100K/epi/epi/epi_sceQTL_result.allpairs.withsymbol.withFDR.FDR0.05.txt",header=T)
variant_pattern <- paste0("^chr", chr, ":", bp, "(:|$)")
sceQTL <- subset(sceQTL,  
                          gene_id == geneid &  
                          grepl(variant_pattern, variant_id))  
##
exp=fread("/data/gc_sceqtl/sceqtl/seuratall/allfilter_orig/epi/epi/round1_residual_HP_df.xls")
exp=as.data.frame(exp[,c("sampleid",geneid),with=F])
dim(exp)
#[1] 203   2

#
dosage=fread("/data1/gy/sceQTL/genotype/epi/epi_CHN100K_filter05.dosage.txt")
dosage=subset(dosage,dosage$'#CHROM' == chr & dosage$POS ==bp)
#
REF <- dosage$REF[1]  
ALT <- dosage$ALT[1] 
#
dosage=as.data.frame(dosage[,c(3,10:ncol(dosage)), with = FALSE])
rownames(dosage)=dosage$ID
dosage=dosage[,-1]
dosage=as.data.frame(t(dosage))
dosage <- dosage %>%   
  as.data.frame() %>%   
  tibble::rownames_to_column("sampleid") 
dim(dosage)
#[1] 203   2

##
input=merge(dosage,exp)
dim(input)
#[1] 203   3
input=input[,-1]
input=na.omit(input)
dim(input)
#[1] 203   2

############################################Plot
library(ggplot2)  
##
SNP = colnames(input)[1]  
geneid = colnames(input)[2]  
outdir = "/data1/gy/ATAC_for_review/Figure5F-I/output/"

##########
beta = sceQTL$slope[1]  
P_value = sceQTL$pval_nominal[1]  
#
input[, SNP] <- ifelse(input[, SNP] == 0, paste0(REF, REF),   
                            ifelse(input[, SNP] == 1, paste0(REF, ALT),   
                                   ifelse(input[, SNP] == 2, paste0(ALT, ALT), NA)))  
#
count_00 <- sum(input[, SNP] == paste0(REF, REF))  
count_01 <- sum(input[, SNP] == paste0(REF, ALT)) 
count_11 <- sum(input[, SNP] == paste0(ALT, ALT))
count_00
#[1] 101
count_01
#[1] 82
count_11
#[1] 20
#
if (count_00 < count_11) {  
  input[, SNP] <- ifelse(input[, SNP] == paste0(REF, ALT), paste0(ALT,REF), input[, SNP])  
}  
#
if (count_00 < count_11) { 
  input[, SNP] <- factor(input[, SNP], levels = c(paste0(ALT, ALT), paste0(ALT, REF), paste0(REF, REF)))  
}else{
  input[, SNP] <- factor(input[, SNP], levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
}

#
if (count_00 < count_11) { 
  beta_new = -beta  
}else{
  beta_new = beta
}

##
##############NIPAL1
p <- ggplot(data = input, aes(x = input[, SNP], y = input[, geneid])) +
    geom_boxplot(
        aes(color = input[, SNP]), fill = "white",
        width = 0.6, alpha = 1,
        outlier.shape = NA
    ) +
    geom_jitter(aes(color = input[, SNP]), width = 0.2, size = 1.3) +
    scale_color_manual(values = c("#4EB3D3", "#737373", "#EF3B2C"), name = "") +
    labs(  
    x = "rs875179 genotype",  
    y = "Normalized expression",  
    title = symbol,
    subtitle = bquote(  
      atop(
        "Epithelium sc-eQTL",
        italic(beta) ~ "= " ~ .(beta_new) * ", " ~ italic(P) ~ "= " ~ .(P_value)
       )  
     )  
    )
p <- p +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, colour = "black", fill = NA),
        plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "black", hjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),
        axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
        axis.ticks.length.x = unit(-0.1, "cm"),
        axis.ticks.length.y = unit(-0.1, "cm"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.line = element_blank(),
        legend.position = "none"
    )

#
ggsave(filename = paste0(outdir, "rs875179", "~", symbol, ".epi_sceQTL.scatterplot.pdf"),   
         plot = p, width = 3, height = 3.5)  ##Figure5I
q()
