##############################
###########Figure3F###########
##############################
#####ref: PMID 39848247
##################
###1.input file###
##################
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input
#bed file
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(openxlsx)
##############################
##
caSNP_inpeak=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.bed")
asSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_ASCA_FDR0.1.SNP.bed")
caSNP_inpeak_nosig=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_caQTL_nosig.SNP.in_corresponding_peak.bed")
asSNP_nosig=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/motifbreakR/input/95sample_ASCA_nosig.SNP.bed")
##
fwrite(caSNP_inpeak[,1:4],"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.bed",col.names=F,row.names=F,quote=F,sep="\t")
fwrite(asSNP[,1:4],"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/95sample_ASCA_FDR0.1.SNP.bed",col.names=F,row.names=F,quote=F,sep="\t")
fwrite(caSNP_inpeak_nosig[,1:4],"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/95sample_caQTL_nosig.SNP.in_corresponding_peak.bed",col.names=F,row.names=F,quote=F,sep="\t")
fwrite(asSNP_nosig[,1:4],"/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/95sample_ASCA_nosig.SNP.bed",col.names=F,row.names=F,quote=F,sep="\t")
##############################bg
SNPs=fread("/data1/gy/ATACseq_RWAS/genotype/from_bam_STITCH/STITCH_QC_foranalysis/plink/95sample.stitch.chr1_22.INFO0.5.maf0.05_geno0.95_HWE0.000001.convert_GT.withSNPID.recode.for_analysis.withchr.bim")
SNPs=SNPs[,c(1,4,2)]
SNPs$start=SNPs$V4-1
SNPs_bed=SNPs[,c(1,4,2,3)]
#
set.seed(123) 
# 
target_size <- 98492
#
if (nrow(SNPs_bed) > target_size) {
  #
  sampled_indices <- sample(nrow(SNPs_bed), target_size, replace = FALSE)
  background_snps <- SNPs_bed[sampled_indices, ]
} else {
  background_snps <- SNPs_bed
}

#
setorder(background_snps, V1, start)
background_snps$V1 <- paste0("chr", background_snps$V1)

#
fwrite(background_snps, "/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/SNPs_background.bed", 
       sep = "\t", col.names = FALSE, quote = FALSE)


#############
###2.homer###
#############
mkdir -p /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/homer_output
###########motif enrichment
cd /Public/gaoyun/software/homer 
#perl configureHomer.pl -install hg19
export PATH=$PATH:/Public/gaoyun/software/homer/bin/  
cd /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/homer_output
##
mkdir -p caSNP_in_corresponding_peak
findMotifsGenome.pl \
  /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/95sample_caQTL_FDR0.1.SNP.in_corresponding_peak.bed \
  hg19 \
  caSNP_in_corresponding_peak \
  -bg /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/SNPs_background.bed \
  -size 20 \
  -p 10
##
mkdir -p asSNP
findMotifsGenome.pl \
  /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/95sample_ASCA_FDR0.1.SNP.bed \
  hg19 \
  asSNP \
  -bg /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/SNPs_background.bed \
  -size 20 \
  -p 10
##
mkdir -p caSNP_nosig_in_corresponding_peak
findMotifsGenome.pl \
  /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/95sample_caQTL_nosig.SNP.in_corresponding_peak.bed \
  hg19 \
  caSNP_nosig_in_corresponding_peak \
  -bg /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/SNPs_background.bed \
  -size 20 \
  -p 10
##
mkdir -p asSNP_nosig
findMotifsGenome.pl \
  /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/95sample_ASCA_nosig.SNP.bed \
  hg19 \
  asSNP_nosig \
  -bg /data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/input/SNPs_background.bed \
  -size 20 \
  -p 10

######################knownmotif
#
cat caSNP_in_corresponding_peak/knownResults/*.motif > caSNP_in_corresponding_peak/knownMotifs.motifs
cat asSNP/knownResults/*.motif > asSNP/knownMotifs.motifs
cat caSNP_nosig_in_corresponding_peak/knownResults/*.motif > caSNP_nosig_in_corresponding_peak/knownMotifs.motifs
cat asSNP_nosig/knownResults/*.motif > asSNP_nosig/knownMotifs.motifs

#####################################
###3.calculate Motif Content Score###
#####################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(openxlsx)
########################################1.caSNP_in_corresponding_peak
caSNP_in_corresponding_peak=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/homer_output/caSNP_in_corresponding_peak/knownResults.txt")
caSNP_in_corresponding_peak_q0.05=subset(caSNP_in_corresponding_peak,caSNP_in_corresponding_peak$'q-value (Benjamini)'<0.05)
###
df <- caSNP_in_corresponding_peak_q0.05
total_elements=6039   ##elements filtered and retained by HOMER's analysis
#
df[, target_num := as.numeric(gsub("%", "", `% of Target Sequences with Motif`))]
df[, bg_num     := as.numeric(gsub("%", "", `% of Background Sequences with Motif`))]

#
df[, diff_val := target_num - bg_num]

# 
sum_diff <- sum(df$diff_val, na.rm = TRUE)
n_motifs <- nrow(df)

final_tf_content_score.caSNP_in_corresponding_peak <-  (sum_diff * n_motifs) / sqrt(total_elements)
final_tf_content_score.caSNP_in_corresponding_peak
#[1] 66.59906

########################################2.caSNP_nosig_in_corresponding_peak
caSNP_nosig_in_corresponding_peak=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/homer_output/caSNP_nosig_in_corresponding_peak/knownResults.txt")
caSNP_nosig_in_corresponding_peak_q0.05=subset(caSNP_nosig_in_corresponding_peak,caSNP_nosig_in_corresponding_peak$'q-value (Benjamini)'<0.05)
###
df <- caSNP_nosig_in_corresponding_peak_q0.05
total_elements=92453   ##elements filtered and retained by HOMER's analysis
#
df[, target_num := as.numeric(gsub("%", "", `% of Target Sequences with Motif`))]
df[, bg_num     := as.numeric(gsub("%", "", `% of Background Sequences with Motif`))]

#
df[, diff_val := target_num - bg_num]

#
sum_diff <- sum(df$diff_val, na.rm = TRUE)
n_motifs <- nrow(df) 

final_tf_content_score.caSNP_nosig_in_corresponding_peak <-  (sum_diff * n_motifs) / sqrt(total_elements)
final_tf_content_score.caSNP_nosig_in_corresponding_peak
#[1] 9.578346

########################################3.asSNP
asSNP=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/homer_output/asSNP/knownResults.txt")
asSNP_q0.05=subset(asSNP,asSNP$'q-value (Benjamini)'<0.05)
###
df <- asSNP_q0.05
total_elements=10470  ##elements filtered and retained by HOMER's analysis
#
df[, target_num := as.numeric(gsub("%", "", `% of Target Sequences with Motif`))]
df[, bg_num     := as.numeric(gsub("%", "", `% of Background Sequences with Motif`))]

#
df[, diff_val := target_num - bg_num]

#
sum_diff <- sum(df$diff_val, na.rm = TRUE) 
n_motifs <- nrow(df)  

final_tf_content_score.asSNP <-  (sum_diff * n_motifs) / sqrt(total_elements)
final_tf_content_score.asSNP
#[1] 44.90446

########################################4.asSNP_nosig
asSNP_nosig=fread("/data1/gy/ATACseq_RWAS/caQTL_vs_ASCA_STITCH/TF_content_score/homer_output/asSNP_nosig/knownResults.txt")
asSNP_nosig_q0.05=subset(asSNP_nosig,asSNP_nosig$'q-value (Benjamini)'<0.05)
###
df <- asSNP_nosig_q0.05
total_elements=88012 ##elements filtered and retained by HOMER's analysis
# 
df[, target_num := as.numeric(gsub("%", "", `% of Target Sequences with Motif`))]
df[, bg_num     := as.numeric(gsub("%", "", `% of Background Sequences with Motif`))]

#
df[, diff_val := target_num - bg_num]

#
sum_diff <- sum(df$diff_val, na.rm = TRUE)
n_motifs <- nrow(df)

final_tf_content_score.asSNP_nosig <-  (sum_diff * n_motifs) / sqrt(total_elements)
final_tf_content_score.asSNP_nosig
#[1] 8.761371


###################################
###4.Motif Content Score barplot###
###################################
library(ggplot2)
library(scales)

data <- c(
  log2(final_tf_content_score.caSNP_nosig_in_corresponding_peak),
  log2(final_tf_content_score.caSNP_in_corresponding_peak),
  log2(final_tf_content_score.asSNP_nosig),
  log2(final_tf_content_score.asSNP)
)
df <- data.frame(
  group = factor(c("noncaSNP","caSNP","nonasSNP","asSNP"),
                 levels=c("noncaSNP","caSNP","nonasSNP","asSNP")),
  score  = data,
  x     = c(1, 2, 4, 5) 
)

#
p <- ggplot(df, aes(x=x, y=score, fill=group)) +
  #
  geom_col(width=0.8, alpha=0.75) +  
  
  #
  geom_text(aes(label = round(score, 2)),
            vjust = -0.5, size = 4) +
  
  #
  scale_x_continuous(breaks = df$x, 
                     labels = df$group,
                     expand = expansion(mult=c(0.05, 0.05))) +
  
  #
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    limits = c(0, NA) 
  ) +
  
  #
  scale_fill_manual(values=c("#BEBEBE", "#094D92", "#7A7A7A", "#238444")) +
  
  #
  labs(x="", y="log2(TF content score)", title="", fill=NULL) +
  
  #
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),            
    axis.line = element_line(color="black"), 
    axis.ticks = element_line(color="black"),
    axis.ticks.length = unit(3, "pt"),
    axis.text.x = element_text(color="black", size=12), 
    axis.text.y = element_text(color="black"),
    legend.position = "none" 
  )

ggsave("TF_motif_content_score.barplot.pdf",p,width=6,height=7.5) ##Figure3F
