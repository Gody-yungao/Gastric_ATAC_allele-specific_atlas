############################
#########FigureS10B#########
############################
########################################################################
########1.RWAS conditioning test of OCR chr4:48072992-48073492##########
########################################################################
conda activate fusion
mkdir /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/RWAS_postprocess
cd /data1/gy/ATACseq_RWAS/GC_candidate_region_final_STITCH/chr4:48072992-48073492/RWAS_postprocess
##
chr=4
infile=/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter/sig/RWAS.fdr0.1.sig.result.txt  
outfile=RWAS.fdr0.1.sig.chr4:48072992-48073492.result.txt  

##
head -n1 "$infile" > "$outfile"  
grep 'chr4:48072992-48073492' "$infile" >> "$outfile"  
model=$(tail -n+2 "$outfile" | cut -f16)  

##test
Rscript /data1/gy/software/cwas-master/fusion/FUSION.post_process.R \
--input RWAS.fdr0.1.sig.chr4:48072992-48073492.result.txt \
--out chr4:48072992-48073492 \
--sumstats /data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/input/gwas/GWAS_QCandFilter.new.with_SNPid_RWAS.without_ambiguous.for_RWAS.QC.sumstats.gz \
--ref_ld_chr /data1/gy/ATACseq_RWAS/RWAS_STITCH/ldref/per_chr/1000G.EAS. \
--plot TRUE \
--plot_eqtl TRUE \
--plot_corr TRUE \
--eqtl_model $model \
--save_loci TRUE \
--report TRUE \
--plot_individual TRUE \
--chr $chr

################################################################################################################
########2.GWAS results conditioning on best RWAS model genetic predictor of OCR chr4:48072992-48073492##########
################################################################################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
rwas <- fread("RWAS.fdr0.1.sig.chr4:48072992-48073492.result.txt")

##chromosome region
upper <- rwas$P1
lower <- rwas$P0
min.limit <- lower - 200000
max.limit <- upper + 200000

##GWAS results of conditioning test 
library(dplyr)
library(tidyr)
dat <- read.csv("chr4:48072992-48073492.loc_1.cond.csv") %>%
    drop_na(GWAS.LOGP) %>%
    filter(POS >= min.limit) %>%
    filter(POS <= max.limit)

##
rawgwas <- dat %>%
    select(POS, logpv=GWAS.LOGP) %>%
    mutate(note = "Main GWAS")
condgwas <- dat %>%
    select(POS, logpv= GWAS_cond.LOGP) %>%
    mutate(note = "Conditioned on\npredictor")
dat1 <- rbind(rawgwas,condgwas)

##plot
library(ggplot2)
GWASp0 <- ggplot(data = dat1, aes(x = POS/1e6, y = logpv)) +
    geom_point(
        aes(color = note), pch = 16,
        size = 1.2, alpha = 1
    ) +
    scale_color_manual(
        values = rev(c("grey", "#F08E59")),
        name = ""
    ) +
    xlab("chr4 (Mb)") +
    ylab(expression("GWAS"~"-log"[10]*"("*italic("P")*")")) +
    scale_x_continuous(expand = c(0.005, 0.005)) +
    scale_y_continuous(expand = c(0.02, 0.02))
GWASp1 <- GWASp0 +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),

    #
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),

    axis.text.x = element_text(size = 8, color = "black", margin = margin(t = 6)),
    axis.text.y = element_text(size = 8, color = "black", margin = margin(r = 6)),

    #
    axis.ticks.length.x = unit(0.15, "cm"),
    axis.ticks.length.y = unit(0.15, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.3),

    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),

    legend.position = c(0.85, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.5
    )
  ) +
  #
  coord_cartesian(clip = "off")
ggsave(GWASp1,filename="chr4:48072992-48073492.loc_1.cond.scatterplot.pdf",
       width=4,height=2.5) ##FigureS10B
q()
