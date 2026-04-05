##########################################
############FigureS6I-J###################
##########################################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATAC_for_review/FigureS6I-J/output")
library(ggplot2)
results=read.delim("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR",header=T)
results_sig=read.delim("/data1/gy/ATACseq_RWAS/ASCA_STITCH/ASCA_result/region_250bp/95sample_stratas_output_files_250bp_combined_withFDR_FDR0.1",header=T)

# Allele fraction for significantly imbalanced SNPs  
sig.df <- data.frame(AF = results$ALL.AF[results$C0.BBINOM.FDR < 0.1])  

p1 <- ggplot(sig.df, aes(x = AF)) +  
  geom_histogram(binwidth = 0.025, fill = "#0073C2B2", color = "#e9ecef", alpha = 0.7) +  
  labs(  
    title = "Histogram of AF for significantly imbalanced SNPs",  
    x = "Allele fraction",  
    y = "Frequency"  
  ) +  
  theme_minimal() +  
  theme(  
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 10),  
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")  
  ) +  
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  

# Allele fraction for balanced SNPs  
notsig.df <- data.frame(AF = results$ALL.AF[results$C0.BBINOM.FDR >= 0.1 & results$N.READS > 0])

p2 <- ggplot(notsig.df, aes(x = AF)) +  
  geom_histogram(binwidth = 0.025, fill = "#0073C2B2", color = "#e9ecef", alpha = 0.7) +  
  labs(  
    title = "Histogram of AF for balanced SNPs",  
    x = "Allele fraction",  
    y = "Frequency"  
  ) +  
  theme_minimal() +  
  theme(  
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 10),  
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")  
  ) +  
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  

#
ggsave(p1, file = "ASCA_sig_snp_AF_distribution.barplot.pdf", width = 6, height = 6)  ##FigureS6I
ggsave(p2, file = "ASCA_nosig_snp_AF_distribution.barplot.pdf", width = 6, height = 6)  ##FigureS6J
