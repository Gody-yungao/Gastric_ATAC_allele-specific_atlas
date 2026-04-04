#######################
#######FigureS1B#######
#######################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATAC_for_review/FigureS1B/output")
library(openxlsx)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(patchwork)
sample102=read.xlsx("/data1/gy/ATACseq_RWAS/ATACseq/QC/102sample_QC.xlsx")
sample95=subset(sample102,PBC1 >= 0.7 & PBC2 >= 1 & NRF >= 0.7 & TSS_enrichment_score >=6 & FRiP >=0.2)
dim(sample95)
#[1] 95  7
#
sample102$group="ATACseq(n=102)"
sample95$group="ATACseq(n=95)"
#col <-c("#5CB85C","#337AB7")

#####quality anno
library(dplyr)
sample102 <- sample102 %>%
  mutate(
    QC = if_else(SampleID %in% sample95$SampleID,
                 "high_quality", "low_quality")
  )

#########plot
col_bfQC="#337AB7"
#####################################1.TSS enrichment
p1 <- ggplot(sample102, aes(x = group, y = TSS_enrichment_score, fill = NULL, color=group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(
    aes(shape = QC),
    width = 0.2, size = 3, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  scale_y_continuous(limits = c(0, max(sample102$TSS_enrichment_score) + 0.5)) +
  geom_hline(yintercept = 6, linetype = "dashed", color = "grey", size = 0.5) +
  scale_color_manual(values = col_bfQC, guide = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(y = "TSS enrichment") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "right"  )


#####################################2.FRiP
p2 <- ggplot(sample102, aes(x = group, y = FRiP, fill = NULL, color=group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(
    aes(shape = QC),
    width = 0.2, size = 3, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey", size = 0.5) +
  scale_color_manual(values = col_bfQC, guide = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(y = "FRiP") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "right"  )

#####################################3.NRF
p3 <- ggplot(sample102, aes(x = group, y = NRF, fill = NULL, color=group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(
    aes(shape = QC),
    width = 0.2, size = 3, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  scale_y_continuous(limits = c(0,max(sample102$NRF)+0.05)) +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "grey", size = 0.5) +
  scale_color_manual(values = col_bfQC, guide = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(y = "NRF") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "right"  )

#####################################4.PBC1
p4 <- ggplot(sample102, aes(x = group, y = PBC1, fill = NULL, color=group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(
    aes(shape = QC),
    width = 0.2, size = 3, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  scale_y_continuous(limits = c(0,max(sample102$PBC1)+0.05)) +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "grey", size = 0.5) +
  scale_color_manual(values = col_bfQC, guide = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),  
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 12)) + 
  labs(y = "PBC1") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position = "right"  ) 

#####################################5.PBC2
p5 <- ggplot(sample102, aes(x = group, y = PBC2, fill = NULL, color=group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  geom_jitter(
    aes(shape = QC),  
    width = 0.2, size = 3, show.legend = TRUE
  ) +
  scale_shape_manual(
    name = "QC",
    values = c("high_quality" = 20, "low_quality" = 17),
    labels = c("High quality", "Low quality")
  ) +
  scale_y_continuous(limits = c(0,max(sample102$PBC2)+0.05)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey", size = 0.5) + 
  scale_color_manual(values = col_bfQC, guide = "none") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 0.5), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 12)) + 
  labs(y = "PBC2") +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"), 
        legend.position = "right"  ) 
		
#
fix_x_blank <- theme(axis.text.x = element_blank(), axis.title.x = element_blank())
combined_plot <- plot_grid(
  p1 + fix_x_blank,
  p2 + fix_x_blank,
  p4 + fix_x_blank,
  p5 + fix_x_blank,
  p3+ fix_x_blank,
  ncol = 1,    
  rel_heights = rep(1, 5),
  align = "v"
)
ggsave("102sample_QC.boxplot.pdf",combined_plot,width=4,height=20) ##FigureS1B



