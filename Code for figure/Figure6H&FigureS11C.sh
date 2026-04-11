##########################################
###########Figure6H&FigureS11C############
##########################################
##############################
######1.AGS ATAC vs. WGS######
##############################
/Public/gaoyun/software/R-4.2.0/bin/R
#install.packages("ggalluvial")
library(forcats)
library(ggalluvial)
library(dplyr)

# 2x2 matrix
allele_table <- matrix(
  c(13, 15,   # WGS: G, C
    392, 131), # ATAC-seq: G, C (GSE264550: GSM219908)
  nrow = 2,
  byrow = TRUE
)
rownames(allele_table) <- c("WGS", "ATAC-seq")
colnames(allele_table) <- c("G", "C")

##fisher test
fisher_result <- fisher.test(allele_table)
fisher_result$p.value
#[1] 0.001756192

#
sankey_ATAC <- data.frame(
  SeqType = rep(c("WGS", "ATAC-seq"), each = 2),
  Allele = rep(c("G", "C"), 2),
  Count = c(13, 15, 392, 131)
)

#
sankey_ATAC <- sankey_data |>
  dplyr::group_by(SeqType) |>
  dplyr::mutate(Fraction = Count / sum(Count))

#
sankey_ATAC$Allele <- factor(sankey_ATAC$Allele, levels = c("C", "G"))

#
allele_colors <- c("G" = "#C86B85", "C" = "#769FCD")

#########sankeyplot
setwd("/data1/gy/ATAC_for_review/Figure6H&FigureS11C/output")
p_ATAC = ggplot(sankey_ATAC, aes(x = SeqType, y = Fraction, fill = Allele, stratum = Allele, alluvium = Allele)) +
    geom_alluvium(width=0.5,alpha =0.3,knot.pos=0.2,color='white') +
    geom_col(width = 0.5,color='white') +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(limits = c("WGS", "ATAC-seq")) +
    scale_fill_manual(values = allele_colors, breaks = c("G","C")) +
    xlab("rs875179 allele") + ylab("Alelle fraction(%)") + 
    theme_classic() + 
    theme(
      #
      axis.text.x = element_text(  
      size = 10,  
      color = "black",  
      angle = 45,
      hjust = 1
      ),  
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      plot.margin = margin(15, 15, 15, 15),
      panel.grid = element_blank()
)
ggsave(p_ATAC,file="AGS_ATAC.ASCA.sankeyplot.pdf",width=6,height=5) ##FigureS11C
q()

#############################################
######2.AGS GATA4/GATA6 ChIPseq vs. WGS######
#############################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(forcats)
library(ggalluvial)
library(dplyr)
library(ggplot2)

allele_counts <- data.frame(
  SeqType = rep(c("WGS", "GATA4", "GATA6"), each = 2),
  Allele  = rep(c("G", "C"), times = 3),
  Count   = c(13, 15,  ##AGS WGS
               20, 3,  ##AGS GATA4 ChIP-seq (GSE51705: GSM1250896)
               13, 2)  ##AGS GATA6 ChIP-seq (GSE51705: GSM1250897)
)

#
wgs <- allele_counts %>% filter(SeqType == "WGS")

##fisher test
compare_results <- lapply(unique(allele_counts$SeqType)[-1], function(sample) {
  x <- allele_counts %>% filter(SeqType == sample)
  mat <- matrix(
    c(wgs$Count[wgs$Allele=="G"], wgs$Count[wgs$Allele=="C"],
      x$Count[x$Allele=="G"], x$Count[x$Allele=="C"]),
    nrow = 2,
    byrow = TRUE
  )
  rownames(mat) <- c("WGS", sample)
  colnames(mat) <- c("G", "C")
  fisher_p <- fisher.test(mat)$p.value
  data.frame(Comparison = paste("WGS vs", sample), P_value = fisher_p)
})

##
fisher_summary <- do.call(rbind, compare_results)
print(fisher_summary)
#    Comparison     P_value
#1 WGS vs GATA4 0.003331345
#2 WGS vs GATA6 0.020243685

###
sankey_data <- data.frame(
  SeqType = rep(c("WGS", "GATA4", "GATA6"), each = 2),
  Allele  = rep(c("G", "C"), times = 3),
  Count   = c(13, 15,   #AGS WGS
              20, 3,    #AGS GATA4 ChIP-seq (GSE51705: GSM1250896)
              13, 2)    #AGS GATA6 ChIP-seq (GSE51705: GSM1250897)
)

#
sankey_data <- sankey_data %>%
  group_by(SeqType) %>%
  mutate(Fraction = Count / sum(Count))

#
sankey_data$Allele <- factor(sankey_data$Allele, levels = c("C", "G"))

#
allele_colors <- c("G" = "#C86B85", "C" = "#769FCD")

#########sankeyplot
setwd("/data1/gy/ATAC_for_review/Figure6H&FigureS11C/output")
p_sankey <- ggplot(
  sankey_data,
  aes(x = SeqType, y = Fraction, fill = Allele, stratum = Allele, alluvium = Allele)
) +
  geom_alluvium(width = 0.5, alpha = 0.3, knot.pos = 0.2, color = "white") +
  geom_col(width = 0.5, color = "white") +
  scale_x_discrete(limits = c("WGS", "GATA4", "GATA6"), expand = c(0.2, 0.2)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = allele_colors, breaks = c("G", "C")) +
  labs(
    title = "",
    x = "Sample type",
    y = "Allele fraction"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm"),
    plot.margin = margin(15, 15, 15, 15),
    panel.grid = element_blank()
  )
ggsave(p_sankey, file = "AGS_GATA4_GATA6.ASCA.sankeyplot.pdf", width = 5, height = 5) ##Figure6H
q()
