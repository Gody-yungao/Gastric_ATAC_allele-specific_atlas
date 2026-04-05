##############################
#########FigureS6A-B##########
##############################
#####################################
###1.20sample dosage cor calculate###
#####################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
#
indir="/data1/gy/ATACseq_RWAS/genotype/ASAonchip_vs_bam_checkgt/checkgt"
files <- list.files(indir, pattern = "\\.compare\\.snps\\.txt$", full.names = TRUE)
bn <- basename(files) 
sampleIDs <- sub("\\..*$", "", bn) 

# GT to DS func
library(dplyr)
library(stringr)
dosage_from_gt <- function(x) {
  x <- as.character(x)
  out <- rep(NA_real_, length(x))
  bad <- is.na(x) | x %in% c(".", "./.", ".|.")
  ok  <- !bad
  alleles <- strsplit(x[ok], "[/|]")
  out[ok] <- vapply(alleles, function(a) {
    if (length(a) != 2) return(NA_real_)
    if (!all(a %in% c("0","1"))) return(NA_real_)
    sum(as.integer(a))
  }, numeric(1))
  out
}

##
res <- data.frame(
  sample = character(0),
  cor = numeric(0),
  stringsAsFactors = FALSE
)
for (sampleID in sampleIDs){
  t = as.data.frame(fread(paste0(indir,"/",sampleID,".compare.snps.txt"),header=FALSE))
  # require that a genotype call was made in the gold-standard vcf
  t = t[t[,5]!="./.",]
  
  #ASA GT to DS
  t <- t %>%
  mutate(
    V5_dosage = dosage_from_gt(V5) 
  )
  
  # cor analysis
  r <- cor(t$V4, t$V5_dosage, method = "pearson", use = "complete.obs")
  res <- rbind(res, data.frame(sample = sampleID, cor = r))
}
write.table(res,
            file = "/data1/gy/ATACseq_RWAS/genotype/ASAonchip_vs_bam_checkgt/plot/20sample.ATACbam.vs.ASA.dosage_cor.res", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################
###2.20sample dosage cor distribution barplot###
################################################
library(ggplot2)
library(dplyr)
library(scales)
library(grid)
setwd("/data1/gy/ATAC_for_review/FigureS6A-B/output")
p <- ggplot(res, aes(x = cor)) +
  geom_histogram(
    aes(fill = ..x..),
    bins = 15, color = "white", linewidth = 0.4
  ) +
  scale_fill_gradientn(
    colours = c("#6BAED6", "#3182BD", "#6A51A3"),
    values  = rescale(c(min(res$cor, na.rm = TRUE), mu, max(res$cor, na.rm = TRUE))),
    guide   = "none"
  ) +
  geom_vline(xintercept = mu, linetype = "dashed", color = "#2C3E50", linewidth = 0.8) +
  annotate(
    "label", x = mu, y = Inf, label = sprintf("Mean = %.3f", mu),
    vjust = 1.3, hjust = -0.05, size = 3.5, label.size = 0,
    fill = "white", color = "#2C3E50"
  ) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = seq(x_min, x_max, by = 0.02),
    minor_breaks = seq(x_min, x_max, by = 0.01),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(6),
    minor_breaks = waiver(),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = "Genotype correlation coefficient\n(imputation vs. array)",
    y = "Number of individuals",    
    title = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey25", linewidth = 0.5),
    axis.ticks = element_line(color = "grey25", linewidth = 0.4),
    axis.ticks.length = unit(3, "pt"),
    axis.text = element_text(color = "black")
  )

ggsave("20sample.ATACbam.vs.ASA.dosage_cor.plot.pdf", plot = p, width = 5, height = 6) ##FigureS6A


####################
###3.20sample ROC###
####################
# summarize accuracy of imputation for het snps, with array-based genotypes as ground truth
library(data.table)
library(dplyr)
library(stringr)
#
indir="/data1/gy/ATACseq_RWAS/genotype/ASAonchip_vs_bam_checkgt/checkgt"
files <- list.files(indir, pattern = "\\.compare\\.snps\\.txt$", full.names = TRUE)
bn <- basename(files)
sampleIDs <- sub("\\..*$", "", bn)
#
for (sampleID in sampleIDs){
  t = as.data.frame(fread(paste0(indir,"/",sampleID,".compare.snps.txt"),header=FALSE))
  # require that a genotype call was made in the gold-standard vcf
  t = t[t[,5]!="./.",]

  x = t[,5] == "0/1" | t[,5] == "1/0" # find true hets 
  y = t[,4] # gene dosage from sequencing data

  res <- lapply(seq(0, 1, 0.01), function(s) {
    keep <- abs(y - 1) < s
    sens <- sum(x[keep]) / sum(x)
    spec <- mean(x[keep])
    data.frame(s = s, sensitivity = sens, precision = spec)
  })
  res <- do.call(rbind, res)
  fwrite(res,paste0(indir,"/",sampleID,".compare.ASAhetsnps.ROC.res"),sep="\t",quote=F,row.names=F,col.names=T)
}

#########################
###4.20sample ROC plot###
#########################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

setwd("/data1/gy/ATAC_for_review/FigureS6A-B/output")
ROCfiles <- list.files(indir, pattern = "\\.compare\\.ASAhetsnps\\.ROC\\.res$", full.names = TRUE)

#
plotfile <- "20sample.ASAhetSNP.ROCplot.pdf"
statfile <- str_replace(plotfile, "plot\\.pdf$", ".stats")

#
read_one <- function(f) {
  df <- as.data.frame(fread(f, header = TRUE))
  #
  if (!all(c("sensitivity", "specificity") %in% names(df))) {
    if (ncol(df) < 3) stop("at least need 3 cols：threshold, sensitivity, specificity")
    colnames(df)[2] <- "sensitivity"
    colnames(df)[3] <- "specificity"
  }
  df$sample <- str_replace(basename(f), "\\.compare\\.ASAhetsnps\\.ROC\\.res$", "")
  df
}

roc_df <- rbindlist(lapply(ROCfiles, read_one), use.names = TRUE, fill = TRUE) %>%
  mutate(FPR = 1 - specificity,
         TPR = sensitivity) 

#
get_sens_at_spec <- function(spec, sens, target = 0.9, method = c("nearest", "which_max")) {
  method <- match.arg(method)
  if (method == "nearest") {
    idx <- which.min(abs(spec - target))
    return(sens[idx])
  } else {
    idxs <- which(spec > target)
    if (length(idxs) == 0) return(NA_real_)
    return(sens[max(idxs)])
  }
}

stats <- roc_df %>%
  group_by(sample) %>%
  summarize(
    sens80 = get_sens_at_spec(specificity, sensitivity, target = 0.8, method = "which_max"),
    sens90 = get_sens_at_spec(specificity, sensitivity, target = 0.9, method = "which_max"),
    .groups = "drop"
  )

#
stats <- stats %>%
  mutate(
    sens80 = ifelse(is.na(sens80), 0, sens80),
    sens90 = ifelse(is.na(sens90), 0, sens90)
  )

m90 <- mean(stats$sens90)

#####plot
p <- ggplot(roc_df, aes(x = FPR, y = sensitivity, color = sample)) +
  geom_line(alpha = 0.85, linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = m90, color = "red3", linetype = "dashed") +
  annotate("text", x = 0.98, y = m90, label = sprintf("mean sens@spec=0.9 = %.3f", m90),
           hjust = 1, vjust = -0.2, size = 3.2, color = "red3") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_color_discrete(guide = guide_legend(ncol = 1, title = "Sample")) +
  labs(
    title = "ROC for het SNP detection",
    x = "1-Specificity",
    y = "Sensitivity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    #
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #
    plot.title = element_text(face = "bold", hjust = 0.5),
    #
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(3, "pt"),
    #
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.6),
    legend.position = "right",
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.4, "cm")
  )

#
ggsave(plotfile, plot = p, width = 7, height = 6) ##FigureS6B
write.table(stats %>% rename(sens_at_spec_0.8 = sens80, sens_at_spec_0.9 = sens90),
            file = statfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)












