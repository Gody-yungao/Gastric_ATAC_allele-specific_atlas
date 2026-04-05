##########################
#########FigureS4C########
##########################
/Public/gaoyun/software/R-4.2.0/bin/R
library(factoextra)
library(ggbiplot)
library(data.table)
library(irlba) 
library(Rtsne)
library(ggplot2)
library(ggrepel) 
#####TMM log2 matrix
TMM_log2_mat=read.delim("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/final_TMM/pantissue.pantissue_OCR.IterativeOverlapPeakSet.final.log2TMM",row.names=1)

###top 250000 variable OCR
library(matrixStats) 
#row variance
variances <- rowVars(as.matrix(TMM_log2_mat))  
#sort
idx_sorted <- order(variances, decreasing = TRUE)  
#top 250k
top_250k_idx <- idx_sorted[1:250000]  
TMM_log2_mat_top250k <- TMM_log2_mat[top_250k_idx, ]  
##transpose
TMM_log2_mat_top250k_transpose=t(TMM_log2_mat_top250k)


###
set.seed(123)  # seed
tsne_out <- Rtsne(  
  TMM_log2_mat_top250k_transpose,
  perplexity = 20,
  max_iter = 10000,
  pca = TRUE,
  theta = 0.5
)  
# tsne_out$Y (x, y)
df_tsne <- data.frame(  
  x = tsne_out$Y[, 1],  
  y = tsne_out$Y[, 2],  
  sample = colnames(TMM_log2_mat_top250k))

###TSNE plot
library(openxlsx)
metadata=read.xlsx("/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/QCsample_metadata/pantissue_allsample_afterQC_metadata.xlsx")
df_tsne=merge(df_tsne,metadata,by.x="sample",by.y="sample_new")
write.xlsx(df_tsne,"/data1/gy/ATACseq_RWAS/ATACseq/pantissue/pantissue_OCR_readcount/final_TMM/pantissue.pantissue_OCR.IterativeOverlapPeakSet.final.log2TMM.tsne.xlsx")

my22colors <- c(
  '#f3715c', '#8e7437', '#74905d', '#494e8f', '#bd6758',
  '#fab27b', '#b2d235', '#DC050C', '#E8601C', '#F1932D',
  '#F6C141', '#F7EE55', '#CAE0AB', '#7BAFDE', '#5289C7',
  '#1965B0', '#D6C1DE', '#B178A6', '#882E72', '#8A2BE2',
  '#FF69B4', '#7CFC00'
)

#
col <- c(
  "adrenal-gland" = '#f3715c',
  "brain" = '#494e8f',
  "colon" = '#fab27b',
  "esophagus" = '#b2d235',
  "fallopian-tube" = '#DC050C',
  "fat" = '#E8601C',
  "heart" = '#F1932D',
  "kidney" = '#F6C141',
  "liver" = '#F7EE55',
  "lung" = '#CAE0AB',
  "muscle" = '#7BAFDE',
  "nerve" = '#5289C7',
  "ovary" = '#1965B0',
  "pancreas" = '#D6C1DE',
  "Retina" = '#B178A6',
  "skin" = '#882E72',
  "spleen" = '#8A2BE2',
  "thyroid-gland" = '#FF69B4',
  "stomach" = '#8e7437'
)

#
library(ggplot2)
library(dplyr) 
setwd("/data1/gy/ATAC_for_review/FigureS4C/output")

# label position
centers <- df_tsne %>%
  group_by(tissue) %>%
  summarize(x = mean(x), y = mean(y))

# t-SNE plot
p=ggplot(df_tsne, aes(x = x, y = y, color = tissue)) +  
  geom_point(size = 2) +  
  labs(  
    title = "t-SNE on Top 250k Variable OCRs",  
    x = "tSNE Dimension 1",  
    y = "tSNE Dimension 2"  
  ) +  
  scale_color_manual(values = col) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) +
  #
  geom_text(data = centers, aes(x = x, y = y, label = tissue), 
            size = 5, 
            vjust = -1, 
            color = "black")
ggsave("pantissue.pantissue_OCR.IterativeOverlapPeakSet.final.log2TMM.tSNEplot.pdf",p,width=7.5,height=6) ##FigureS4C
