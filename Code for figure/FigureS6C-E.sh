#################################
###########FigureS6C-E###########
#################################
#####################################################
###########1.caOCR number across 0-15 peer###########
#####################################################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATAC_for_review/FigureS6C-E/output")
library(data.table)
library(ggplot2)
statistics=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/all_model_result_statistics/allmodel_type_95sample_caQTL_nominal_200kb_result.statistics")
############
statistics$peer_num=c(0:(nrow(statistics)-1))
p=ggplot(statistics, aes(x = peer_num, y = num_uniquesigpeak)) +   
  geom_point(size = 3, color = "#0073C2B2") +
  labs(  
    x = "PEER factors",  
    y = "OCRs with caQTLs"  
  ) +  
  theme_minimal() +  
  theme(  
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 10),  
    strip.text = element_text(size = 12, face = "bold"),  
    #
    panel.grid = element_blank(),  
    #
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")  
  )

# 
ggsave(p,file="all_model_peer_factor.Scatterplot.pdf", width = 5, height = 5) ##FigureS6C

##########################################################################################################################
###########2.Plot the distribution of distances to caOCR center for the top caSNPs from the 7peer model results###########
##########################################################################################################################
caQTL_0.1=fread("/data1/gy/ATACseq_RWAS/caQTL_STITCH/caQTL_model_result/TMM_exp/15.sex_age_Hp_PC1_5_peer1_7/200kb_nominal/95sample_caQTL_200kb_nominal_result.withFDR.FDR0.1.txt")
caQTL_0.1_top=subset(caQTL_0.1,caQTL_0.1$top_proximal_var==1)
dim(caQTL_0.1_top)
#[1] 13665    15
######
caQTL_0.1_top$dist_from_center_kb=(caQTL_0.1_top$var_start-(caQTL_0.1_top$pheno_start+caQTL_0.1_top$pheno_end)/2)/1000

######within 40kb
nrow(subset(caQTL_0.1_top,abs(dist_from_center_kb)<=40))
#[1] 6691

####distribution histogram
p1 <- ggplot(caQTL_0.1_top, aes(x = dist_from_center_kb)) +  
  geom_histogram(binwidth = 20, boundary = -200, fill = "#0073C2B2", color = "#e9ecef", alpha = 0.7) +
  labs(x = "Distance from OCRs center (kb)",
       y = "Top caQTL number") +  
  theme_minimal() +  
  theme(  
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 10),  
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )+ 
  scale_x_continuous(limits = c(-200, 200), breaks = seq(-200, 200, by = 40))
ggsave(p1,file="caQTL_top_SNP_num.dist_from_OCR_center.Histogram.pdf") ##FigureS6D

####P value scatterplot
p2 <- ggplot(caQTL_0.1_top, aes(x = dist_from_center_kb, y = -log10(p_nominal))) +  
  geom_point(color = "#0073C2B2") +
  labs(x = "Distance from OCRs center (kb)",  
       y = expression(-log[10](P~value))) +  
  theme_minimal() +  
  theme(  
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 10),  
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +  
  scale_x_continuous(limits = c(-200, 200), breaks = seq(-200, 200, by = 40))

ggsave(p2,file="caQTL_top_SNP_Pvalve.dist_from_OCR_center.Scatterplot.pdf") ##FigureS6E
