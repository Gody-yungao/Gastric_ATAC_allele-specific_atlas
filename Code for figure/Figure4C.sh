################################
############Figure4C############
################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
setwd("/data1/gy/ATAC_for_review/Figure4C/output")
#############################################################
#######1.colocalized caOCR-eGene distance distribution#######
#############################################################
####input
caQTL_eQTL_coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv")
#
caQTL_eQTL_coloc <- caQTL_eQTL_coloc %>%  
  mutate(Distance_kb = abs(Distance / 1000)) 
  
#
caQTL_eQTL_coloc <- caQTL_eQTL_coloc %>%  
  mutate(class = ifelse(Distance_kb <= 1, "Proximal", "Distal"))  

#
count_within_1kb <- sum(caQTL_eQTL_coloc$Distance_kb <= 1)  
total_count <- nrow(caQTL_eQTL_coloc)  
proportion_within_1kb <- count_within_1kb / total_count  

#
cat("Distance < 1 kb:", proportion_within_1kb * 100, "%\n")  
#Distance < 1 kb: 14.66147 %

#
caQTL_eQTL_coloc <- caQTL_eQTL_coloc %>%  
  mutate(Distance_kb_offset = ifelse(Distance_kb <= 0, 0.001, Distance_kb))  

#
mycols <- c("Proximal" = "#FB6A4A", "Distal" = "#FE9929")  

#
breaks <- 10^seq(-3, 3, by = 0.2) 

# barplot
p1 <- ggplot(caQTL_eQTL_coloc, aes(x = Distance_kb_offset, fill = class)) +  
  geom_histogram(breaks = breaks, color = "white") + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),  
                labels = c("0", "0.01", "0.1", "1", "10", "100", "1000")) +  
  scale_fill_manual(values = mycols) +  
  labs(x = "Distance to TSS (Kb)", y = "Number of colocalized caOCR-eGene Pairs") +  
  theme_bw() +  
  theme(panel.grid = element_blank()) 
##
ggsave("caOCR-eGene_coloc_PPH4_0.5.distance_distribution_histogram.pdf", plot = p1, width = 6.5, height = 4) ##Figure4C (bottom)

####################################################################################################
#######2.colocalized Proximal(TSS distance<=1kb) & Distal(TSS distance>1kb) pairs' proportion#######
####################################################################################################
##caQTL-eQTL coloc result
coloc=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/coloc_result/caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.chr1_22.csv")
coloc=tidyr::separate(coloc,Peak,into=c("chr","region"),remove=F,sep=":")
coloc=tidyr::separate(coloc,region,into=c("Start","End"),remove=T,sep="-")
dim(coloc)
#[1] 2053   16
##
coloc_promoter=subset(coloc,abs(Distance)<=1000)
dim(coloc_promoter)
#[1] 301  16
coloc_enhancer=subset(coloc,abs(Distance)>1000)
dim(coloc_enhancer)
#[1] 1752   16

##
coloc_promoter$anno="Proximal"
coloc_enhancer$anno="Distal"
coloc_conbined=rbind(coloc_promoter,coloc_enhancer)
#
element_counts=as.data.frame(table(coloc_conbined$anno))
colnames(element_counts) <- c("class", "count") 
#
element_counts$prop <- element_counts$count / sum(element_counts$count)
#
element_counts$lab.ypos <-  cumsum(element_counts$prop) - 0.5*element_counts$prop

element_counts$class <- factor(element_counts$class, levels = c("Proximal","Distal"))  
#
mycols <- c("Proximal" = "#FB6A4A", "Distal" = "#FE9929")  
#pieplot
library(ggplot2)
p2 <- ggplot(element_counts, aes(x = 2, y = prop, fill = class)) +  
  geom_bar(  
    stat = "identity",   
    color = "white", 
    width = 0.8,  
    size = 0.2   
  ) +
  coord_polar(theta = "y", start = 0) + 
  geom_text(aes(y = lab.ypos, label = paste(count, "\n(", scales::percent(prop, accuracy = 0.01),")")), color = "black") +
  scale_fill_manual(values = mycols) + 
  theme_void() +
  xlim(0.5, 2.5) + 
  annotate("text", x = 0.5, y = 0, label = "colocalized\ncaOCR-eGene", size = 6, fontface = "bold", hjust = 0.5, vjust = 0.5) 
ggsave(p2,file="caQTL_eQTL_coloc_200kb_1MB_PPH4_0.5.proximal_distal_anno.pieplot.pdf") ##Figure4C (top)
