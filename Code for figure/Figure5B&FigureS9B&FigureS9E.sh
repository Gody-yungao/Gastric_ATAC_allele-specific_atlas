##############################################
#########Figure5B&FigureS9B&FigureS9E#########
##############################################
#####################################
#####1.manhatton plot(coloc PP4)#####FigureS9E
#####################################
/Public/gaoyun/software/R-4.2.0/bin/R
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggrepel)
library(data.table)
library(cowplot)
smy_all = read.csv('/data1/gy/ATACseq_RWAS/RWAS_STITCH/caQTL_GWAS_coloc/coloc.abf_result/RWAS_sig_Peak.caQTL_GWAS_coloc_200kb_withRegion_PPH4.chr1_22.with_rsID.csv')[,-c(2:3)]
Peak = data.frame(fread('/data1/gy/ATACseq_RWAS/caQTL_STITCH/exp/95sample_IterativeOverlapPeakSet.final.TMM_qnorm.bed.gz'))
Peak$X.Chr=substr(Peak$X.Chr,4,nchar(Peak$X.Chr))
#
coloc = smy_all[order(smy_all$PPH4),]
coloc = coloc[which(coloc$HitSNP >= 30),]
coloc$order = 1:nrow(coloc)
#
Peak = Peak[,c('X.Chr','Geneid','Start','End')]
coloc = merge(coloc,Peak,by.x='Peak',by.y='Geneid')
#
d = coloc[,c('X.Chr','Start','End','Peak','PPH4','HitSNP')]
colnames(d)[1] = 'CHR'

#hg19
hg19_size=fread("/data1/gy/public/genome/hg19/hg19.chrom.sizes")
hg19_autosomes <- hg19_size[
  grepl("^chr[0-9]+$", V1)
][
  order(as.numeric(sub("chr", "", V1)))
]
#
hg19_autosomes$CHR <- as.numeric(sub("chr", "", hg19_autosomes$V1))

#
hg19_autosomes$V2 <- as.numeric(hg19_autosomes$V2)
hg19_autosomes$cumlen <- cumsum(hg19_autosomes$V2)
hg19_autosomes$offset <- c(0, head(hg19_autosomes$cumlen, -1))

#
tick_info <- hg19_autosomes %>%
  dplyr::mutate(center = offset + V2 / 2)

ticks <- tick_info$center
tick_labels <- paste0("chr", tick_info$CHR)

ticklim <- c(0, max(hg19_autosomes$cumlen))

#
d <- d %>%
  dplyr::filter(CHR %in% 1:22) %>%
  dplyr::arrange(CHR, Start)

#
d$CHR=as.numeric(d$CHR)
d <- d %>%
  dplyr::left_join(hg19_autosomes[, c("CHR","offset")], by="CHR") %>%
  dplyr::mutate(pos = Start + offset)

##
d_pos<-d[d$PPH4 >= 0.5,]

##
chr_labs<-as.character(unique(d$CHR))

##
ylimit=1

###################################plot
setwd("/data1/gy/ATAC_for_review/Figure5B&FigureS9B&FigureS9E/output")
p <- ggplot(d, aes(x = pos, y = PPH4)) +  
  geom_point(
             aes(fill = "#f08e59"           
                 ,   
             size = HitSNP),
             colour = "black", shape = 21) +
  scale_x_continuous(name = "Chromosome", breaks = ticks, labels = c(1:22), expand = c(0.01, 0.01)) +  
  scale_y_continuous(  
    name = "PP4",
    limits = c(0, ylimit),   
    breaks = seq(0, ylimit, by = 0.25),
    expand = c(0.02, 0.02)
  ) +   
  geom_hline(yintercept = 0.5, colour = "black", linetype = "dashed",linewidth=0.5) +  
  scale_fill_identity()

##add label
p<-p+geom_text_repel(data=d_pos, aes(x=pos,y=PPH4, label=Peak), colour='black', nudge_y=0.03, size=2.5, force=5, segment.alpha=0.25)

##
p <- p +   
  theme_cowplot() +  
  theme(  
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    axis.text.x = element_text(size = 7, hjust = 0.5),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

ggsave(p,filename="RWAS_sig_peaks.coloc_result.manhattan.pdf",width=7.5,height=4) ##FigureS9E
q()

########################################
#####2.manhatton plot(RWAS P-value)#####Figure5B
########################################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATAC_for_review/Figure5B&FigureS9B&FigureS9E/output")
library(data.table)
library(ggrepel)
library(ggplot2))
library(cowplot))
library(dplyr)
#RWAS result
dataframe=as.data.frame(fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter/cv/RWAS.withFDR.result.txt",header=T))

#RWAS with coloc
coloc=read.csv("/data1/gy/ATACseq_RWAS/RWAS_STITCH/caQTL_GWAS_coloc/coloc.abf_result/RWAS_sig_Peak.caQTL_GWAS_coloc_200kb_withRegion_PPH4.csv")
coloc=coloc[,c(1,4)]
d=left_join(dataframe,coloc,by=c("ID"="Peak"))
dim(d)
#[1] 8441   22

#hg19
hg19_size=fread("/data1/gy/public/genome/hg19/hg19.chrom.sizes")
hg19_autosomes <- hg19_size[
  grepl("^chr[0-9]+$", V1)
][
  order(as.numeric(sub("chr", "", V1)))
]
#
hg19_autosomes$CHR <- as.numeric(sub("chr", "", hg19_autosomes$V1))

#
hg19_autosomes$V2 <- as.numeric(hg19_autosomes$V2)
hg19_autosomes$cumlen <- cumsum(hg19_autosomes$V2)
hg19_autosomes$offset <- c(0, head(hg19_autosomes$cumlen, -1))

#
tick_info <- hg19_autosomes %>%
  dplyr::mutate(center = offset + V2 / 2)

ticks <- tick_info$center
tick_labels <- paste0("chr", tick_info$CHR)

ticklim <- c(0, max(hg19_autosomes$cumlen))

# 
d <- d %>%
  dplyr::filter(CHR %in% 1:22) %>%
  dplyr::arrange(CHR, P0)

#
d <- d %>%
  dplyr::left_join(hg19_autosomes[, c("CHR","offset")], by="CHR") %>%
  dplyr::mutate(pos = P0 + offset)

#
d$TWAS.P.log10=-log10(d$TWAS.P)
d_sig<-d[which(d$TWAS.fdr < 0.1),]
d_sig<-d_sig[rev(order(d_sig$TWAS.P.log10)),]
d_sig<-d_sig[!duplicated(d_sig$ID),]
dim(d_sig)
#[1] 58 24
Sig_P_log10_Thresh=-log10(0.0007)
d$Sig_P_log10_Thresh<-Sig_P_log10_Thresh

##
mycols=rep(c("black","gray"),60)

##
d_sig_pos<-d_sig[d_sig$PPH4 >= 0.5,]

#
chr_labs<-as.character(unique(d$CHR))

#
ylimit=max(d$TWAS.P.log10,na.rm=T)+0.5

###################################plot
p <- ggplot(d, aes(x = pos, y = TWAS.P.log10, colour = factor(CHR))) +  
  geom_point(size = 0.5) +  
  scale_x_continuous(name = "Chromosome", breaks = ticks, labels = chr_labs, expand = c(0.01, 0.01)) +  
  scale_y_continuous(  
    name = expression('-log'[10]*'(P)'),
    limits = c(0, ylimit),   
    breaks = seq(0, ylimit, by = 5),
    expand = c(0, 0)
  ) +  
  scale_colour_manual(values = mycols, guide = FALSE) +  
  geom_hline(yintercept = Sig_P_log10_Thresh, colour = "black", linetype = "dashed",linewidth=0.5) +  
  #add coloc PP4 to RWAS significant OCRs
  geom_point(data = d_sig,   
             aes(x = pos, y = TWAS.P.log10,   
                 fill = case_when(  
                   PPH4 >= 0.75 ~ "#CB181D", 
                   PPH4 >= 0.5 ~ "#F16913",
                   PPH4 >= 0.1 ~ "#238443",
                   TRUE ~ "#67A9CF"
                 )),   
             colour = "black", size = 1.5, shape = 21) +
  scale_fill_identity()

##add RWAS significant OCR label
p<-p+geom_text_repel(data=d_sig_pos, aes(x=pos,y=TWAS.P.log10, label=ID), colour='black', nudge_y=1, size=2.5, force=5, segment.alpha=0.25, ylim=c(Sig_P_log10_Thresh+0.1,NA))
	
##
p <- p +   
  theme_cowplot() +  
  theme(  
    axis.text.x = element_text(size = 7, hjust = 0.5),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

#
legend_box = data.frame(  
  x = 0.9,
  y = seq(0.90, 0.81, -0.03),
  color = c("#CB181D", "#F16913", "#238443", "#67A9CF"),
  label = c("0.75-1", "0.5-0.75", "0.1-0.5", "0-0.1")
)  

#
p_withlegend = ggdraw(p) +  
  # 
  geom_point(data = legend_box,   
             aes(x = x, y = y, fill = color),
             size = 1.5,
             shape = 21,
             colour = "black") +
  scale_fill_identity() +
  #
  geom_text(data = legend_box,   
            aes(x = x + 0.01, y = y, label = label),
            hjust = 0,
            size = 2.5) +
  #
  draw_label("Coloc PP4",   
             x = legend_box$x[1],
             y = legend_box$y[1] + 0.03,
             hjust = 0,
             size = 7)

##
ggsave(p_withlegend,filename="RWAS_result.Pvalue.withPPH4.manhattan.pdf",width=8.5,height=4) ##Figure5B
q()

########################################
#####3.manhatton plot(RWAS Z-score)#####FigureS9B
########################################
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATAC_for_review/Figure5B&FigureS9B&FigureS9E/output")
library(data.table)
library(ggrepel)
library(ggplot2))
library(cowplot))
library(dplyr)
#RWAS result
dataframe=as.data.frame(fread("/data1/gy/ATACseq_RWAS/RWAS_STITCH/fusion/output/GC_GWAS/merged_filter/cv/RWAS.withFDR.result.txt",header=T))

#RWAS with coloc
coloc=read.csv("/data1/gy/ATACseq_RWAS/RWAS_STITCH/caQTL_GWAS_coloc/coloc.abf_result/RWAS_sig_Peak.caQTL_GWAS_coloc_200kb_withRegion_PPH4.csv")
coloc=coloc[,c(1,4)]
d=left_join(dataframe,coloc,by=c("ID"="Peak"))
dim(d)
#[1] 8441   22

#hg19
hg19_size=fread("/data1/gy/public/genome/hg19/hg19.chrom.sizes")
hg19_autosomes <- hg19_size[
  grepl("^chr[0-9]+$", V1)
][
  order(as.numeric(sub("chr", "", V1)))
]
#
hg19_autosomes$CHR <- as.numeric(sub("chr", "", hg19_autosomes$V1))

#
hg19_autosomes$V2 <- as.numeric(hg19_autosomes$V2)
hg19_autosomes$cumlen <- cumsum(hg19_autosomes$V2)
hg19_autosomes$offset <- c(0, head(hg19_autosomes$cumlen, -1))

#
tick_info <- hg19_autosomes %>%
  dplyr::mutate(center = offset + V2 / 2)

ticks <- tick_info$center
tick_labels <- paste0("chr", tick_info$CHR)

ticklim <- c(0, max(hg19_autosomes$cumlen))

#
d <- d %>%
  dplyr::filter(CHR %in% 1:22) %>%
  dplyr::arrange(CHR, P0)

#
d <- d %>%
  dplyr::left_join(hg19_autosomes[, c("CHR","offset")], by="CHR") %>%
  dplyr::mutate(pos = P0 + offset)

#
Sig_Z_Thresh=3.402
d$Sig_Z_Thresh<-Sig_Z_Thresh
d_sig<-d[which(abs(d$TWAS.Z) > d$Sig_Z_Thresh),]
d_sig<-d_sig[rev(order(abs(d_sig$TWAS.Z))),]
d_sig<-d_sig[!duplicated(d_sig$ID),]
dim(d_sig)
#[1] 58 25

##
mycols=rep(c("gray35","gray72"),60)

##
if(sum(d_sig$TWAS.Z > 0) > 0){
	d_sig_pos<-d_sig[d_sig$PPH4 >= 0.5 & d_sig$TWAS.Z > 0,]
}
	
if(sum(d_sig$TWAS.Z < 0) > 0){
	d_sig_neg<-d_sig[d_sig$PPH4 >= 0.5 & d_sig$TWAS.Z < 0,]
}

##
chr_labs<-as.character(unique(d$CHR))

##
ylimit=max(abs(dataframe$'TWAS.Z'),na.rm=T)+0.5

############################plot
p <- ggplot(d, aes(x = pos, y = TWAS.Z, colour = factor(CHR))) +  
  geom_point(size = 0.5) +  
  scale_x_continuous(name = "Chromosome", breaks = ticks, labels = chr_labs) +  
  scale_y_continuous(  
    name = 'Z score',   
    limits = c(-ylimit, ylimit),   
    breaks = seq(0, ylimit, by = 5) %>% c(-rev(seq(0, ylimit, by = 5)), seq(0, ylimit, by = 5))
  ) +  
  scale_colour_manual(values = mycols, guide = FALSE) +  
  geom_hline(yintercept = 0, colour = "black") +  
  geom_hline(yintercept = Sig_Z_Thresh, colour = "black", linetype = "dashed",linewidth=0.5) +  
  geom_hline(yintercept = -Sig_Z_Thresh, colour = "black", linetype = "dashed",linewidth=0.5) +  
  #
  geom_point(data = d_sig,   
             aes(x = pos, y = TWAS.Z,   
                 fill = case_when(  
                   PPH4 >= 0.75 ~ "#CB181D",
                   PPH4 >= 0.5 ~ "#F16913",
                   PPH4 >= 0.1 ~ "#238443",
                   TRUE ~ "#67A9CF"
                 )),   
             colour = "black", size = 1.5, shape = 21) +
  scale_fill_identity()

##
if(sum(d_sig$TWAS.Z > 0) > 0){
	p<-p+geom_text_repel(data=d_sig_pos, aes(x=pos,y=TWAS.Z, label=ID), colour='black', nudge_y=1, size=2.5, force=5, segment.alpha=0.25, ylim=c(Sig_Z_Thresh+0.1,NA))
}
	
if(sum(d_sig$TWAS.Z < 0) > 0){
	p<-p+geom_text_repel(data=d_sig_neg, aes(x=pos,y=TWAS.Z, label=ID), colour='black', nudge_y=-1, size=2.5, force=5, segment.alpha=0.25, ylim=c(NA,-Sig_Z_Thresh-0.1))
}

##
p <- p +   
  theme_cowplot() +  
  theme(  
    axis.text.x = element_text(size = 7, hjust = 0.5),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

#
legend_box = data.frame(  
  x = 0.9,
  y = seq(0.90, 0.81, -0.03),
  color = c("#CB181D", "#F16913", "#238443", "#67A9CF"),
  label = c("0.75-1", "0.5-0.75", "0.1-0.5", "0-0.1")
)  

#
p_withlegend = ggdraw(p) +  
  #
  geom_point(data = legend_box,   
             aes(x = x, y = y, fill = color),
             size = 1.5,
             shape = 21,
             colour = "black") +
  scale_fill_identity() +
  #
  geom_text(data = legend_box,   
            aes(x = x + 0.01, y = y, label = label), 
            hjust = 0,
            size = 2.5) +
  #
  draw_label("Coloc PP4",   
             x = legend_box$x[1],
             y = legend_box$y[1] + 0.03,
             hjust = 0,
             size = 7)

##
ggsave(p_withlegend,filename="RWAS_result.zscore.withPPH4.manhattan.pdf",width=7.5,height=4) ##FigureS9B
q()

