#############################
#########FigureS9C-D#########
#############################
################################################
#####1.manhatton plot(Z-score) of bulk TWAS#####
################################################
/Public/gaoyun/software/R-4.2.0/bin/R
suppressMessages(library(data.table))
suppressMessages(library(ggrepel))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
#
d=as.data.frame(fread("/data1/gy/ATACseq_RWAS/TWAS_262/result/TWAS.withFDR.result.txt",header=T))
#
d$P0=d$P1-1000000
#
d=d[order(d$CHR, d$P0),]
d=d[!is.na(d$TWAS.P),]

#
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
Sig_Z_Thresh=3.27
d$Sig_Z_Thresh<-Sig_Z_Thresh
d_sig<-d[which(abs(d$TWAS.Z) > d$Sig_Z_Thresh),]
d_sig<-d_sig[rev(order(abs(d_sig$TWAS.Z))),]
d_sig<-d_sig[!duplicated(d_sig$ID),]
dim(d_sig)
#[1] 45 24

#
mycols=rep(c("gray35","gray72"),60)

##
if(sum(d_sig$TWAS.Z > 0) > 0){
	d_sig_pos<-d_sig[d_sig$TWAS.Z > 0,]
}
	
if(sum(d_sig$TWAS.Z < 0) > 0){
	d_sig_neg<-d_sig[d_sig$TWAS.Z < 0,]
}

##
chr_labs<-as.character(unique(d$CHR))

##
ylimit=max(abs(d$'TWAS.Z'),na.rm=T)+1

##plot
setwd("/data1/gy/ATAC_for_review/FigureS9B-C/output")
p <- ggplot(d, aes(x = pos, y = TWAS.Z, colour = factor(CHR))) +  
  geom_point(size = 0.5) +  
  scale_x_continuous(name = "Chromosome", breaks = ticks, labels = chr_labs) +  
  scale_y_continuous(  
    name = 'Z score',   
    limits = c(-ylimit, ylimit),   
    breaks = seq(0, ylimit, by = 5) %>% c(-rev(seq(0, ylimit, by = 5)), seq(0, ylimit, by = 5))  # 设置Y轴刻度以5为间隔  
  ) +  
  scale_colour_manual(values = mycols, guide = FALSE) +  
  geom_hline(yintercept = 0, colour = "black") +  
  geom_hline(yintercept = Sig_Z_Thresh, colour = "black", linetype = "dashed",linewidth=0.5) +  
  geom_hline(yintercept = -Sig_Z_Thresh, colour = "black", linetype = "dashed",linewidth=0.5) +  
  #
  geom_point(data = d_sig,   
             aes(x = pos, y = TWAS.Z,   
                 fill = "#145390"),   
             colour = "black", size = 1.5, shape = 21) +
  labs(title = "") +
  scale_fill_identity()

##
if(sum(d_sig$TWAS.Z > 0) > 0){
	p<-p+geom_text_repel(data=d_sig_pos, aes(x=pos,y=TWAS.Z, label=gene_name), colour='black', nudge_y=1, size=2.5, force=5, segment.alpha=0.25, ylim=c(Sig_Z_Thresh+0.1,NA),max.overlaps = 50)
}
	
if(sum(d_sig$TWAS.Z < 0) > 0){
	p<-p+geom_text_repel(data=d_sig_neg, aes(x=pos,y=TWAS.Z, label=gene_name), colour='black', nudge_y=-1, size=2.5, force=5, segment.alpha=0.25, ylim=c(NA,-Sig_Z_Thresh-0.1),max.overlaps = 50)
}

##
p <- p +   
  theme_cowplot() +  
  theme(  
    axis.text.x = element_text(size = 7, hjust = 0.5),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
	axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

##
ggsave(p,filename="TWAS_result.zscore.manhattan.pdf",width=7.5,height=4) ##FigureS9C
q()

######################################################
#####2.manhatton plot(Z-score) of epithelium TWAS#####
######################################################
/Public/gaoyun/software/R-4.2.0/bin/R
suppressMessages(library(data.table))
suppressMessages(library(ggrepel))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
#
d=as.data.frame(fread("/data1/gy/ATACseq_RWAS/epiTWAS_203/result/epiTWAS.withFDR.result.txt",header=T))
#
d$P0=d$P1-1000000
#
d=d[order(d$CHR, d$P0),]
d=d[!is.na(d$TWAS.P),]

#
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
Sig_Z_Thresh=3.30
d$Sig_Z_Thresh<-Sig_Z_Thresh
d_sig<-d[which(abs(d$TWAS.Z) > d$Sig_Z_Thresh),]
d_sig<-d_sig[rev(order(abs(d_sig$TWAS.Z))),]
d_sig<-d_sig[!duplicated(d_sig$ID),]
dim(d_sig)
#[1] 10 24

##
mycols=rep(c("gray35","gray72"),60)

##
if(sum(d_sig$TWAS.Z > 0) > 0){
	d_sig_pos<-d_sig[d_sig$TWAS.Z > 0,]
}
	
if(sum(d_sig$TWAS.Z < 0) > 0){
	d_sig_neg<-d_sig[d_sig$TWAS.Z < 0,]
}

##
chr_labs<-as.character(unique(d$CHR))

##
ylimit=max(abs(d$'TWAS.Z'),na.rm=T)+1

##plot
setwd("/data1/gy/ATAC_for_review/FigureS9B-C/output")
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
                 fill = "#145390"),   
             colour = "black", size = 1.5, shape = 21) +
  labs(title = "") +
  scale_fill_identity()

##
if(sum(d_sig$TWAS.Z > 0) > 0){
	p<-p+geom_text_repel(data=d_sig_pos, aes(x=pos,y=TWAS.Z, label=gene_name), colour='black', nudge_y=1, size=2.5, force=5, segment.alpha=0.25, ylim=c(Sig_Z_Thresh+0.1,NA),max.overlaps = 50)
}
	
if(sum(d_sig$TWAS.Z < 0) > 0){
	p<-p+geom_text_repel(data=d_sig_neg, aes(x=pos,y=TWAS.Z, label=gene_name), colour='black', nudge_y=-1, size=2.5, force=5, segment.alpha=0.25, ylim=c(NA,-Sig_Z_Thresh-0.1),max.overlaps = 50)
}

##
p <- p +   
  theme_cowplot() +  
  theme(  
    axis.text.x = element_text(size = 7, hjust = 0.5),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
	axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

##
ggsave(p,filename="epiTWAS_result.zscore.manhattan.pdf",width=7.5,height=4) ##FigureS9D
q()
