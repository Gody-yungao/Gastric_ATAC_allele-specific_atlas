################################
############Figure4A############
################################
############################################
#######1.caSNP&non-caSNP eQTL signal#######
############################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(qqman)
################
caqtl=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/95sample_caQTL_200kb_nominal_result.withFDR.for_caQTL_eQTL_coloc.txt")
eqtl=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/262sample_ciseQTL_1MB_result.cis_qtl_pairs.chr1_22.with_symbol.for_caQTL_eQTL_coloc.txt")
####
filtered_caqtl <- caqtl %>%
  group_by(var_id) %>%    
  filter(p_nominal == min(p_nominal, na.rm = TRUE)) %>% 
  ungroup()
filtered_caqtl=as.data.table(filtered_caqtl)
dim(filtered_caqtl)
#[1] 3768444      16
filtered_caqtl=filtered_caqtl[,c("var_id","FDR")]
caSNP=subset(filtered_caqtl,FDR<0.1)
dim(caSNP)
#[1] 261464      2
noncaSNP=subset(filtered_caqtl,FDR>=0.1)
dim(noncaSNP)
#[1] 3506980       2
#
filtered_eqtl <- eqtl %>%
  group_by(variant_id) %>%      
  filter(pval_nominal == min(pval_nominal, na.rm = TRUE)) %>% 
  ungroup()
filtered_eqtl=as.data.table(filtered_eqtl)
dim(filtered_eqtl)
#[1] 4377101      13
filtered_eqtl=filtered_eqtl[,c("variant_id","pval_nominal")]
####
snplist=fread("/data1/gy/ATACseq_RWAS/colocalization_final_262_STITCH/caQTL_eQTL_coloc/input/caQTL_eQTL_intersected_SNPlist.txt")
snplist=snplist[,c("caqtlID","eqtlID")]

#############################caSNP observed eQTL P-value
caSNP=inner_join(snplist,caSNP,by=c("caqtlID"="var_id"))
caSNP=inner_join(caSNP,filtered_eqtl,by=c("eqtlID"="variant_id"))
dim(caSNP)
#[1] 261464      4
noncaSNP=inner_join(snplist,noncaSNP,by=c("caqtlID"="var_id"))
noncaSNP=inner_join(noncaSNP,filtered_eqtl,by=c("eqtlID"="variant_id"))
dim(noncaSNP)
#[1] 3506980       4

############################expected P-value
###caSNP
caSNP = caSNP[order(caSNP$pval_nominal,decreasing = F),]
caSNP$obs=-log10(caSNP$pval_nominal)
caSNP$exp=-log10(ppoints(length(caSNP$pval_nominal)))
caSNP=caSNP[,c('obs','exp')]
caSNP$group="caSNP"
###noncaSNP
#
set.seed(123)
#
sampled_noncaSNP <- noncaSNP %>%
  sample_n(size = nrow(caSNP))
#
sampled_noncaSNP = sampled_noncaSNP[order(sampled_noncaSNP$pval_nominal,decreasing = F),]
sampled_noncaSNP$obs=-log10(sampled_noncaSNP$pval_nominal)
sampled_noncaSNP$exp=-log10(ppoints(length(sampled_noncaSNP$pval_nominal)))
sampled_noncaSNP=sampled_noncaSNP[,c('obs','exp')]
sampled_noncaSNP$group="sampled_noncaSNP"

######
data=rbind(caSNP,sampled_noncaSNP)
data$group=factor(data$group,level=c("caSNP","sampled_noncaSNP"))

######################
#######2.QQplot#######
######################
##
setwd("/data1/gy/ATAC_for_review/Figure4A/output")
library(ggplot2)
# 
p <- ggplot(data, aes(x = exp, y = obs)) +  
  geom_point(aes(color = group), alpha = 0.8, size = 3, shape = 16) + 
  geom_abline() +  
  scale_x_continuous(expand = c(0.01, 0.01)) +  
  #scale_y_continuous(breaks = seq(0, 11, by = 2.5)) +  
  labs(x = expression(paste("Expected -log"[10], "(", italic('P'), ")")),  
       y = expression(paste("Observed -log"[10], "(", italic('P'), ")"))) +  
  theme_bw() +  
  theme(panel.grid = element_blank(),  
        panel.border = element_blank(), 
        axis.line = element_line(color = "black"),  
        axis.text = element_text(size = 12, color = 'black', face = 'plain'),  
        axis.title = element_text(size = 12, color = 'black', face = 'plain'),  
        legend.text = element_text(size = 12), 
        legend.title = element_blank(), 
        # legend.position = c(0.2, 0.8)
        legend.position = "none") + 
        scale_color_manual(values = c("caSNP" = "#EF3B2C", "sampled_noncaSNP" = "black"))  
ggsave(p,file="caQTL.eQTL_signal.QQplot.tiff", dpi = 600)  ##Figure4A

##KS test
ks.test(caSNP$obs, sampled_noncaSNP$obs)
#p-value < 2.2e-16
