################################
############FigureS2A###########
################################
#############################################
###1.ENCODE stomach DNase/ATAC preparation###
#############################################
mkdir -p /data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_raw
cd /data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_raw
######DNase
##53Y ENCSR006IMH
wget https://www.encodeproject.org/files/ENCFF731YGM/@@download/ENCFF731YGM.bed.gz --no-check-certificate
gunzip ENCFF731YGM.bed.gz
mv ENCFF731YGM.bed ENCBS515JRU_DNase_hg38.bed
##37Y ENCSR177NIJ
wget https://www.encodeproject.org/files/ENCFF221TMD/@@download/ENCFF221TMD.bed.gz --no-check-certificate
gunzip ENCFF221TMD.bed.gz
mv ENCFF221TMD.bed ENCBS384ZDL_DNase_hg38.bed
##51Y ENCSR641ZPF
wget https://www.encodeproject.org/files/ENCFF151DGI/@@download/ENCFF151DGI.bed.gz --no-check-certificate
gunzip ENCFF151DGI.bed.gz
mv ENCFF151DGI.bed ENCBS329RLP_DNase_hg19.bed
##54Y ENCSR163PKT
wget https://www.encodeproject.org/files/ENCFF229RBJ/@@download/ENCFF229RBJ.bed.gz --no-check-certificate
gunzip ENCFF229RBJ.bed.gz
mv ENCFF229RBJ.bed ENCBS294WHX_DNase_hg19.bed
##3Y ENCSR246PXX
wget https://www.encodeproject.org/files/ENCFF701IYS/@@download/ENCFF701IYS.bed.gz --no-check-certificate
gunzip ENCFF701IYS.bed.gz
mv ENCFF701IYS.bed ENCBS578HBL_DNase_hg19.bed
##34Y ENCSR782SSS
wget https://www.encodeproject.org/files/ENCFF941HUI/@@download/ENCFF941HUI.bed.gz --no-check-certificate
gunzip ENCFF941HUI.bed.gz
mv ENCFF941HUI.bed ENCBS441WEO_DNase_hg19.bed
######ATACseq
##53Y ENCBS296KUN
wget https://www.encodeproject.org/files/ENCFF880EHI/@@download/ENCFF880EHI.bed.gz --no-check-certificate
gunzip ENCFF880EHI.bed.gz
mv ENCFF880EHI.bed ENCBS296KUN_ATAC_hg38.bed
##54Y ENCBS424ANS
wget https://www.encodeproject.org/files/ENCFF028CAI/@@download/ENCFF028CAI.bed.gz --no-check-certificate
gunzip ENCFF028CAI.bed.gz
mv ENCFF028CAI.bed ENCBS424ANS_ATAC_hg38.bed
######rm bed.gz
rm -r *.bed.gz

##simple bed
mkdir ../input_new
for file in /data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_raw/*.bed; do  
  filename=$(basename "$file")  
  awk -F'\t' '{print $1 "\t" $2 "\t" $3}' "$file" > "../input_new/${filename}"  
done  

##hg38 liftover to hg19
cd /data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_new
######DNase
##53Y 
/data1/gy/software/LiftOver/liftOver \
  ENCBS515JRU_DNase_hg38.bed \
  /data1/gy/software/LiftOver/hg38ToHg19.over.chain \
  ENCBS515JRU_DNase_hg19.bed \
  ENCBS515JRU_DNase_hg19.unMapped.bed
##37Y
/data1/gy/software/LiftOver/liftOver \
  ENCBS384ZDL_DNase_hg38.bed \
  /data1/gy/software/LiftOver/hg38ToHg19.over.chain \
  ENCBS384ZDL_DNase_hg19.bed \
  ENCBS384ZDL_DNase_hg19.unMapped.bed
######ATACseq
##53Y
/data1/gy/software/LiftOver/liftOver \
  ENCBS296KUN_ATAC_hg38.bed \
  /data1/gy/software/LiftOver/hg38ToHg19.over.chain \
  ENCBS296KUN_ATAC_hg19.bed \
  ENCBS296KUN_ATAC_hg19.unMapped.bed
##54Y
/data1/gy/software/LiftOver/liftOver \
  ENCBS424ANS_ATAC_hg38.bed \
  /data1/gy/software/LiftOver/hg38ToHg19.over.chain \
  ENCBS424ANS_ATAC_hg19.bed \
  ENCBS424ANS_ATAC_hg19.unMapped.bed
#
rm -r *_hg38.bed
rm -r *_hg19.unMapped.bed

######keep chr1-22
mkdir ../input_chr1_22_rmdup
#
for bed_file in /data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_new/*.bed; do  
    #
    base_name=$(basename "$bed_file" .bed)  
    
    #
    output_file="/data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_chr1_22_rmdup/${base_name}_chr1_22.rmdup.bed"  

     #
    awk -F'\t' '($1 ~ /^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$/) {  
        #
        print $1 "\t" $2  "\t" $3 "\t" $1 ":" ($2 + 1) "-" $3  
    }' "$bed_file" | sort -u -k1,1 -k2n | awk -F'\t' '{ print $1 "\t" $2 "\t" $3 "\t" $4 }' > "$output_file" 

######rm blacklist
conda activate /Public/gaoyun/miniconda3/envs/ATACseq_ASE
cd /data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_chr1_22_rmdup
mkdir ../input_final
#
for bed_file in /data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_chr1_22_rmdup/*_chr1_22.rmdup.bed; do 
    #
    sample_name=$(basename "$bed_file" .bed|cut -d "_" -f 1-2)
    bash /data1/gy/code/Bash/iterative_peak_filtering/remove_peaks_overlapping_blacklisted_regions.sh \
    $bed_file \
    ../input_final/${sample_name}.bed \
    /data1/gy/public/genome/hg19/hg19.chrom.sizes \
    /data1/gy/public/blacklist/consensusBlacklist.bed
done

######################################
###2.overlap with ENCODE sample OCR###
######################################
######
conda activate /Public/gaoyun/miniconda3/envs/ATACseq
mkdir ../95sample_IterativeOverlapPeakSet_overlap
cd /data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/input_final
bedtools intersect \
  -a /data1/gy/ATACseq_RWAS/ATACseq/macs2/iterative_peak_filtering/7.combined_extend_peak_501bp_rmBlacklist_Iteratively_filter_Nomarlized_Iteratively_filter_2rep_SPM4/95sample/95sample_IterativeOverlapPeakSet.sort.bed \
  -b *.bed \
  -C -filenames \
> ../95sample_IterativeOverlapPeakSet_overlap/95sample_IterativeOverlapPeakSet.intersected_with_ENCODE_stomach_DNase_ATAC
#
awk 'BEGIN{FS=OFS="\t"} \
{split($6, a, "."); $6=a[1]; print}' \
../95sample_IterativeOverlapPeakSet_overlap/95sample_IterativeOverlapPeakSet.intersected_with_ENCODE_stomach_DNase_ATAC \
> ../95sample_IterativeOverlapPeakSet_overlap/95sample_IterativeOverlapPeakSet.intersected_with_ENCODE_stomach_DNase_ATAC.final  
##
rm -r ../95sample_IterativeOverlapPeakSet_overlap/95sample_IterativeOverlapPeakSet.intersected_with_ENCODE_stomach_DNase_ATAC


###################################################################
###3.overlap with overall and single-sample ENCODE OCR & barplot###
###################################################################
######
/Public/gaoyun/software/R-4.2.0/bin/R
setwd("/data1/gy/ATAC_for_review/FigureS2A/output")
library(data.table)  
library(dplyr)  
library(tidyr)  
OCR <- fread("/data1/gy/ATACseq_RWAS/ATACseq/public_stomach_OCR/95sample_IterativeOverlapPeakSet_overlap/95sample_IterativeOverlapPeakSet.intersected_with_ENCODE_stomach_DNase_ATAC.final", header = F)  
colnames(OCR) <- c("CHROM", "START", "END", "ID", "SCORE", "FEATURE", "VALUE")  
# 
OCR_wide <- OCR %>%  
  pivot_wider(names_from = FEATURE, values_from = VALUE, values_fill = list(VALUE = 0))  

#####################
######i.overall######
#####################
#
sample_columns <- c("ENCBS294WHX_DNase", "ENCBS296KUN_ATAC",   
                    "ENCBS329RLP_DNase", "ENCBS384ZDL_DNase",   
                    "ENCBS424ANS_ATAC", "ENCBS441WEO_DNase",   
                    "ENCBS515JRU_DNase", "ENCBS578HBL_DNase")  

#
# rowSums(...) > 0 
overlap_count <- sum(rowSums(OCR_wide[, sample_columns] > 0) > 0)  
overlap_count
#[1] 82946
novel_count <- nrow(OCR_wide) - overlap_count  
novel_count
#[1] 26241
total_count <- nrow(OCR_wide)  

#
overlap_pct <- overlap_count / total_count * 100 
overlap_pct
#[1] 75.96692 
novel_pct   <- novel_count / total_count * 100  
novel_pct
#[1] 24.03308

#
df_summary <- data.frame(  
  Category   = c("Overlap", "Novel"),  
  Count      = c(overlap_count, novel_count),  
  Percentage = c(overlap_pct, novel_pct)  
)  

##plot
library(ggplot2)
p1=ggplot(df_summary, aes(x = factor(1), y = Percentage, fill = Category)) +  
  geom_bar(stat = "identity", width = 0.5) +  
  #
  scale_fill_manual(values = c("Overlap" = "#878787", "Novel" = "#FDAE61")) +  
  #
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  
  labs(  
    x = NULL,  
    y = "Percent Overlap\n stomach ATAC vs ENCODE stomach DNase/ATAC",  
    fill = NULL  
  ) +  
  #
  theme_minimal(base_size = 14) +  
  theme(  
    panel.grid = element_blank(),      
    panel.background = element_blank(), 
    plot.background = element_blank(),
    panel.border = element_blank(),   
    axis.line = element_line(size = 0.5),
    axis.ticks.y = element_line(size = 0.5),   
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()  
  )  
ggsave(p1,filename="95sample_IterativeOverlapPeakSet.intersected_with_ENCODE_stomach_DNase_ATAC.total.barplot.pdf",width=4,height=6) ##FigureS2A (left)


#########################
######ii.per sample######
#########################
#  
dnase_cols <- c("ENCBS294WHX_DNase", "ENCBS329RLP_DNase",   
                "ENCBS384ZDL_DNase", "ENCBS441WEO_DNase",   
                "ENCBS515JRU_DNase", "ENCBS578HBL_DNase") 

atac_cols <- c("ENCBS296KUN_ATAC", "ENCBS424ANS_ATAC")  

#
calc_sample_overlap <- function(df, sample_name) {  
  overlap_count <- sum(df[[sample_name]] > 0)  
  total_count   <- nrow(df)  
  return(overlap_count / total_count * 100)  
}  

#
df_proportions <- data.frame(  
  sample      = c(dnase_cols, atac_cols),  
  proportion  = c(  
    sapply(dnase_cols, calc_sample_overlap, df = OCR_wide),  
    sapply(atac_cols, calc_sample_overlap, df = OCR_wide)  
  )  
) %>%  
  mutate(  
    type = ifelse(sample %in% dnase_cols, "DNase", "ATAC")  
  )  
#
df_proportions <- df_proportions %>%  
  rowwise() %>%  
  mutate(count = sum(OCR_wide[[sample]] > 0)) %>%  
  ungroup()  
  

#
df_proportions$sample <- factor(df_proportions$sample, levels = c(dnase_cols, atac_cols))  

#
p2 <- ggplot(df_proportions, aes(x = sample, y = proportion, fill = type)) +  
  geom_col(width = 0.7) +  
  geom_text(  
    aes(label = paste0(count, "\n", round(proportion, 2), "%")),  
    vjust = -0.3,  
    size = 3  
  ) +  
  scale_fill_manual(values = c("DNase" = "#67A9CF", "ATAC" = "#78C679")) +  
  scale_y_continuous(  
    limits = c(0, 100),  
    breaks = seq(0, 100, 20),  
    labels = function(x) paste0(x, "%"),  
    expand = c(0, 0)  
  ) +  
  labs(  
    x  = NULL,  
    y  = "Percent Overlap\n stomach ATAC vs ENCODE stomach DNase/ATAC",  
    fill = NULL  
  ) +  
  theme_minimal(base_size = 14) +  
  theme(  
    panel.grid      = element_blank(),  
    panel.background= element_blank(),  
    plot.background = element_blank(),  
    panel.border    = element_blank(),  
    axis.line       = element_line(size = 0.5),  
    axis.ticks.y    = element_line(size = 0.5),  
    #
    axis.ticks.x    = element_line(size = 0.5),  
    #
    axis.text.x     = element_text(angle = 30, hjust = 1)  
  )  
ggsave(p2,filename="95sample_IterativeOverlapPeakSet.intersected_with_ENCODE_stomach_DNase_ATAC.persample.barplot.pdf",width=7,height=6) ##FigureS2A (right)

###############biosample id converted to ENCODE assay accession with Adobe Illustrator
