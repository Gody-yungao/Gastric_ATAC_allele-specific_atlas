########################
########FigureS5E#######
########################
#######################################
#############1.DEG analysis############
#######################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(edgeR)

gtf_v19 <- rtracklayer::import('/data1/gy/public/gtf/gencode.v19.annotation.gtf')
head(gtf_v19)
gtf_v19_df <- as.data.frame(gtf_v19)
geneid_df <- dplyr::select(gtf_v19_df,c(gene_name,gene_id))#,gene_biotype
geneid_df<-unique(geneid_df)
geneid_df$new <- geneid_df$gene_id
new_gtf_v19<-geneid_df
new_gtf_v19$gene_id<-sapply(stringr::str_split(new_gtf_v19$gene_id, "\\."), function(v)  return(v[1]))
new_gtf_v19<-unique(new_gtf_v19)
new_gtf_v19 <- new_gtf_v19[,c(1,2)]

count_norm <- as.data.frame(fread("/data1/gy/multistage_RNAseq/xuxianfeng/CombineCounts.FilterLowExpression-MergeMutiSample.TMM.tsv"))
count_norm=count_norm[,c(133,1:132)]
#
colnames(count_norm) <- gsub("_DGC|_IGC", "_GC", colnames(count_norm))

##########################################################################################
class_type <- c( "Normal" , "IM" , "GC")

##########################################################################################
##
condition <- factor(sapply(strsplit(colnames(count_norm)[-1] , "_") , "[" , 2))

##
count_norm_mat <- count_norm[,-1]
rownames(count_norm_mat) <- count_norm$gene_id

##########################################################################################
##
result_diff <- c()

for( k in 1:(length(class_type)-1) ){
  class1 <- class_type[k]
  print(class1)
  
  for( j in (k+1):length(class_type) ){
    class2 <- class_type[j]
    print(class2)
    
    ## Run the Wilcoxon rank-sum test for each gene
    use_sample <- grep( paste0("_",class1 , "|" ,"_", class2) , colnames(count_norm_mat) , value = T)
    count_norm_use <- count_norm_mat[,use_sample]
    conditions_use <- grep( paste0(class1 , "|" , class2) , condition , value = T)
    
    pvalues <- sapply(1:nrow(count_norm_use),function(i){
      data<-cbind.data.frame(gene=as.numeric(t(count_norm_use[i,])),conditions_use)
      p <- wilcox.test(gene~conditions_use, data)$p.value
      return(p)
    })
    fdr=p.adjust(pvalues,method = "BH")
    
    ## Calculate the fold-change for each gene
    conditionsLevel <- c(class1 , class2)
    dataCon1 <- count_norm_use[,c(which(conditions_use==conditionsLevel[1]))]
    dataCon2 <- count_norm_use[,c(which(conditions_use==conditionsLevel[2]))]
    foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
    
    ## Output results based on FDR threshold
    outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
    rownames(outRst)=rownames(count_norm_use)
    outRst=na.omit(outRst)
    
    ##
    outRst$class1 <- class2
    outRst$class2 <- class1
    outRst$meanTMM_class1 <- rowMeans(dataCon2)
    outRst$meanTMM_class2 <- rowMeans(dataCon1)
    outRst$gene_id <- rownames(outRst)
    
    result_diff <- rbind( result_diff , outRst)
  }
}

##########################################################################################
colnames(result_diff) <- c("log2FoldChange" , "pvalue" , "padj" , "class1" , "class2" , "meanTMM_class1" , "meanTMM_class2" , "gene_id")
result_diff=merge(new_gtf_v19,result_diff)

write.table( result_diff , "/data1/gy/multistage_RNAseq/xuxianfeng/diff_gene_result.txt" , row.names = F , sep = "\t" , quote = F )

#############################
#############2.GSEA##########FigureS5E
#############################
conda activate GSVA
R
library(data.table)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
setwd("/data1/gy/ATAC_for_review/FigureS5E/output")
#####################################
#############i.IMvsNormal############
#####################################
##
data <- read.delim("/data1/gy/multistage_RNAseq/xuxianfeng/diff_gene_result.txt",header=T)
data=subset(data,class1=="IM" & class2=="Normal")
head(data)
dim(data)
#[1] 14750     9
data <- na.omit(data)
dim(data)
#[1] 14750     9
##
IDs_ENSEMBL_and_ENTREZID <- bitr(geneID = data$gene_id, 
                                fromType = "ENSEMBL", 
                                toType = "ENTREZID", 
                                OrgDb = "org.Hs.eg.db")
#
IDs_ENSEMBL_and_SYMBOL <- bitr(geneID = data$gene_id, 
                                fromType = "ENSEMBL", 
                                toType = "SYMBOL", 
                                OrgDb = "org.Hs.eg.db")                                

##
data$EZTREZID <- IDs_ENSEMBL_and_ENTREZID[match(data$gene_id,IDs_ENSEMBL_and_ENTREZID$ENSEMBL),2]
data$SYMBOL <- IDs_ENSEMBL_and_SYMBOL[match(data$gene_id,IDs_ENSEMBL_and_SYMBOL$ENSEMBL),2]
dim(data)
#[1] 14750    11
data <- na.omit(data)
dim(data)
#[1] 13786    11
########
#
duplicated_rows <- duplicated(data$EZTREZID) | duplicated(data$SYMBOL)
#
data_unique <- data[!duplicated_rows, ]
dim(data_unique)
#[1] 13773    11

####
data_unique <- data_unique[order(data_unique$log2FoldChange,decreasing = T),]
####
gene_list <- data_unique$log2FoldChange
names(gene_list) <- data_unique$EZTREZID
gene_list <- gene_list[!is.na(names(gene_list))]
length(gene_list)
#[1] 13773

#########################HALLMARK pathway
Hallmark <- read.gmt("/data1/gy/public/GSEA_geneset_new/h.all.v2024.1.Hs.entrez.gmt")
Hallmark_list = split(Hallmark$gene, Hallmark$term)
# 
GSEA_enrichment_hallmark <- GSEA(gene_list, 
                        TERM2GENE = Hallmark,
                        pvalueCutoff = 1,
                        minGSSize = 10, 
                        maxGSSize = 500, 
                        eps = 0, 
                        pAdjustMethod = "BH")
GSEA_enrichment_hallmark_final = setReadable(GSEA_enrichment_hallmark, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")                

##
GSEA_enrichment_hallmark_final_table <- GSEA_enrichment_hallmark_final@result
dim(GSEA_enrichment_hallmark_final_table)
#[1] 50 11
write.csv(GSEA_enrichment_hallmark_final_table,"/data1/gy/multistage_RNAseq/xuxianfeng/GSEA/IMvsNormal.HALLMARK.GSEA_result.csv",row.names=F)

##choose pathway for plot
##################1.HALLMARK_E2F_TARGETS
p1=gseaplot2(GSEA_enrichment_hallmark_final, geneSetID = "HALLMARK_E2F_TARGETS", pvalue_table = T,color = "#66A61E",subplots = 1:2,rel_heights = c(1.5, 0.3))
ggsave(p1,filename="IMvsNormal.HALLMARK_E2F_TARGETS.GSEAplot.pdf",height=2.5,width=4)
##################2.HALLMARK_G2M_CHECKPOINT
p2=gseaplot2(GSEA_enrichment_hallmark_final, geneSetID = "HALLMARK_G2M_CHECKPOINT", pvalue_table = T,color = "#66A61E",subplots = 1:2,rel_heights = c(1.5, 0.3))
ggsave(p2,filename="IMvsNormal.HALLMARK_G2M_CHECKPOINT.GSEAplot.pdf",height=2.5,width=4)


#####################################
#############ii.GCvsNormal############
#####################################
##
data <- read.delim("/data1/gy/multistage_RNAseq/xuxianfeng/diff_gene_result.txt",header=T)
data=subset(data,class1=="GC" & class2=="Normal")
head(data)
dim(data)
#[1] 14750     9
data <- na.omit(data)
dim(data)
#[1] 14750     9
##
IDs_ENSEMBL_and_ENTREZID <- bitr(geneID = data$gene_id, 
                                fromType = "ENSEMBL", 
                                toType = "ENTREZID", 
                                OrgDb = "org.Hs.eg.db")
#
IDs_ENSEMBL_and_SYMBOL <- bitr(geneID = data$gene_id, 
                                fromType = "ENSEMBL", 
                                toType = "SYMBOL", 
                                OrgDb = "org.Hs.eg.db")                                

#
data$EZTREZID <- IDs_ENSEMBL_and_ENTREZID[match(data$gene_id,IDs_ENSEMBL_and_ENTREZID$ENSEMBL),2]
data$SYMBOL <- IDs_ENSEMBL_and_SYMBOL[match(data$gene_id,IDs_ENSEMBL_and_SYMBOL$ENSEMBL),2]
dim(data)
#[1] 14750    11
data <- na.omit(data)
dim(data)
#[1] 13786    11
########
#
duplicated_rows <- duplicated(data$EZTREZID) | duplicated(data$SYMBOL)
#
data_unique <- data[!duplicated_rows, ]
dim(data_unique)
#[1] 13773    11

####
data_unique <- data_unique[order(data_unique$log2FoldChange,decreasing = T),]
####
gene_list <- data_unique$log2FoldChange
names(gene_list) <- data_unique$EZTREZID
gene_list <- gene_list[!is.na(names(gene_list))]
length(gene_list)
#[1] 13773

#########################HALLMARK pathway
Hallmark <- read.gmt("/data1/gy/public/GSEA_geneset_new/h.all.v2024.1.Hs.entrez.gmt")
Hallmark_list = split(Hallmark$gene, Hallmark$term)
#
GSEA_enrichment_hallmark <- GSEA(gene_list,  
                        TERM2GENE = Hallmark,
                        pvalueCutoff = 1, 
                        minGSSize = 10,  
                        maxGSSize = 500, 
                        eps = 0,        
                        pAdjustMethod = "BH")  
GSEA_enrichment_hallmark_final = setReadable(GSEA_enrichment_hallmark, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")                
GSEA_enrichment_hallmark_final_table <- GSEA_enrichment_hallmark_final@result
dim(GSEA_enrichment_hallmark_final_table)
#[1] 50 11
write.csv(GSEA_enrichment_hallmark_final_table,"/data1/gy/multistage_RNAseq/xuxianfeng/GSEA/GCvsNormal.HALLMARK.GSEA_result.csv",row.names=F)

##choose pathway for plot
##################1.HALLMARK_E2F_TARGETS
p1=gseaplot2(GSEA_enrichment_hallmark_final, geneSetID = "HALLMARK_E2F_TARGETS", pvalue_table = T,color = "#66A61E",subplots = 1:2,rel_heights = c(1.5, 0.3))
ggsave(p1,filename="GCvsNormal.HALLMARK_E2F_TARGETS.GSEAplot.pdf",height=2.5,width=4)
##################2.HALLMARK_G2M_CHECKPOINT
p2=gseaplot2(GSEA_enrichment_hallmark_final, geneSetID = "HALLMARK_G2M_CHECKPOINT", pvalue_table = T,color = "#66A61E",subplots = 1:2,rel_heights = c(1.5, 0.3))
ggsave(p2,filename="GCvsNormal.HALLMARK_G2M_CHECKPOINT.GSEAplot.pdf",height=2.5,width=4)

