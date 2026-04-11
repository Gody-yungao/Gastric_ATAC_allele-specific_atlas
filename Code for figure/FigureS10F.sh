#################################
###########FigureS10F############
#################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(data.table)
library(dplyr)
library(forestploter)
library(meta)
##
chr=4
bp=48073300
variant_pattern <- paste0("^", chr, ":", bp, "(:|$)")  

#
data_dir <- "/data1/gy/EAS_GWAS_meta/per_GWAS/"  
pre_files <- list.files(  
  path = data_dir,   
  pattern = "_inINDEL_preMETA$",
  full.names = TRUE  
)  

#
subset_list <- lapply(pre_files, function(f) {  
  tryCatch({  
    #
    gwas_id <- sub("_inINDEL_preMETA$", "", basename(f))  
    
    #
    dt <- fread(f, showProgress = FALSE)[  
      grepl(variant_pattern, MarkerName)  
    ]  
    
    #
    if (nrow(dt) > 0) {  
      dt[, GWAS := gwas_id]  
    } else {  
      warning("No variants found in: ", f)  
      dt <- data.table(GWAS = character())  
    }  
    
    return(dt)  
  }, error = function(e) {  
    message("Error processing ", f, ": ", e$message)  
    return(data.table())
  })  
})  

#
merged_data <- rbindlist(subset_list, use.names = TRUE, fill = TRUE) 

#
gwas_order <- c("BJ", "NJ", "NCI", "ONCO", "BF", "NF", "BBJ")  

#
merged_data[, GWAS := factor(GWAS, levels = gwas_order)]  
setorder(merged_data, GWAS) 

#
merged_data[, `:=`(  
  OR = exp(BETA),  
  CI_lower = exp(BETA - 1.96 * SE),  
  CI_upper = exp(BETA + 1.96 * SE)  
)]  

#
merged_data[, c("OR", "CI_lower", "CI_upper") :=   
  lapply(.SD, round, 2),   
.SDcols = c("OR", "CI_lower", "CI_upper")]  

##
res <- metagen(
       data = merged_data,
       TE = log(merged_data$OR),
       lower = log(merged_data$CI_lower),
       upper = log(merged_data$CI_upper),
       sm = "OR",
       studlab = paste0(merged_data$GWAS," GWAS"),
       title = "rs875179",
       random = F) 

## forestplot
setwd("/data1/gy/ATAC_for_review/FigureS10F/output")
pdf('rs875179_7GWAS_meta.forestplot.pdf',width = 10,height = 7) ##FigureS10F
forest(res, xlim=c(0.75,1.5),col.square = "black", col.diamond = "grey", col.diamond.lines = "black", hetstat = T, spacing = 1.5)
dev.off()
q()
