#####################################################
################Figure2A&FigureS5A###################
#####################################################
/Public/gaoyun/software/R-4.2.0/bin/R
library(relaimpo)
library(dplyr)
library(matrixStats)
library(openxlsx)
library(data.table)
###
system("mkdir -p /data1/gy/ATACseq_RWAS/ATACseq/PVE_analysis/result")
OUTPUTS_dir="/data1/gy/ATACseq_RWAS/ATACseq/PVE_analysis/result"
###
exp=as.data.frame(fread("/data1/gy/ATACseq_RWAS/ATACseq/peak_exp/95sample/TMM/95sample_IterativeOverlapPeakSet.final.TMM_qnorm"))
peakset=exp$Geneid
baseline=read.xlsx("/data1/gy/ATACseq_RWAS/ATACseq/baseline/95sample_baseline.final.xlsx")
rownames(baseline)=baseline$ID
###
na_rows <- which(rowSums(is.na(baseline)) > 0)  
na_rows
#BYES02070

#
na_ids <- baseline$ID[na_rows]  
na_ids
#[1] "BYES02070"

#remove BYES02070
baseline=baseline[-na_rows,]
exp=exp[, !colnames(exp) %in% na_ids]  
##
rownames(exp)=exp$Geneid
exp=as.matrix(exp[,-c(1:6)])

###################################################
####age
baseline <- baseline %>%  
  mutate(age_group = case_when(  
    is.na(age) ~ NA_character_, # 保持 NA 为 NA  
    age < 50 & age >= 40 ~ "0",  
    age < 60 & age >= 50 ~ "1",
    age < 70 & age >= 60 ~ "2",
  ))

####BMI
baseline <- baseline %>%  
  mutate(BMI = case_when(  
    is.na(bmi_calc) ~ "NA",
    bmi_calc <= 23.9 ~ "0",  
    bmi_calc > 23.9 ~ "1"  
  ))

####Gastric_lesion
baseline <- baseline %>%  
  mutate(Gastric_lesion_status = case_when(  
    is.na(Gastric_lesion) ~ NA_character_, 
    Gastric_lesion == 0|Gastric_lesion == 1 ~ "0",  
    Gastric_lesion == 2 ~ "1",
    Gastric_lesion == 3 ~ "2"
  ))
  
####smoke_status
baseline <- baseline %>%  
  mutate(smoke_status = case_when(  
    is.na(smoke) ~ "NA",
    smoke < 2 ~ "0",  
    smoke >= 2 ~ "1"  
  ))
####drink_status
baseline <- baseline %>%  
  mutate(drink_status = case_when(  
    is.na(drink) ~ "NA", 
    drink < 3 ~ "0",  
    drink >= 3 ~ "1"  
  ))
####Tea_consumption
baseline <- baseline %>%  
  mutate(Tea_consumption = case_when(  
    is.na(Tea_consumption) ~ "NA",
    Tea_consumption < 2 ~ "0",  
    Tea_consumption >= 2 ~ "1"  
  ))

####Speed of meal
baseline <- baseline %>%  
  mutate(Speed_of_meal = case_when(  
    is.na(Speed_of_meal) ~ "NA",
    Speed_of_meal == 0 ~ "0",  
    Speed_of_meal == 1 ~ "1",
    Speed_of_meal == 2 ~ "2"
  ))

####Fried foods intake
baseline <- baseline %>%  
  mutate(Fried_foods_intake = case_when(  
    is.na(Fried_foods_intake) ~ "NA",
    Fried_foods_intake < 2 ~ "0",  
    Fried_foods_intake >= 2 ~ "1"
  ))
####Pickled foods intake
baseline <- baseline %>%  
  mutate(Pickled_foods_intake = case_when(  
    is.na(Pickled_foods_intake) ~ "NA",
    Pickled_foods_intake < 2 ~ "0",  
    Pickled_foods_intake >= 2 ~ "1"
  ))
####Consumption of fresh vegetables
baseline <- baseline %>%  
  mutate(Consumption_of_fresh_vegetables = case_when(  
    is.na(Consumption_of_fresh_vegetables) ~ "NA",
    Consumption_of_fresh_vegetables < 4 ~ "0",  
    Consumption_of_fresh_vegetables >= 4 ~ "1"
  ))
####Consumption of fresh vegetables
baseline <- baseline %>%  
  mutate(Consumption_of_fresh_fruits = case_when(  
    is.na(Consumption_of_fresh_fruits) ~ "NA",
    Consumption_of_fresh_fruits < 2 ~ "0",  
    Consumption_of_fresh_fruits >= 2 ~ "1"
  ))
####Spicy foods intake
baseline <- baseline %>%  
  mutate(Spicy_foods_intake = case_when(  
    is.na(Spicy_foods_intake) ~ "NA",
    Spicy_foods_intake < 2 ~ "0",  
    Spicy_foods_intake >= 2 ~ "1"
  ))
####Meat intake
baseline <- baseline %>%
  mutate(Meat_intake = case_when(
    is.na(Meat_intake)                          ~ "NA",
    Meat_intake == 1 | Meat_intake == 2         ~ "1",
    Meat_intake < 1 | Meat_intake > 2           ~ "0"
  ))
####Cereal fiber intake
baseline <- baseline %>%  
  mutate(Cereal_fiber_intake = case_when(  
    is.na(Cereal_fiber_intake) ~ "NA",
    Cereal_fiber_intake < 2 ~ "0",  
    Cereal_fiber_intake >= 2 ~ "1"
  ))

####
baseline=baseline[,c(21,3,22,4,23:25,19,9:18)] 
baseline=baseline[,-12] ##remove a factor with one cluster

#### 
baseline <- as.data.frame(lapply(baseline, function(x) as.numeric(as.character(x))), row.names = row.names(baseline)) 
####
#
baseline_row_names <- row.names(baseline)  
#
exp_column_names <- colnames(exp)  
#
if (all(baseline_row_names == exp_column_names)) {  
    print("TRUE")  
} else {  
    print("FALSE")  
}  
#[1] "TRUE"

####
peak_model = lm(exp[1,] ~ ., data = baseline)
rel_impo = suppressWarnings(calc.relimp(peak_model, type = "lmg"))
##
importances = data.frame(matrix(nrow = nrow(exp), ncol = length(rel_impo$lmg)))
rownames(importances) <- rownames(exp)
colnames(importances) <- names(rel_impo$lmg)

## loop over all peaks
for(i in 1:nrow(exp)){
    #
    if(i%%100 == 0) print(i)
    ## model ~ age_group + sex + smoke + drink + HP_C14 + diagnosis
    peak_model = lm(exp[i,] ~ ., data = baseline)

     out <- tryCatch(
      {
      rel_impo = suppressWarnings(calc.relimp(peak_model, type = "lmg"))
      importances[i,] = rel_impo$lmg
          }, 
      error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      importances[i,] = "NA"
        }
      )  
}
#
write.table(cbind("PeakID"=rownames(importances),importances), paste0(OUTPUTS_dir,"/95sample_peakset.PVE_perpeak.txt"), quote = FALSE,sep="\t",col.names=T,row.names=F)

## calculate col means
importances=read.delim(paste0(OUTPUTS_dir,"/95sample_peakset.PVE_perpeak.txt"),row.names=1)
means <- as.data.frame(colMeans(as.matrix(importances)))
colnames(means)[1] <- "means"
rownames(means) <- colnames(importances)
write.table(cbind("PeakID"=rownames(means),means), paste0(OUTPUTS_dir,"/95sample_peakset.PVE_means.txt"), quote = FALSE,sep="\t",col.names=T,row.names=F)

########plot
library(reshape2) 
library(ggplot2) 
#install.packages("ggridges",repos = c(CRAN="https://mirrors.bfsu.edu.cn/CRAN/"))   
library(ggridges)
setwd("/data1/gy/ATAC_for_review/Figure2A&FigureS5A/output")
######################################
##########1.7 major factors###########Figure2A
######################################
vars=c("age_group","sex","HP_C14","Gastric_lesion_status","smoke_status","drink_status","Tea_consumption")
data=cbind(ID=rownames(importances),importances)
data=data[,c("ID",vars)]
#
long_data <- melt(data, id.vars = "ID", variable.name = "variable", value.name = "value")  
#
custom_colors <- rev(c( "#377EB8","#A65628", "#FF7F00", "#E41A1C", "#984EA3", "#FFFF33", "#4DAF4A"))
#
long_data$variable <- factor(long_data$variable, levels = rev(levels(factor(vars,levels=vars))))
#
mean_data <- long_data %>%  
  group_by(variable) %>%  
  summarise(mean_value = mean(value))  
#
p <- ggplot(long_data, aes(x = value, y = variable, fill = variable)) +  
  geom_density_ridges(alpha = 0.7, scale = 1.8) +  
  labs(x = "percent of variance explained (PVE)", y = "") +  
  scale_fill_manual(values = custom_colors) +  
  guides(fill = "none") +
  scale_x_continuous(limits = c(-0.01, 0.2)) +
  theme(  
    panel.background = element_blank(),  
    panel.grid = element_blank(),  
    panel.border = element_blank(),
    #panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black", size = 12),  
    axis.title = element_text(size = 12),  
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 12)  
  )  +
  #
  geom_segment(data = mean_data, aes(x = mean_value, xend = mean_value,   
               y = variable, yend = as.numeric(variable) + 0.2),   
               color = "black", size = 1) +  
  geom_point(data = mean_data, aes(x = mean_value, y = as.numeric(variable) + 0.2),   
             color = "black", size = 3) + 
  #  
  geom_hline(data = mean_data, aes(yintercept = as.numeric(variable)),   
             color = "black", linetype = "solid", size = 0.5)  
ggsave(p,file="/95sample_peakset.PVE_distribution.7main_vars.ridge_plot.pdf", height = 6, width = 6.5)  ##Figure2A

#######################################
##########2.extra 10 factors###########FigureS5A
#######################################
vars=colnames(importances)[!(colnames(importances) %in% c("age_group","sex","HP_C14","Gastric_lesion_status","smoke_status","drink_status","Tea_consumption"))]
data=cbind(ID=rownames(importances),importances)
data=data[,c("ID",vars)]
# 转换为长格式  
long_data <- melt(data, id.vars = "ID", variable.name = "variable", value.name = "value")  
#颜色设置
custom_colors <- c(
  "#8DD3C7",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#FDB462",
  "#B3DE69",
  "#FCCDE5",
  "#D9D9D9",
  "#BC80BD",
  "#CCEBC5"
)
#
long_data$variable <- factor(long_data$variable, levels = rev(levels(factor(vars,levels=vars))))
#
mean_data <- long_data %>%  
  group_by(variable) %>%  
  summarise(mean_value = mean(value))  
#
p <- ggplot(long_data, aes(x = value, y = variable, fill = variable)) +  
  geom_density_ridges(alpha = 0.7, scale = 1.8) +  
  labs(x = "percent of variance explained (PVE)", y = "") +  
  scale_fill_manual(values = custom_colors) + 
  guides(fill = "none") + 
  scale_x_continuous(limits = c(-0.01, 0.2)) +  
  theme(  
    panel.background = element_blank(),  
    panel.grid = element_blank(),  
    panel.border = element_blank(), 
    #panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black", size = 12),  
    axis.title = element_text(size = 12),  
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 12)  
  )  +
  # 
  geom_segment(data = mean_data, aes(x = mean_value, xend = mean_value,   
               y = variable, yend = as.numeric(variable) + 0.2),   
               color = "black", size = 1) + 
  geom_point(data = mean_data, aes(x = mean_value, y = as.numeric(variable) + 0.2),   
             color = "black", size = 3) + 
  #
  geom_hline(data = mean_data, aes(yintercept = as.numeric(variable)),   
             color = "black", linetype = "solid", size = 0.5)  
ggsave(p,file="/95sample_peakset.PVE_distribution.10extra_vars.ridge_plot.pdf", height = 8.5, width = 7) ##FigureS5A
