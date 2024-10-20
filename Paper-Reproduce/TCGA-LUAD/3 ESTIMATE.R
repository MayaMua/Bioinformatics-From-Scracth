# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge,
#                  dependencies=TRUE)
# install.packages("survminer")
# install.packages("survival")

library(utils)
library(estimate)
library(tidyverse)

####ESTIMATE####
# Calculate patient immune score and tumor purity
# Stromal score, Immune score, Tumor score, Tumor purity

exp <- read.table("TCGAdata/tpms01A_log2.txt",
                  sep = "\t",
                  row.names = 1,
                  check.names = F,
                  stringsAsFactors = F,
                  header = T)

# Get the common genes
# Common genes are filtered to ensure that the analysis 
# only includes genes shared across the dataset, 
# ensuring accuracy and comparability in the results.

filterCommonGenes(input.f = "TCGAdata/tpms01A_log2.txt",   
                  output.f = "ESTIMATE/tpms01A_log2.gct", 
                  id = "GeneSymbol")  

estimateScore("ESTIMATE/tpms01A_log2.gct",   # input
              "ESTIMATE/tpms01A_log2_estimate_score.txt",  # output
              platform="affymetrix")

# The ESTIMATE score is used to quantify the stromal and immune content in tumor samples and to infer tumor purity.
# ImmuneScore indicates immune cell infiltration.
# StromalScore indicates the presence of stromal cells.
# ESTIMATEScore provides an overall estimate of tumor purity.
# Tumor purity and the tumor microenvironment (immune and stromal content) a
# re essential for understanding cancer behavior, 
# predicting prognosis, and tailoring treatments such as immunotherapy.

ESTIMATE_result <- read.table("ESTIMATE/tpms01A_log2_estimate_score.txt", sep = "\t", 
                              row.names = 1, check.names = F,
                              stringsAsFactors = F, header = T)

# Remove the first column
ESTIMATE_result <- ESTIMATE_result[,-1]  
# Set the column names
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]   
# Transpose the matrix
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
# Set the row names (. -> -)
rownames(ESTIMATE_result) <- colnames(exp)
# Save the result
write.table(ESTIMATE_result, file = "ESTIMATE/ESTIMATE_result.txt",
            sep = "\t", row.names = T,
            col.names = NA, quote = F) 


#### Survival data process ####
# Which group has longer survival time?
# https://xenabrowser.net/datapages/
# Download survival data

# In terminal:
# mkdir Survival_data && cd Survival_data
# wget -c https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FLUAD_survival.txt
# mv *survival*.txt LUAD_survival.txt
survival <- read.delim("Survival_data/LUAD_survival.txt", row.names = 1)
# Only need OS and OS.time
# 0: dead 1: alive
survival <- survival[, 2:3]
# Add a new column 'sample'
survival <- survival %>% rownames_to_column('sample')
# Add 'a' to the end of the sample name and save the result in a new column 'name'
survival$name <- paste0(survival$sample, 'A') # paste0: seamless
table(substring(survival$name, 14, 16))
rownames(survival) <- survival$name
survival <- survival[, 2:3]


# Combine survival information with gene expression profile
tpms01A_log2 <- read.table("TCGAdata/tpms01A_log2.txt", sep = "\t",
                           row.names = 1, check.names = F, header = T) # Tumor patient gene expression profile
# Get the common samples
intersec_samples <- intersect(colnames(tpms01A_log2), rownames(survival))
table(substr(intersec_samples, 14, 16))
exp_01A <- tpms01A_log2[, intersec_samples]
surv_01A <- survival[intersec_samples,]
exp_01A <- exp_01A %>% t() %>% as.data.frame()
identical(rownames(exp_01A),rownames(surv_01A))
exp_surv_01A <- cbind(surv_01A,exp_01A)

write.table(exp_surv_01A,"Survival_data/exp_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#合并生存信息与ESTIMATE
ESTIMATE_result <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
identical(rownames(ESTIMATE_result),rownames(surv_01A))
ESTIMATE_result_surv_01A <- cbind(surv_01A,ESTIMATE_result)

##保存文件##
write.table(ESTIMATE_result_surv_01A,"ESTIMATE_result_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####根据ESTIMATE_result高低组做生存分析####
setwd("TCGA-LUAD")
setwd("survival")
surv <- read.table("ESTIMATE_result_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv$OS.time <- surv$OS.time/365

#median 
#ImmuneScore
surv$group <- ifelse(surv$ImmuneScore > median(surv$ImmuneScore), "High", "Low")
class(surv$group)
surv$group <- factor(surv$group, levels = c("Low", "High")) 
class(surv$group)
table(surv$group)
# install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ", round(pValue, 3))))
# install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           palette = "jco", # color palette jco
           legend.labs = c("Low", "High"), 
           size = 1,
           xlim = c(0,20), # x axis length
           break.time.by = 5, # x step length is 5
           legend.title = "ImmuneScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

#StromalScore
surv$group <- ifelse(surv$StromalScore > median(surv$StromalScore),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,20), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "StromalScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

#ESTIMATEScore
surv$group <- ifelse(surv$ESTIMATEScore > median(surv$ESTIMATEScore),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,20), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "ESTIMATEScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()