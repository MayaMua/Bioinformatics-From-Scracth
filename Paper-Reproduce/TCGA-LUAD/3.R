
####ESTIMATE####
# Calculate patient immune score and tumor purity
# Matrix score, Immune score, Tumor score, Tumor purity
# setwd("ESTIMATE")
library(utils)
# rforge <- "http://r-forge.r-project.org"
install.packages("estimate",
                 dependencies=TRUE)
library(estimate)
library(tidyverse)

exp <- read.table("TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#计算免疫评分
filterCommonGenes(input.f = "tpms01A_log2.txt",   #输入文件名
                  output.f = "tpms01A_log2.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("tpms01A_log2.gct",   #刚才的输出文件名
              "tpms01A_log2_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#提取结果并整理
ESTIMATE_result <- read.table("tpms01A_log2_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ESTIMATE_result <- ESTIMATE_result[,-1]   
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]   
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
rownames(ESTIMATE_result) <- colnames(exp)
#保存结果
write.table(ESTIMATE_result, file = "ESTIMATE_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F) 

####生存信息整理####
# Which group has longer survival time?
#xena官网：https://xenabrowser.net/datapages/
#下载生存信息
setwd("TCGA-LUAD")
setwd("Survival_data")
library(tidyverse)

survival <- read.delim("OS.txt", row.names = 1)

survival <- survival[, 2:3]
survival <- survival %>% rownames_to_column('sample')
survival$name <- paste0(survival$sample, 'A')#paste 0: seamless
table(substring(survival$name, 14, 16))
rownames(survival) <- survival$name
survival <- survival[, 2:3]
#合并生存信息与表达谱
tpms01A_log2 <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T) # 肿瘤患者基因表达谱
a <- intersect(colnames(tpms01A_log2),rownames(survival))
table(substr(a,14,16))
exp_01A <- tpms01A_log2[,a]
surv_01A <- survival[a,]
exp_01A <- exp_01A %>% t() %>% as.data.frame()
identical(rownames(exp_01A),rownames(surv_01A))
exp_surv_01A <- cbind(surv_01A,exp_01A)
##保存文件##
write.table(exp_surv_01A,"exp_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

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
           legend.title = "ESTIMATEScore",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()