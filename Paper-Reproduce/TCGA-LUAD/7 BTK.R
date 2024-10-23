
#### Differences in BTK expression in tumor and normal samples ####
#### Bar chart ####
library(tidyverse)
library(survival)
library(survminer)

tpms01A_log2 <- read.table("TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms11A_log2 <- read.table("TCGAdata/tpms11A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

gene <- "BTK"
a <- tpms01A_log2[gene,]
b <- tpms11A_log2[gene,]

# a <- a %>% t() %>% as.data.frame()
# b <- b %>% t() %>% as.data.frame()
# write.csv(a, file = "BTK/BTK_01A.csv")
# write.csv(b, file = "BTK/BTK_11A.csv")
# Combine the data into a single dataframe
BTK_data <- data.frame(expression = c(as.numeric(a), as.numeric(b)),
                       group = factor(c(rep("Normal", length(a)), rep("Tumor", length(b)))))

# Perform a t-test to calculate the p-value between the two groups
p_value <- t.test(as.numeric(a), as.numeric(b))$p.value

# Create the boxplot with jitter points
ggplot(BTK_data, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) + # Create the boxplot
  geom_jitter(shape = 16, position = position_jitter(0.2), color = "black") + # Add jitter points for each sample
  scale_fill_manual(values = c("blue", "red")) + # Set the colors (blue for normal, red for tumor)
  labs(y = "BTK expression", x = "", title = "BTK Expression in Normal vs Tumor") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(BTK_data$expression) + 0.5, label = paste0("p = ", signif(p_value, digits = 3)), size = 6)


#### paired plot ####

tpms01A_log2 <- tpms01A_log2 %>% t() %>% as.data.frame()
tpms11A_log2 <- tpms11A_log2 %>% t() %>% as.data.frame()
rownames(tpms01A_log2) <- substring(rownames(tpms01A_log2), 1, 12)
rownames(tpms11A_log2) <- substring(rownames(tpms11A_log2), 1, 12)
# Intersect the row names to get the common samples
a <- intersect(rownames(tpms01A_log2),rownames(tpms11A_log2))
tpms01A_log2 <- tpms01A_log2[a,]
tpms11A_log2 <- tpms11A_log2[a,]
paired_data <- cbind(tpms11A_log2[,gene], tpms01A_log2[,gene]) 
paired_data <- as.data.frame(paired_data)
write.csv(paired_data, file = "BTK/paired.csv")

# Rename the columns
colnames(paired_data) <- c("Normal", "Tumor")

# Create a long-format dataframe for ggplot
paired_data_long <- data.frame(
  Value = c(paired_data$Normal, paired_data$Tumor),
  Group = rep(c("Normal", "Tumor"), each = nrow(paired_data)),
  Patient = rep(1:nrow(paired_data), times = 2)
)

# Perform paired t-test for statistical significance
p_value <- t.test(paired_data$Normal, paired_data$Tumor, paired = TRUE)$p.value

# Create the paired plot
ggplot(paired_data_long, aes(x = Group, y = Value, group = Patient)) +
  geom_point(aes(color = Group), size = 3) + # Add individual data points
  geom_line(aes(group = Patient), color = "black") + # Connect paired samples
  labs(y = "Expression", x = "", title = "Paired Plot: Normal vs Tumor") +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(paired_data_long$Value) + 0.5, 
           label = paste0("***\nP = ", signif(p_value, 3)), size = 6) +
  scale_color_manual(values = c("blue", "orange"))

dev.off()

#### Survival analysis according to BTK high and low groups ####

surv <- read.table("Survival_data/exp_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv$OS.time <- surv$OS.time/365


# BTK
surv$group <- ifelse(surv$BTK > median(surv$BTK),"High","Low")
class(surv$group)
surv$group <- factor(surv$group, levels = c("Low", "High")) 
class(surv$group)
table(surv$group)

fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)


#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))

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
           legend.title = "BTK",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()



#### 不同分期BTK的表达 ####

clinical.expr01A = read.table("Clinical/clinical.expr01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- "BTK"
clinical_BTK <- cbind(clinical.expr01A[,1:6], clinical.expr01A[,gene])
write.csv(clinical_BTK, file = "BTK/clinical_BTK.csv")

