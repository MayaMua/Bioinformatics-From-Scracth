# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos=rforge,
#                  dependencies=TRUE)
# install.packages("survminer")
# install.packages("survival")

library(utils)
library(estimate)
library(tidyverse)
library(survival)
library(survminer)

#### ESTIMATE####
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

# Combine survival information with ESTIMATE
ESTIMATE_result <- read.table("ESTIMATE/ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
identical(rownames(ESTIMATE_result),rownames(surv_01A))
ESTIMATE_result_surv_01A <- cbind(surv_01A, ESTIMATE_result)

##保存文件##
write.table(ESTIMATE_result_surv_01A, "Survival_data/ESTIMATE_result_surv_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#### Survival Analysis ####

surv <- read.table("Survival_data/ESTIMATE_result_surv_01A.txt",
                   sep = "\t", row.names = 1,
                   check.names = F,
                   stringsAsFactors = F,header = T)
# Covert time to year
surv$OS.time <- surv$OS.time/365
time_col = "OS.time"
status_col = "OS"
# Define the function to plot survival curves for different scores
plot_survival_curve <- function(surv, score_column, score_name, time_col, status_col) {
  
  # Step 1: Split the data into High and Low groups based on the median of the specified score
  med <- median(surv[[score_column]], na.rm = TRUE)
  surv$group <- ifelse(surv[[score_column]] > med, "High", "Low")
  surv$group <- factor(surv$group, levels = c("Low", "High"))
  
  # Step 2: Perform log-rank test (survdiff) to compare survival between the groups
  fitd <- survdiff(Surv(as.numeric(surv[[time_col]]), as.numeric(surv[[status_col]])) ~ group, 
                   data = surv, na.action = na.exclude)
  pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ", round(pValue, 3))))
  
  # Step 3: Fit the survival curve
  fit <- survfit(Surv(as.numeric(surv[[time_col]]), as.numeric(surv[[status_col]])) ~ group, data = surv)
  
  # Step 4: Plot the survival curve using ggsurvplot
  surv_plot <- ggsurvplot(fit,                      # The fitted survival curve object (output of survfit())
                          data = surv,               # The dataset containing the survival data (same as used to fit the model)
                          pval = p.lab,              # Displays the p-value on the plot (typically for log-rank test)
                          conf.int = TRUE,           # Adds confidence interval for the survival curves
                          risk.table = TRUE,         # Shows a risk table below the plot, displaying the number of subjects at risk at different time points
                          risk.table.col = "strata", # Colors the risk table according to different strata (groups) in the data
                          palette = "jco",           # Specifies the color palette (in this case, "jco" is a predefined palette from the Journal of Clinical Oncology)
                          legend.labs = c("Low", "High"), # Custom labels for the legend, indicating different groups (Low and High, for example, based on immune score levels)
                          size = 1,                  # Sets the line size (thickness) for the survival curves
                          xlim = c(0, 20),           # Sets the limits of the x-axis to range from 0 to 20 (typically time in years or months)
                          break.time.by = 5,         # Sets the intervals for the x-axis ticks, breaking time every 5 units
                          legend.title = "ImmuneScore", # Sets the title of the legend, describing what the legend labels represent
                          surv.median.line = "hv",   # Adds a horizontal and vertical line at the median survival time (h = horizontal, v = vertical)
                          ylab = "Survival probability (%)", # Label for the y-axis, showing the probability of survival in percentages
                          xlab = "Time (Years)",     # Label for the x-axis, showing time in years
                          ncensor.plot = TRUE,       # Adds a plot showing the number of censoring events over time
                          ncensor.plot.height = 0.25, # Sets the height of the censor plot (as a fraction of the main plot)
                          risk.table.y.text = FALSE  # Hides the y-axis labels of the risk table for a cleaner look
  )
  
  # Return the plot
  return(surv_plot)
}

# Example usage for ImmuneScore
plot_survival_curve(surv = surv, score_column = "ImmuneScore", 
                    score_name = "ImmuneScore", 
                    time_col = "OS.time", 
                    status_col = "OS")

# Example usage for StromalScore
plot_survival_curve(surv = surv, score_column = "StromalScore", 
                    score_name = "StromalScore", 
                    time_col = "OS.time", 
                    status_col = "OS")

# Example usage for ESTIMATEScore
plot_survival_curve(surv = surv, score_column = "ESTIMATEScore", 
                    score_name = "ESTIMATEScore", 
                    time_col = "OS.time", 
                    status_col = "OS")

