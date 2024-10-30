library(tidyverse)


#### TCGA Clinical Data####

load("TCGAdata/luad.gdc_2022.rda")

# Extract clinical data
clinical <- as.data.frame(expquery2@colData) %>%   
  .[!duplicated(.$sample),]
clinical <-clinical[, c("gender","age_at_index","ajcc_pathologic_stage",
                       "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m")]

class(clinical$gender)
class(clinical$age_at_index)
class(clinical$ajcc_pathologic_stage)
class(clinical$ajcc_pathologic_t)
class(clinical$ajcc_pathologic_n)
class(clinical$ajcc_pathologic_m)

table(clinical$gender)
table(clinical$age_at_index)
table(clinical$ajcc_pathologic_stage)
table(clinical$ajcc_pathologic_t)
table(clinical$ajcc_pathologic_n)
table(clinical$ajcc_pathologic_m)


clinical$ajcc_pathologic_stage <- gsub("A","",clinical$ajcc_pathologic_stage) # substitute
clinical$ajcc_pathologic_stage <- gsub("B","",clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_t <- gsub("a","",clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub("b","",clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_m <- gsub("a","",clinical$ajcc_pathologic_m)
clinical$ajcc_pathologic_m <- gsub("b","",clinical$ajcc_pathologic_m)

# Truncate the row name 1-16 index
rownames(clinical) <- substring(rownames(clinical), 1, 16)

# Merge clinical data with gene expression data
exp01A <- read.table("TCGAdata/tpms01A_log2.txt", sep = "\t", row.names = 1,
                     check.names = F, stringsAsFactors = F, header = T)

# Find the common samples
clinical01A <- clinical[colnames(exp01A),]   
exp01A <- exp01A %>% t() %>% as.data.frame()

identical(rownames(clinical01A),rownames(exp01A))

clinical.expr01A <- cbind(clinical01A, exp01A)

write.table(clinical.expr01A,"Clinical/clinical.expr01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# Merge clinical data with ESTIMATE result
ESTIMATE_result <- read.table("ESTIMATE/ESTIMATE_result.txt",sep = "\t",
                              row.names = 1, check.names = F,
                              stringsAsFactors = F, header = T)

identical(rownames(clinical01A),rownames(ESTIMATE_result))

clinical.ESTIMATE_result01A <- cbind(clinical01A, ESTIMATE_result)

write.csv(clinical.ESTIMATE_result01A,file = "Clinical/clinical.ESTIMATE_result01A.csv")

# Visualization tool: https://www.helixlife.cn/class
# or check jupyter notebook under 'Clinical' folder