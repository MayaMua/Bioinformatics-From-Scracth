library(DESeq2)
library(tidyverse)
library(pheatmap)

# Function to perform differential analysis and prepare data for heatmap plotting
run_differential_expression <- function(counts_file, estimate_file, score_type, logFC_cutoff = 1, save_file) {
  
  # Step 1: Load the counts data (01A) and ESTIMATE result
  counts_01A <- read.table(counts_file, sep = "\t", row.names = 1, 
                           check.names = F, stringsAsFactors = F, header = T)
  
  estimate <- read.table(estimate_file, sep = "\t", row.names = 1, 
                         check.names = F, header = T)
  
  # Step 2: Process data based on the provided score (ImmuneScore or StromalScore)
  med <- as.numeric(median(estimate[, score_type], na.rm = TRUE))
  estimate <- as.data.frame(t(estimate))
  
  # Ensure that the sample names match between counts and estimate data
  if (!identical(colnames(counts_01A), colnames(estimate))) {
    stop("Column names in counts data and estimate data do not match!")
  }
  
  # Step 3: Define conditions for the differential analysis
  conditions <- data.frame(sample = colnames(counts_01A),
                           group = factor(ifelse(estimate[score_type, ] > med, "high", "low"),
                                          levels = c("low", "high"))) %>%
    column_to_rownames("sample")
  
  # Step 4: Check if the differential expression results already exist
  if (file.exists(save_file)) {
    message("Loading existing differential expression results from file...")
    load(save_file)
  } else {
    message("Running DESeq differential analysis...")
    
    # Step 5: Perform differential expression analysis
    dds <- DESeqDataSetFromMatrix(countData = counts_01A, colData = conditions, design = ~ group)
    dds <- DESeq(dds)
    res <- results(dds)
    
    # Save the differential expression results
    save(res, file = save_file)
    print(paste("Saved differential expression results to", res_file))
    
  }
  
  # Step 6: Add differential expression status (UP, DOWN, NOT)
  DEG <- as.data.frame(res)
  type1 <- (DEG$padj < 0.05) & (DEG$log2FoldChange < -logFC_cutoff)
  type2 <- (DEG$padj < 0.05) & (DEG$log2FoldChange > logFC_cutoff)
  print("Differential expression status added to DEG data frame.")
  
  DEG$change <- ifelse(type1, "DOWN", ifelse(type2, "UP", "NOT"))
  
  # Step 7: Extract differentially expressed genes
  gene_up <- filter(DEG, change == "UP")
  gene_down <- filter(DEG, change == "DOWN")
  gene_up_down <- rbind(gene_up, gene_down)
  gene_up_down.rownames <- rownames(gene_up_down)
  print("=============================================")
  print(paste("Number of upregulated genes:", nrow(gene_up)))
  print(paste("Number of downregulated genes:", nrow(gene_down)))
  print("=============================================")
  
  # Step 8: Load TPM expression data
  exp <- read.table("TCGAdata/tpms01A_log2.txt", sep = "\t", row.names = 1, 
                    check.names = F, stringsAsFactors = F, header = T)
  
  exp_diff <- exp[gene_up_down.rownames, ]
  
  # Step 9: Set up the grouping information and adjust column order
  annotation_col <- conditions
  level_high <- filter(annotation_col, group == 'high')
  level_low <- filter(annotation_col, group == 'low')
  exp_diff_high <- exp_diff[, rownames(level_high)]
  exp_diff_low <- exp_diff[, rownames(level_low)]
  exp_diff <- cbind(exp_diff_high, exp_diff_low)
  
  print(paste("Finished diff express function for", score_type))
  # Return the results for further use in plotting
  return(list(exp_diff = exp_diff, annotation_col = annotation_col))
}

# Example usage for ImmuneScore
immune_results <- run_differential_expression(
  counts_file = "TCGAdata/counts01A.txt", 
  estimate_file = "ESTIMATE/ESTIMATE_result.txt", 
  score_type = "ImmuneScore", 
  save_file = "DEG Plots/DEG_ImmuneScore.rda"
)

# Plot the heatmap for ImmuneScore
pheatmap(immune_results$exp_diff,
         annotation_col = immune_results$annotation_col,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize = 10,
         fontsize_row = 3,
         fontsize_col = 3,
         main = "Differential Expression Heatmap: ImmuneScore" 
)


# Example usage for StromalScore
stromal_results <- run_differential_expression(
  counts_file = "TCGAdata/counts01A.txt", 
  estimate_file = "ESTIMATE/ESTIMATE_result.txt", 
  score_type = "StromalScore", 
  save_file = "DEG Plots/DEG_StromalScore.rda"
)

# Plot the heatmap for StromalScore
pheatmap(stromal_results$exp_diff,
         annotation_col = stromal_results$annotation_col,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize = 10,
         fontsize_row = 3,
         fontsize_col = 3,
         main = "Differential Expression Heatmap: StromalScore"
         )



# library(BiocManager)
# # BiocManager::install('DESeq2')
# library(DESeq2)
# library(tidyverse)
# # install.packages("pheatmap")
# library(pheatmap)
# 
# #### Differential Expression ####
# #### Get Immune Score ####
# 
# # TCGA differential analysis done with counts Read 01A 
# # since it's grouping patient 01A for differential analysis
# counts_01A <- read.table("TCGAdata/counts01A.txt",sep = "\t",
#                          row.names = 1, check.names = F,
#                          stringsAsFactors = F, header = T)
# 
# # Immunization score grouping
# estimate <- read.table("ESTIMATE/ESTIMATE_result.txt", sep = "\t",
#                        row.names = 1, check.names = F, header = T)
# # Process Data
# # Median 
# x <- "ImmuneScore"
# med <- as.numeric(median(estimate[, x])) 
# estimate <- as.data.frame(t(estimate)) # Transpose
# 
# identical(colnames(counts_01A), colnames(estimate))
# 
# # information of separate groups
# conditions <- data.frame(sample = colnames(counts_01A),
#                       group = factor(ifelse(
#                       estimate[x,] > med, "high", "low"),
#                       levels = c("low","high"))) %>% 
#                       column_to_rownames("sample")
# 
# #conditions <- data.frame(sample=colnames(counts_01A),
# #                      group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high")))
# #conditions <- column_to_rownames(conditions,"sample")
# 
# #### Differential Analysis ####
# # Two files are needed to do the variance analysis, 
# # one is the counts expression spectrum file 
# # and the other is the classification file.
# # Step 1: Create a DESeqDataSet object
# dds <- DESeqDataSetFromMatrix(
#   # Only need counts data and conditions
#   countData = counts_01A,
#   colData = conditions,
#   design = ~ group)
# 
# # Step 2: Perform differential analysis
# dds <- DESeq(dds)
# 
# # Step 3: Extract results
# resultsNames(dds)
# res <- results(dds)
# save(res, file = "DEG Plots/DEG_ImmuneScore.rda")
# 
# #### Plot Heatmap ####
# # log2FoldChange: The log2 fold change of the gene expression between the two groups.
# # Positive values indicate upregulation, and negative values indicate downregulation. 
# # High absolute values indicate high fold changes.
# # Low absolute values indicate low fold changes, 
# # but high upregulation in the low group.
# # padj: The adjusted p-value of the gene expression.
# 
# DEG <- as.data.frame(res)
# 
# # TPMS
# exp <- read.table("TCGAdata/tpms01A_log2.txt",
#                   sep = "\t", row.names = 1,
#                   check.names = F, stringsAsFactors = F, header = T)
# # Set the up and down regulation threshold
# logFC_cutoff <- 1
# type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
# type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
# # If the gene is down-regulated (change < -1), it is marked as DOWN.
# # Else if the gene is up-regulated (change > 1), it is marked as UP.
# # If the gene is not differentially expressed, it is marked as NOT.
# DEG$change = ifelse(type1,"DOWN", ifelse(type2, "UP", "NOT"))
# table(DEG$change)
# # DOWN   NOT    UP 
# # 556 18199   733 
# 
# # Extract differentially expressed gene expression spectrum
# gene_up <- filter(DEG, change == 'UP')
# gene_down <- filter(DEG, change == 'DOWN')
# gene_up_down <- rbind(gene_up, gene_down)
# gene_up_down.rownames <- rownames(gene_up_down)
# exp_diff <- exp[gene_up_down.rownames, ]
# 
# # Set the grouping information
# annotation_col <- conditions
# # Process the order of the columns in exp_diff
# level_high <- filter(annotation_col, group == 'high') 
# level_low <- filter(annotation_col, group == 'low')
# exp_diff_high <- exp_diff[, rownames(level_high)]
# exp_diff_low <- exp_diff[, rownames(level_low)]
# exp_diff <- cbind(exp_diff_high, exp_diff_low)
# 
# # Plotting
# pheatmap(exp_diff,
#          annotation_col = annotation_col,
#          scale = "row",
#          show_rownames = F,
#          show_colnames = F,
#          color = colorRampPalette(c("navy", "white", "red"))(50),
#          cluster_cols = F,
#          cluster_rows = T,
#          fontsize = 10,
#          fontsize_row = 3,
#          fontsize_col = 3
#          )
# 
# dev.off()
# 
# #### Stromal Score ####
# 
# x <- "StromalScore"
# med <- as.numeric(median(estimate[, x]))
# estimate <- as.data.frame(t(estimate))
# identical(colnames(counts_01A),colnames(estimate))
# 
# conditions=data.frame(sample=colnames(counts_01A),
#                       group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high"))) %>% 
#   column_to_rownames("sample")
# #差异分析准备工作
# dds <- DESeqDataSetFromMatrix(
#   countData = counts_01A,
#   colData = conditions,
#   design = ~ group)
# 
# #开始差异分析
# dds <- DESeq(dds)
# #这句很重要
# resultsNames(dds)
# #提取结果
# res <- results(dds)
# save(res,file="DEG_StromalScore.Rda")
# 
# ####热图绘制####
# DEG <- as.data.frame(res)
# #读取表达谱
# exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# #添加上下调信息
# logFC_cutoff <- 1
# type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
# type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
# DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
# table(DEG$change)
# 
# library(pheatmap)
# #提取差异基因表达谱
# a <- filter(DEG, change == 'UP')
# b <- filter(DEG, change == 'DOWN')
# c <- rbind(a, b)
# d <- rownames(c)
# exp_diff <- exp[d, ]
# #设置分组信息
# annotation_col <- conditions
# #对exp_diff 列的顺序进行处理
# a <- filter(annotation_col, group == 'high')
# b <- filter(annotation_col, group == 'low')
# exp_diff_high <- exp_diff[, rownames(a)]
# exp_diff_low <- exp_diff[, rownames(b)]
# exp_diff <- cbind(exp_diff_high, exp_diff_low)
# #开始画图
# pheatmap(exp_diff,
#          annotation_col=annotation_col,
#          scale = "row",
#          show_rownames = F,
#          show_colnames =F,
#          color = colorRampPalette(c("navy", "white", "red"))(50),
#          cluster_cols =F,
#          cluster_rows = T,
#          fontsize = 10,
#          fontsize_row=3,
#          fontsize_col=3)
# #保存图片 调整大小
# dev.off()#关闭画板