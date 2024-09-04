
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("GenomicFeatures")
#                      #Installation of DESeq2                     
                     # BiocManager::install("DESeq2")
                     

cnt <- read.csv("Data/counts.csv")
str(cnt)
met <- read.csv("Data/metadata2.csv")
str(met)

# Make sure that the column names of the count matrix are the same as the row names of the metadata
all(colnames(cnt) %in% rownames(met))
# Check order of row names and column names
all(colnames(cnt) == rownames(met))

library(DESeq2)

###
# Two Steps:
# 1. Create a DESeqDataSet object
# 2. Perform DEG analysis
###

dds <- DESeqDataSetFromMatrix(countData = cnt, 
                              colData = met,
                              design = ~dexamethasone)
dds

# Removal of low count reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds

# Setting reference for DEG Analysis
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

deg <- DESeq(dds)
res <- results(deg)
# write.csv(res, "test_udemy.csv")

# Summary statistics
summary(res)

res0.05 <- results(deg, alpha = 0.05)
summary(res0.05)

# BiocManager::install("org.Hs.eg.db") 
library(org.Hs.eg.db)

res0.05.df <- as.data.frame(res0.05)
str(res0.05.df)

res0.05.df$Symbol <- mapIds(org.Hs.eg.db, rownames(res0.05.df), 
                            keytype = "ENSEMBL", column = "SYMBOL")
res0.05.df
str(res0.05.df)
write.csv(res0.05.df, "Data/final_test_udemy.csv")

# PCA plot
vsd <- vst(deg, blind = FALSE)
plotPCA(vsd, intgroup = "dexamethasone")

# Size factor estimation
sizeFactors(deg)

# Disperssion
plotDispEsts(deg)

# MA plot
plotMA(res0.05)

# Best Genes
library(dplyr)
best_genes <- res0.05.df %>% 
  arrange(padj) %>% 
  head(10)

best_genes
write.csv(best_genes, "Outputs/best_genes.csv")

# Volcano plot
library(ggplot2)
# ggplot(res0.05.df, aes(x = log2FoldChange, y = -log10(pvalue))) +
#   geom_point(aes(color = padj < 0.05)) +
#   theme_minimal() +
#   geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#   labs(x = "Log2 Fold Change", y = "-Log10 p-value") +
#   theme(legend.position = "none")

vol <- res0.05.df %>% 
  filter(!is.na(padj))
ggplot(vol, aes(x = log2FoldChange, y = -log10(padj), 
                color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  geom_point() + 
  geom_text(data = best_genes, aes(label = Symbol), hjust = -0.2, vjust = 0.5)



# Heatmap
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
# top_genes <- rownames(best_genes)
# top_genes
# mat <- counts(dds)[top_genes, ]
# mat <- mat / rowSums(mat) * 1e6
# pheatmap(mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)

top_genes <- res0.05.df %>% 
  arrange(padj) %>% 
  head(30)

mat <- counts(deg, normalized = TRUE)[rownames(top_genes), ]
head(mat, 5)
mat.z <- t(apply(mat, 1, scale))
head(mat.z, 5)
met
colnames(mat.z) <- rownames(met)
Heatmap(mat.z, cluster_rows = TRUE, cluster_columns = TRUE, 
        column_labels = colnames(mat.z), row_labels = top_genes$Symbol)

# Pathway Analysis
