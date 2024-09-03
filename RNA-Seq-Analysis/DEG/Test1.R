
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("GenomicFeatures")
#                      #Installation of DESeq2                     
                     # BiocManager::install("DESeq2")
                     

cnt <- read.csv("counts.csv")
str(cnt)
met <- read.csv("metadata2.csv")
str(met)
all(colnames(cnt) %in% rownames(met))

# checking order of row names and column names
all(colnames(cnt) == rownames(met))

library(DESeq2)
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

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

res0.05.df <- as.data.frame(res0.05)
str(res0.05.df)

res0.05.df$Symbol <- mapIds(org.Hs.eg.db, rownames(res0.05.df), 
                            keytype = "ENSEMBL", column = "SYMBOL")
res0.05.df
str(res0.05.df)
write.csv(res0.05.df, "final_test_udemy.csv")

# PCA plot
vsd <- vst(deg, blind = FALSE)
plotPCA(vsd, intgroup = "dexamethasone")

# size factor estimation
sizeFactors(deg)

# disperssion
plotDispEsts(deg)
