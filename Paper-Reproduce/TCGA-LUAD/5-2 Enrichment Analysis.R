#### Enrichment analysis with intersecting differential genes ####

library(tidyverse)
library(DESeq2)
library("BiocManager")

# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# "org.Hs.eg.db" mainly annotated human genes: 
# for translation between different database IDs
library(org.Hs.eg.db)
library(clusterProfiler)

# install.packages("ggnewscale")
library(ggnewscale)

# Load the differential expression results

DEG_final <- read_csv("Immune_Stromal_DEG/DEG_final.csv")
load("Immune_Stromal_DEG/DEG_ImmuneScore.rda")
DEG <- as.data.frame(res)
DEG <- DEG[DEG_final$SYMBOL,]
DEG <- rownames_to_column(DEG, "SYMBOL")
# Convert Gene ID. 
genelist <- bitr(DEG$SYMBOL, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG, genelist, by="SYMBOL")

#### GO #### biological process
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff = 0.05, 
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego, ego_res, file = "Immune_Stromal_DEG/GO_DEG_final.Rda")

#### KEGG #### pathway
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff = 0.1)
kk_res <- kk@result
save(kk,kk_res,file = "Immune_Stromal_DEG/KEGG_DEG_final.Rda")

# Convert Entrez IDs back to gene symbols for KEGG results
kk_res$geneSymbol <- mapIds(org.Hs.eg.db, kk_res$geneID, 
                            keytype = "ENTREZID", column = "SYMBOL")

# Update the KEGG result object to use symbols
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Network Graph

List = DEG$log2FoldChange
names(List) = DEG$ENTREZID
head(List)
# Gene ID: 1770      974      920      695      933     4057 
#         1.485408 1.498562 1.274762 1.416403 1.477717 1.086343 
List = sort(List,decreasing = T)

# GO
cnetplot(ego, foldChange = List, circular = TRUE, colorEdge = TRUE)
# KEGG
cnetplot(kk, foldChange = List, circular = TRUE, colorEdge = TRUE)
