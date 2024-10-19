setwd("TCGAdata")
library(tidyverse)
load("luad.gdc_2022.rda")
load("gene_annotation_2022.rda")

# table: Creates a frequency table of the different gene types.
table(gene_annotation_2022$type)

# Gene name: symbol ENSEMBL
# Differential Expression only need counts.
# tpms = Transcripts Per Million
# unstranded: counts
# tpm_unstranded = tpms

#### counts ####
# Extracting and Transforming Counts Data
# counts: Extracts the counts data from the SummarizedExperiment object.
# expquery2@assays@data@listData[["unstranded"]]: Accesses the unstranded counts data.
# colnames(counts): Sets the column names of the counts matrix.
# rownames(counts): Sets the row names of the counts matrix.
counts <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES

# Converts the counts matrix to a data frame.
# Adds the row names as a new column named "ENSEMBL".
# Joins the counts data with gene annotation data on the "ENSEMBL" column.
# Removes duplicated rows based on the "symbol" column.
# %>% (The Pipe Operator)
counts_annot <- counts %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]

# make the 'symbol' column as row name
# rownames(counts_annot) <- NULL # Delete row names
# counts_annot <- column_to_rownames(counts_annot, var = 'symbol')

rownames(counts_annot) <- NULL
counts_annot <- counts_annot %>% column_to_rownames("symbol") 

# Extract mRNA 
#table
table(counts_annot$type)
counts_annot_protein_coding <- counts_annot[counts_annot$type == "protein_coding",]
# Remove the first and last columns
counts_annot_protein_coding <- counts_annot_protein_coding[,-c(1, 
                                            ncol(counts_annot_protein_coding))]

ncol(counts_annot_protein_coding) # # of columns (602)
nrow(counts_annot_protein_coding) # # of rows (19934)

# Only keep 16 characters of TCGA barcode and remove duplicate obs.
# e.g. TCGA-73-4658-01A-01R-1755-07 -> TCGA-73-4658-01A
# Overwrite the column names
colnames(counts_annot_protein_coding) <- substring(
                colnames(counts_annot_protein_coding), 1, 16)
# Remove duplicated columns
counts_annot_protein_coding <- counts_annot_protein_coding[
          ,!duplicated(colnames(counts_annot_protein_coding))]

table(substring(colnames(counts_annot_protein_coding), 14, 16))

# Extract 01A (cancer)
counts01A <- counts_annot_protein_coding[,
             substring(colnames(counts_annot_protein_coding),
                       14,16) == c("01A")]
# Extract 11A (normal)
counts11A <- counts_annot_protein_coding[,
             substring(colnames(counts_annot_protein_coding),
                       14,16) == c("11A")]



#### tpms ####
tpms <- expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(tpms) <- expquery2@colData@rownames
rownames(tpms) <- expquery2@rowRanges@ranges@NAMES
tpms <- tpms %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]

rownames(tpms) <- NULL
tpms <- tpms %>% column_to_rownames("symbol") 
tpms <- tpms[tpms$type == "protein_coding",]
tpms <- tpms[,-c(1,ncol(tpms))]
colnames(tpms) <- substring(colnames(tpms),1,16)
tpms <- tpms[,!duplicated(colnames(tpms))]
tpms01A <- tpms[,substring(colnames(tpms),14,16) == c("01A")]
tpms11A <- tpms[,substring(colnames(tpms),14,16) == c("11A")]

# Check if the row names and column names are the same in 2 data sets
identical(rownames(counts01A),rownames(counts11A))
identical(rownames(tpms01A),rownames(tpms11A))
identical(rownames(counts01A),rownames(tpms01A))
identical(colnames(counts01A),colnames(tpms01A))
identical(colnames(counts11A),colnames(tpms11A))

# Save counts and tpms
# write.table(counts01A, "counts01A.txt", sep = "\t", 
#             row.names = T, col.names = NA, quote = F)
# write.table(counts11A,"counts11A.txt",sep = "\t",
#             row.names = T, col.names = NA, quote = F)
# write.table(tpms01A, "tpms01A.txt", sep = "\t",
#             row.names = T, col.names = NA, quote = F)
# write.table(tpms11A, "tpms11A.txt", sep = "\t",
#             row.names = T, col.names = NA, quote = F)

# Merge
counts <- cbind(counts01A, counts11A)
tpms <- cbind(tpms01A, tpms11A)

write.table(counts,"counts.txt", sep = "\t",
            row.names = T, col.names = NA, quote = F)
write.table(tpms,"tpms.txt",sep = "\t",
            row.names = T, col.names = NA, quote = F)



# Check the range of the data
range(tpms)
# log2 transform. Why need + 1? 
# log2(0) is undefined, so we add 1 to all values to avoid log2(0).
tpms_log2 <- log2(tpms+1)
range(tpms_log2)
tpms01A_log2 <- log2(tpms01A+1)
range(tpms01A_log2)
tpms11A_log2 <- log2(tpms11A+1)
range(tpms11A_log2)

write.table(tpms_log2,"tpms_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A_log2,"tpms01A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A_log2,"tpms11A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

#  Gene Expression Processed Data