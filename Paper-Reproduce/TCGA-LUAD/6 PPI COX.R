#install.packages("survival")
# install.packages("forestplot")
# library(ggplot2)
library(tidyverse)
library(survival)
library(forestplot)

# DEG_final <- read_csv("PPI/DEG_final.csv")

#### PPI ####
#String：https://cn.string-db.org/cgi/input?sessionId=bXmYsv7CnUrH&input_page_active_form=multiple_identifiers

#### The top 30 genes ####

# Step 1: Read the data (assuming it is a CSV or TSV file)
# Replace 'your_file_path.tsv' with the actual path to your file
data <- read.csv("PPI/Nodes.csv", sep = ",", stringsAsFactors = FALSE)

# Step 2: Sort the data by the 'Degree' column in descending order
# Assuming 'degree' is the column that stores the degree values
sorted_data <- data %>% arrange(desc(Degree))

# Step 3: Select the first 30 genes (or rows)
top_30 <- head(sorted_data, 30)

# Step 4: Plot the bar chart using ggplot2
ggplot(top_30, aes(x = reorder(name, Degree), y = Degree)) +
  geom_bar(stat = "identity", fill = "salmon") +
  coord_flip() +  # Flip the coordinates for a horizontal bar chart
  labs(x = "Gene", y = "Degree", title = "Top 30 Genes by Degree") +
  theme_minimal()

dev.off()

#### COX ####

exp_surv_01A <- read.table("Survival_data/exp_surv_01A.txt", sep = "\t", row.names = 1,
                          check.names = F, stringsAsFactors = F, header = T)
# Read DEG_final.txt 
DEG_final <- read.table("COX/DEG_final.txt", header = TRUE, stringsAsFactors = FALSE)

surv.expr <- cbind(exp_surv_01A[, 1:2], exp_surv_01A[, DEG_final$SYMBOL])
# a <- exp_surv_01A[,1:2]
# b <- exp_surv_01A[,DEG_final$SYMBOL]
 
# Cox analysis
# How to change the column names of a data frame in R?
# colnames(surv.expr)[ ] <- ""  # the index of the column you want to change

Coxoutput <- NULL 

for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time, OS) ~ surv.expr[,i], data = surv.expr) # OS.time, OS are the column names of the survival data
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}

Coxoutput <- arrange(Coxoutput, pvalue)

### Select significant genes
gene_sig <- Coxoutput[Coxoutput$pvalue < 0.005,] 
write.csv(gene_sig, file = "COX/gene_sig.csv")

# Read the top 30 genes
topgene <- read.csv("COX/gene_sig.csv", header = TRUE, stringsAsFactors = FALSE)

# 3. Plot the forest plot
## 3.1  Create the table text
tabletext <- cbind(c("Gene", topgene$gene),
                   c("HR", format(round(as.numeric(topgene$HR), 3), nsmall = 3)),
                   c("lower 95%CI", format(round(as.numeric(topgene$lower),3), nsmall = 3)),
                   c("upper 95%CI", format(round(as.numeric(topgene$upper),3), nsmall = 3)),
                   c("pvalue", format(round(as.numeric(topgene$p), 3), nsmall = 3)))
## 3.2 
forestplot(labeltext = tabletext,
           mean = c(NA, as.numeric(topgene$HR)),
           lower = c(NA, as.numeric(topgene$lower)), 
           upper = c(NA, as.numeric(topgene$upper)),
           graph.pos = 5, # Position of the point estimate
           graphwidth = unit(0.25, "npc"), # Graph width
           fn.ci_norm="fpDrawDiamondCI", # CI shape function
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"), # Colors
           boxsize=0.4, # Box size
           lwd.ci=1,
           ci.vertices.height = 0.1, ci.vertices=T, # CI settings
           zero=1, # Reference line
           lwd.zero=1.5, # Reference line width
           xticks = c(0.5,1,1.5), # X-axis ticks
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2), # Font size
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # Horizontal lines
                           "2" = gpar(lwd=1.5, col="black"),
                           "53" = gpar(lwd=2, col="black")),
           lineheight = unit(0.75,"cm"), # Line height
           colgap = unit(0.3,"cm"), # Column gap
           mar=unit(rep(1.5, times = 4), "cm"), # Margins
           new_page = FALSE # Keep on the same page
)
# 30*30
dev.off()
