####cibersort####

# setwd("CIBERSORT")   
install.packages('e1071')
install.packages('parallel')
BiocManager::install("preprocessCore", version = "3.19")

install.packages("ggstatsplot")
install.packages("ggcorrplot")
install.packages("corrplot")
library(ggstatsplot)
library(ggcorrplot)
library(corrplot)

library(e1071)
library(parallel)
library(preprocessCore)
library(tidyverse)
source("CIBERSORT/CIBERSORT.R")   # The proportion of immune cells in the tumor microenvironment can be estimated using CIBERSORT

sig_matrix <- "CIBERSORT/LM22.txt"   # LM22.txt is a gene expression profile of 22 immune cells provided by CIBERSORT
mixture_file = 'TCGAdata/tpms01A_log2.txt'   

res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
res_cibersort <- res_cibersort[,1:22]   
# Removal of cells with all zero abundance
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   
ciber.res <- as.data.frame(ciber.res)
write.table(ciber.res,"CIBERSORT/ciber_new.res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

ciber.res <- read.table("CIBERSORT/ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#### cibersort 
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   


#分组比较图
a <- ciber.res
#读取肿瘤患者01A表达谱
exp <- read.table("TCGAdata/tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
med=median(as.numeric(exp["BTK",]))
exp <- exp %>% t() %>% as.data.frame()
exp <- exp %>% mutate(group=factor(ifelse(exp$BTK>med,"high","low"),levels = c("low","high")))
class(exp$group)
identical(rownames(a),rownames(exp))
a$group <- exp$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=CIBERSORT,value = Fraction,-c(group,sample))
ggboxplot(b, x = "CIBERSORT", y = "Fraction",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()


####相关性热图####

cor<-sapply(ciber.res,function(x,y) cor(x,y,method="spearman"),ciber.res)
rownames(cor)<-colnames(ciber.res)

ggcorrplot(cor, 
           hc.order = TRUE, #使用hc.order进行排序
           type = "upper", #图片位置
           outline.color = "white",#轮廓颜色
           lab = TRUE,#true为在图上添加相关系数
           ggtheme = ggplot2::theme_gray, #指ggplot2函数对象，默认值为thememinimal
           colors = c("#01468b", "white", "#ee0000"))

####基因与cibersort相关性散点图####
exp = read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- exp["BTK",]
ciber = read.table("ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber <- ciber %>% t() %>% as.data.frame()
rownames(ciber) <- gsub(" ",".",rownames(ciber))
identical(colnames(ciber),colnames(exp))
exp_ciber <- rbind(exp,ciber)
exp_ciber <- exp_ciber %>% t() %>% as.data.frame()
#install.packages("ggside")
library(ggstatsplot)
library(ggside)
ggscatterstats(data = exp_ciber, #要分析的数据
               y = BTK, #设置Y轴
               x = B.cells.naive,#设置X轴
               type = "nonparametric", 
               margins = "both",#是否显示 边缘，默认为true                                      
               xfill = "#01468b", #x轴边缘图形的颜色
               yfill = "#ee0000", #y轴边缘图形的颜色
               marginal.type = "densigram")#在图片坐标轴边缘添加图形类型

####韦恩图####