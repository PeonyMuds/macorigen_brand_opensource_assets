# 安装R包
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")


# 加载R包
library(DESeq2)
# 读取count数据表
countsTable <- read.delim("result/5.DESeq2/counts.xls",header = T,row.names = 1)
# 读取样本分组信息
countsDesign <- read.csv("group.csv",sep="\t",header=TRUE)
# 先将读取数据由data.frame转换成matrix类型
countsMatrix <- as.matrix(countsTable)
# 再将countsMatrix转换成DESeq2的数据格式（构建dds）：
dds <- DESeqDataSetFromMatrix(countsMatrix,colData = countsDesign, design = ~ group)
# 将所有样本中count值均为0的基因去除
dds <- dds[rowSums(counts(dds))>0,]
# DESeq差异分析：大小因子的估计，离差的估计，负二项分布的拟合以及计算相应的统计量
dds <- DESeq(dds,parallel = T) 
# 提取组间差异表达基因的分析结果：
group_res <- results(dds, contrast = c("group","A","B"), parallel=T)
# 提取组间显著差异表达的基因（padj<0.01）：
group_Sig_0.01 <- group_res[which(group_res$padj<0.01),]
# 将分析结果写入文件：
write.table(group_Sig_0.01,file = "group_Sig_0.01_DESeq2.xls",sep="\t",quote=FALSE)
