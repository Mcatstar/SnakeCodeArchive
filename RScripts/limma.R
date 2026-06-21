rm(list = ls())

#
## 必要输入文件
expfile="../datasets/GSE72326_exp.csv" # 表达矩阵
group ="../datasets/GSE72326_group.csv" # 分组信息
## 输出文件
degfile ="../datasets/GSE72326_deg.txt"
# 包含字段：
# [1] "log2FoldChange" "AveExpr"        "t"              "P.Value"       
# [5] "padj"           "B"              "regulate"       "gene_symbol"   
##
#

library(tidyverse)
library(limma)
library(pheatmap)
library(ggpubr)
library(glue)

exp_matrix = read.csv(expfile, row.names = 1)
group = read.csv(group, row.names = 1)

exp = exp_matrix %>% mutate(across(everything(), as.numeric)) # 数据类型转换
identical(rownames(group), colnames(exp))

# 因子
group_list = factor(group$group, levels = c("Normal", "Diseased"))
# 比较矩阵
design = model.matrix( ~ group_list)
# 线性模型拟合
fit = lmFit(exp, design)
# 贝叶斯检验
fit2 = eBayes(fit)
deg = topTable(fit2, coef = 2, number = Inf)

# 注释上调下调
library(ggVolcano)
### 设置logFC adj.P.val
logFC_cutoff = 0.5
adj_P_Val = 0.05

DEG = ggVolcano::add_regulate(
  deg,
  log2FC_name = "logFC",
  fdr_name = "adj.P.Val",
  log2FC = logFC_cutoff,
  fdr = adj_P_Val
)
DEG = na.omit(Deg)
DEG$gene_symbol <- rownames(DEG)
rownames(DEG) <- NULL# 清空行名，避免后续混淆
write.table(
  DEG,
  file = degfile,
  sep = "\t",
  row.names = T,
  col.names = NA,
  quote = FALSE
)
table(DEG$regulate)
