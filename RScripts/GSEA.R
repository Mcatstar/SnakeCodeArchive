library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)
library(glue)
library(msigdbr)
library(data.table)

rm(list = ls())
#
## 必要输入文件
degfile ="datasets/GSE72326_deg.txt"
## 输出文件
GSEAresult = "datasets/GSE72326_gsea.csv"
GSEAplot = "plots/GSE72326_gsea.pdf"
##
#

# GSEA需要LogFC, 需导入差异分析结果或差异分析
DEG = read.table(
  file = degfile,
  sep = "\t",
  check.names = F,
  stringsAsFactors = F,
  header = T,
  row.names = 1
)
msigdbr_species()
# 如果下载失败则运行
# options(timeout = 600)
# unlink(file.path(tempdir(), "msigdbr"), recursive = TRUE) # 删除缓存
hallmark_gene_set <- msigdbr(
  species = "Homo sapiens",
  # 人类
  collection = "H"              # H = Hallmark通路（无subcategory）
)
# 按通路分组，提取Symbol基因名
hallmark_list <- split(hallmark_gene_set$gene_symbol, hallmark_gene_set$gs_name)

# 保存为GMT文件（Symbol版，可直接用于GSEA()函数）
gmt_lines <- lapply(names(hallmark_list), function(pathway) {
  genes <- hallmark_list[[pathway]]
  # 拼接成GMT行：通路名 + 空字符串（GMT第二列固定） + 基因列表
  paste(c(pathway, "", genes), collapse = "\t")
})

# 将所有行写入GMT文件
writeLines(
  text = unlist(gmt_lines),
  con = "datasets/hsa_hallmark_symbol.gmt",
  sep = "\n" # 每行一个通路
) + rm(gmt_lines, hallmark_gene_set, hallmark_list)

# 读取基因和LogFC
geneList = DEG[, "log2FoldChange"]
names(geneList) = as.character(DEG[, "gene_symbol"])
geneList = sort(geneList, decreasing = TRUE) # 排序

# GSEA
gmt = read.gmt("datasets/hsa_hallmark_symbol.gmt")
GSEA = GSEA(
  geneList = geneList,
  # 你的排序基因列表
  TERM2GENE = gmt,
  # 你的GMT基因集（hallmark）
  pvalueCutoff = 0.05,
  # 显著性阈值
  pAdjustMethod = "BH",
  # 多重检验校正方法
  minGSSize = 20,
  # 最小基因集大小（过滤小通路）
  maxGSSize = 500,
  # 最大基因集大小（过滤大通路）
  eps = 0,
  # 精准估算极小P值
  # verbose = FALSE               # 关闭冗余输出
)
res = GSEA@result
write.csv(res, file = GSEAresult)

# 作图
library(patchwork)
p = gseaplot2(GSEA, geneSetID = 1:6)
ggsave(
  plot = p[[1]] / p[[2]] / p[[3]]+
    plot_layout(heights = c(2, 1, 1.5)),  # 高度比例,
  file = GSEAplot,
  width = 12,
  height = 8
)
