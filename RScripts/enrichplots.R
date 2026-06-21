library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)
#
## 必要输入文件
EnrichObj = "datasets/GSE72326_enrich.rda" # 总体富集对象
## 输出文件
dotplot ="plots/GSE72326_deg.svg"
circosplot = "plots/GSE72326_deg.svg"
##
#

# 作图
load(EnrichObj)
dotgo = enrichplot::dotplot(GO, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = "free")
ggsave(
  glue("{dotplot}_enrichgo.svg"),
  plot = dotgo,
  width = 8,
  height = 12
)

# 作图
dotkegg = enrichplot::dotplot(KEGG, showCategory = 20) + 
  facet_grid(scale = "free")
ggsave(
  glue("{dotplot}_kegg.svg"),
  plot = dotkegg,
  width = 8,
  height = 12
)

# # 富集圈图
#
# 待完成
#