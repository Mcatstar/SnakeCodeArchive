library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)

rm(list = ls())
#
## 必要输入文件
degfile ="datasets/GSE72326_deg.txt" # 应包含基因名
# 包含字段：
# [1] "log2FoldChange" "AveExpr"        "t"              "P.Value"       
# [5] "padj"           "B"              "regulate"       "gene_symbol"
# 其中"gene_symbol"为核心字段，不可缺失
## 输出文件
# IDtrans = "datasets/GSE72326_go.txt"
Go_result ="datasets/GSE72326_go.txt" # Go结果
KEGG_result ="datasets/GSE72326_kegg.txt" # KEGG结果
EnrichObj = "datasets/GSE72326_enrich.rda" # 总体富集对象
##
#

DEG = read.table(
  file = degfile,
  sep = "\t",
  check.names = F,
  stringsAsFactors = F,
  header = T,
  row.names = 1
)

diff = DEG[DEG$regulate != "Normal", ]$gene_symbol

# 转化ID
gene_entrez = bitr(diff,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
untrans = setdiff(diff, gene_entrez$SYMBOL) # 未成功转化的ID
head(untrans)

library(AnnotationDbi)
entrez_ids = AnnotationDbi::mapIds(
  x = org.Hs.eg.db,
  keys = untrans,
  keytype = "ALIAS",
  column = "ENTREZID",
  multiVals = "first"
)
entrez_df = data.frame(ENTREZID = entrez_ids, check.names = FALSE)
entrez_df = na.omit(entrez_df)
entrez_df = entrez_df %>% rownames_to_column(var = "SYMBOL")
# 合并
gene_entrez = rbind(gene_entrez, entrez_df)


# GO富集
GO = enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  # BP MF CC ALL
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
Go_res = GO@result
write.table(
  Go_res,
  file = Go_result,
  sep = "\t",
  quote = F,
  row.names = F
)

# KEGG富集
KEGG = enrichKEGG(gene=gene_entrez$ENTREZID,
                  organism = "hsa",
                  pvalueCutoff = 0.5,
                  qvalueCutoff = 0.5)
KEGG_res = KEGG@result
write.table(
  KEGG_res,
  file = KEGG_result,
  sep = "\t",
  quote = F,
  row.names = F
)
# 保存富集对象
save(GO, KEGG, file = EnrichObj)
