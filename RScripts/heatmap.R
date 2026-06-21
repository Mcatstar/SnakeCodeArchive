rm(list = ls())
#
## 必要输入文件
expfile="datasets/GSE72326_exp.csv" # 表达矩阵
group ="datasets/GSE72326_group.csv" # 分组信息
## 输出文件
heatplot = "plots/GSE72326_heatmap.png"
##
#

# 热图
DEG <- read.csv(expfile, row.names = 1)
library(pheatmap)
diff = DEG[DEG$regulate != "Normal", ]
exp_diff = exp[diff$gene_symbol, ] # 全部展示
range(exp_diff)
hmap_plot = pheatmap::pheatmap(
  exp_diff,
  annotation = group,
  scale = "row",
  # 对每一行进行 Z-score 标准化
  cluster_rows = TRUE,
  # 行聚类
  cluster_cols = FALSE,
  # 列聚类
  show_rownames = TRUE,
  show_colnames = FALSE,
  # 剪枝聚类树（控制深度）
  # cutree_rows = 5,# 行聚类树剪为5个簇（减少分支）
  # cutree_cols = 3,# 列聚类树剪为3个簇
  # 给剪枝后的簇添加间隙
  # gaps_row = cutree(hclust(dist(t(scale(t(gene_for_rich))))), 5),
  # gaps_col = cutree(hclust(dist(scale(t(gene_for_rich)))), 3),
  # main = "Differentially Expressed Genes Heatmap", # 主标题
  # treeheight_row = 10,   # 行聚类树高度（默认是50，调小到5~10）
  # treeheight_col = 10,   # 列聚类树宽度（同理调小）
  color = colorRampPalette(c("#0066CC", "white", "#FF6600"))(200),
  border_color = NA,
  # breaks = seq(-2, 2, length.out = 200),
  filename = heatplot
)
# 特定展示
topnum = 15 # 展示上下调各15个
uptop = diff %>% top_n(topnum, log2FoldChange)
downtop = diff %>% top_n(-topnum, log2FoldChange)
difftop = c(uptop$gene_symbol, downtop$gene_symbol)
exp_difftop = exp[difftop, ]
hmap_plot = pheatmap::pheatmap(
  exp_difftop,
  annotation = group,
  scale = "row",
  # 对每一行进行 Z-score 标准化
  cluster_rows = TRUE,
  # 行聚类
  cluster_cols = FALSE,
  # 列聚类
  show_rownames = TRUE,
  show_colnames = FALSE,
  # 剪枝聚类树（控制深度）
  # cutree_rows = 5,# 行聚类树剪为5个簇（减少分支）
  # cutree_cols = 3,# 列聚类树剪为3个簇
  # 给剪枝后的簇添加间隙
  # gaps_row = cutree(hclust(dist(t(scale(t(gene_for_rich))))), 5),
  # gaps_col = cutree(hclust(dist(scale(t(gene_for_rich)))), 3),
  # main = "Differentially Expressed Genes Heatmap", # 主标题
  # treeheight_row = 10,   # 行聚类树高度（默认是50，调小到5~10）
  # treeheight_col = 10,   # 列聚类树宽度（同理调小）
  color = colorRampPalette(c("#0066CC", "white", "#FF6600"))(200),
  border_color = NA,
  # breaks = seq(-2, 2, length.out = 200),
  filename = glue("{heatplot}topheatmap.pdf")
)
