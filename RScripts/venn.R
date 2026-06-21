rm(list = ls())

# 设定基本信息
#
## 必要输入文件
exp1="datasets/GSE72326_exp.csv" # 表达矩阵
exp2 ="datasets/GSE72326_exp.csv" # 分组信息
## 输出文件
upplot ="plots/GSE72326_up.png"
downplot = "plots/GSE72326_down.png"
# 包含字段：
# [1] "log2FoldChange" "AveExpr"        "t"              "P.Value"       
# [5] "padj"           "B"              "regulate"       "gene_symbol"   
##
#

# 读取差异分析结果
DEG1 = read.csv(exp1, row.names = 1)
DEG2 = read.csv(exp2, row.names = 1)

# 截取上/下调基因列表
up1 = DEG1[DEG1$regulate == "Up", ]$gene_symbol
up2 = DEG2[DEG2$regulate == "Up", ]$gene_symbol
down1 = DEG1[DEG1$regulate == "Down", ]$gene_symbol
down2 = DEG2[DEG2$regulate == "Down", ]$gene_symbol

# 绘图
library(VennDiagram)
# 全局关闭VennDiagram的消息输出
options(verbose = FALSE)
x_up = list(A = up1, B = up2)
x_down = list(A = down1, B = down2)
venn_up = VennDiagram::venn.diagram(
  x_up,
  filename = NULL,
  # 设为NULL，不直接保存，返回绘图对象
  category.names = c(glue("{gse_one}_Up"), glue("{gse_two}_Up")),
  # 集合标签
  scaled = FALSE,
  # 禁用按数量缩放，两圆大小一致
  ext.text = FALSE,
  # 配合scaled=FALSE，避免标签位置偏移
  # 离散配色核心参数（按集合指定）
  fill = c(
    rgb(178, 202, 246, maxColorValue = 255),
    # 第一个集合（左圆）：淡蓝
    rgb(223, 159, 255, maxColorValue = 255)   # 第二个集合（右圆）：淡紫
  ),
  alpha = 0.7,
  # 透明度
  cat.pos = c(0, 0),
  # 集合标签水平居中
  cat.dist = 0.08,
  # 标签与圆的距离
  margin = 0.15,
  # 图形边缘留白
  sep.dist = 0.02,
  # 两圆间距
  label.col = "black",
  # 数字标签黑色
  cat.col = "black",
  # 集合标签黑色
  cex = 1.2,
  # 数字大小
  cat.cex = 0.9,
  # 集合标签大小
  fontfamily = "Arial",
  # 字体
) +
  ggsave(upplot, venn_up)
venn_down = VennDiagram::venn.diagram(
  x_down,
  filename = NULL,
  category.names = c(glue("{gse_one}_Down"), glue("{gse_two}_Down")),
  scaled = FALSE,
  # 禁用按数量缩放，两圆大小一致
  ext.text = FALSE,
  # 配合scaled=FALSE，避免标签位置偏移
  # 离散配色核心参数（按集合指定）
  fill = c(
    rgb(178, 202, 246, maxColorValue = 255),
    # 第一个集合（左圆）：淡蓝
    rgb(223, 159, 255, maxColorValue = 255)   # 第二个集合（右圆）：淡紫
  ),
  alpha = 0.7,
  # 透明度
  cat.pos = c(0, 0),
  # 集合标签水平居中
  cat.dist = 0.08,
  # 标签与圆的距离
  margin = 0.15,
  # 图形边缘留白
  sep.dist = 0.02,
  # 两圆间距
  label.col = "black",
  # 数字标签黑色
  cat.col = "black",
  # 集合标签黑色
  cex = 1.2,
  # 数字大小
  cat.cex = 0.9,
  # 集合标签大小
  fontfamily = "Arial",
  # 字体
) +
  ggsave(downplot, venn_down)
