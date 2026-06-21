library(STRINGdb)
library(dplyr)
library(igraph)

options(timeout = 300) 
## 必要的输入文件
##
degfile ="/datasets/GSE72326_deg.txt"
## 输出文件
ppi_edges = "/datasets/GSE72326ppi_edges.csv"
ppi_node_features = "/datasets/GSE72326ppi_node_features.csv"
ppi_node_mapping = "/datasets/GSE72326ppi_node_mapping.csv"
##

DEG <- read.table(
  file = degfile,
  sep = "\t",
  header = TRUE,        # 第一行为列名
  row.names = 1,        # 第一列为行名
  quote = "",           # 禁用引号解析
  stringsAsFactors = FALSE  # 避免自动将字符转为因子
)

diff = DEG[DEG$regulate != "Normal", ]

deg_genes = diff[, c("log2FoldChange","padj","gene_symbol")]

gene_list <- deg_genes$gene_symbol

# 物种代码：人类 9606，小鼠 10090，大鼠 10116 等
string_db <- STRINGdb$new(version = "11.5", species = 9606, 
                          score_threshold = 400,  # 中等置信度，可调整
                          input_directory = "/datasets/")
# 将基因名映射到 STRING 内部 ID
# 注意：map() 函数返回的数据框中，`STRING_id` 是内部 ID，`queryItem` 是原始基因名
genes_mapped <- string_db$map(data.frame(gene = gene_list), 
                              "gene", 
                              removeUnmappedRows = TRUE)

subgraph <- string_db$get_subnetwork(genes_mapped$STRING_id)
# 如果 subgraph 为空，可能是因为置信度阈值太高或基因太少，可以降低阈值重试

# igraph 对象的边列表
edge_list <- as_edgelist(subgraph, names = TRUE)
# 转换为数据框，并添加边属性（如结合分数）
# STRINGdb 的 igraph 对象中，边的属性存储在 edge_attr 中
edge_attrs <- edge_attr(subgraph)
if ("combined_score" %in% names(edge_attrs)) {
  edges_df <- data.frame(
    source = edge_list[,1],
    target = edge_list[,2],
    combined_score = edge_attrs$combined_score,
    stringsAsFactors = FALSE
  )
} else {
  edges_df <- data.frame(
    source = edge_list[,1],
    target = edge_list[,2],
    stringsAsFactors = FALSE
  )
}

# 去重（无向图，去除重复边）
edges_df <- distinct(edges_df)

# 6. 构建节点属性表
# 获取子网络中所有节点的 STRING ID 和对应的基因名
# 节点在 igraph 中的名称是 STRING ID，需要映射回基因名
node_ids <- V(subgraph)$name
# 从 genes_mapped 中查找对应的基因名
node_info <- genes_mapped %>%
  filter(STRING_id %in% node_ids) %>%
  select(STRING_id, gene) %>%
  rename(gene = gene,string_id = STRING_id)

# 合并差异表达信息
nodes_df <- node_info %>%
  left_join(deg_genes, by = "gene_symbol") %>%
  # 对于没有差异表达信息的基因（可能不在 deg_data 中），填充 NA 或 0
  mutate(
    logFC = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    pvalue = ifelse(is.na(padj), 1, padj)
  ) %>% 
  select(-log2FoldChange, -padj) %>%
  rename(gene = gene_symbol)

# 添加节点拓扑属性（可选，如度中心性）
degrees <- degree(subgraph)
nodes_df$degree <- degrees[match(node_info$string_id, names(degrees))]

# 7. 保存为 CSV 文件（兼容 Cytoscape 和 GNN）
# 边文件：至少两列（source, target），其他属性作为附加列
write.csv(edges_df, file = ppi_edges, row.names = FALSE, quote = FALSE)

# 节点文件：第一列应与边文件中的节点名对应（基因名）
write.csv(nodes_df, file = ppi_node_features.csv, row.names = FALSE, quote = FALSE)

# 如果希望同时保留 STRING ID 和基因名，可单独保存映射表
write.csv(node_info, file = ppi_node_mapping.csv, row.names = FALSE, quote = FALSE)
