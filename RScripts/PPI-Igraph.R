rm(list = ls())

## 必要的输入文件
## 输入文件
net_edges = "datasets/GSE72326ppi_edges.csv"
net_node_features = "datasets/GSE72326ppi_node_features.csv"
net_node_mapping = "datasets/GSE72326ppi_node_mapping.csv"

##

# 1.节点注释
library(tidyverse)

edges = read_csv(net_edges)
nodes = read_csv(net_node_features)
mapfile = read_csv(net_node_mapping)

edges_table <- edges %>%
  # 1. 匹配 Source
  left_join(
    mapfile %>% select(string_id, gene), 
    by = c("source" = "string_id"),
    relationship = "many-to-many"
  ) %>% 
  rename(source_gene = gene) %>%
  # 2. 匹配 Target
  left_join(
    mapfile %>% select(string_id, gene), 
    by = c("target" = "string_id"),
    relationship = "many-to-many"
  ) %>% 
  rename(target_gene = gene) %>%
  # 防止后续画图因为 NA 报错
  filter(!is.na(source_gene) & !is.na(target_gene)) %>%
  # 防止因为一个ID对应多个别名导致边重复
  distinct() %>% 
  filter(combined_score >= 700) %>% 
  select(source_gene, target_gene, everything(), -source, -target) %>% 
  rename(source = source_gene, target = target_gene)

all_genes = c(edges_table$source, edges_table$target) %>% unique()
nodes_table <- nodes %>%
  select(-string_id) %>% 
  rename(symbol = gene) %>% 
  filter(symbol %in% all_genes)
write.csv(edges_table,
            file = "datasets/GSE72326ppi_edges_anno.csv",
            row.names = FALSE,   # 关闭行号
            quote = FALSE)       # 关闭引号

write.csv(nodes_table,
          file = "datasets/GSE72326ppi_node_anno_features.csv",
          row.names = FALSE,   # 关闭行号
          quote = FALSE)       # 关闭引号
