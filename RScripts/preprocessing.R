# ===系统控制，勿动
library(GEOquery)
library(tidyverse)
library(limma)
library(stringr)
library(glue)
set_dir = "datasets"
desc = "GSE43292"
id = "GSE43292"
exp_dir = glue("{set_dir}/{id}_exp.csv")
group_dir = glue("{set_dir}/{id}_group.csv")
cli_dir = glue("{set_dir}/{id}_cli.csv")

# ===
gset = GEOquery::getGEO(id,
                        destdir = set_dir,
                        AnnotGPL = TRUE,
                        getGPL = FALSE)
class(gset)
gset[[1]]

exp = gset[[1]]@assayData$exprs
# 标准化
exp = normalizeBetweenArrays(exp)
boxplot(exp,
        outline = FALSE,
        notch = TRUE,
        las = 2)
# log2变换
range(exp)
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { exp[which(exp <= 0)] = NaN
exp = log2(exp + 1)
print("log2 transform finished")}else{print("log2 transform not needed")}
range(exp)
boxplot(exp,
        outline = FALSE,
        notch = TRUE,
        las = 2)
# 基因注释
gpl_id <- gset[[1]]@annotation # 获取平台
gpl = getGEO(gpl_id, destdir = set_dir) # 下载注释文件
gpl_table = gpl@dataTable@table  # 获取注释矩阵
  # filter(`Species Scientific Name` == "Homo sapiens") %>%  # 筛选人类基因
ids = gpl_table %>% 
  # filter(`Species Scientific Name` == "Homo sapiens") %>%  # 筛选人类基因
  select(ID, gene_assignment)%>% # 截取探针symbol对照
  rename(probe_id = ID, symbol = gene_assignment) %>% # 列名重命名
  # 切割多基因符号并去空格
  mutate(symbol = trimws(str_split(symbol, "//", simplify = TRUE)[, 2])) %>%
  # 删除空白symbol仅保留表达矩阵中存在的探针（对应原筛选行）
  filter(symbol != "") %>%
  filter(probe_id %in% rownames(exp))
ids$probe_id <- as.character(ids$probe_id) %>% trimws()
rownames(exp) <- as.character(rownames(exp)) %>% trimws()
exp = exp[ids$probe_id, ] # 只保留有对应symbol的探针
table(rownames(exp) == ids$probe_id)
exp2 = cbind(ids, exp) # 合并为大矩阵
exp2$symbol = as.character(exp2$symbol)
# 基因去重
exp3 = exp2 %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames(var = "symbol")

# 临床数据
rt1 = pData(gset[[1]])
write.csv(rt1,
          file = cli_dir,
          row.names = TRUE)
# 疾病分组
colnames(rt1)
table(rt1$`disease state:ch1`)
group_list = ifelse(str_detect(rt1$`characteristics_ch1`, "Healthy"), "Normal", "Diseased")
rt2 = rt1 %>% 
  mutate(group = group_list) %>% # 添加分组
  dplyr::select("group") %>% arrange(desc(group)) # 分组排列

# 数据保存
exp3 = exp3[, rownames(rt2)] # 行列对应
identical(rownames(rt2), colnames(exp3))
group = rt2
expr_matrix = exp3

write.csv(group,
          file = group_dir,
          row.names = TRUE)
write.csv(expr_matrix,
          file = exp_dir,
          row.names = TRUE)
# rm(list=ls())

