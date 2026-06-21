rm(list = ls())

#
## 必要输入文件
expfile="datasets/GSE72326_exp.csv" # 表达矩阵
group ="datasets/GSE72326_group.csv" # 分组信息
## 输出文件
dataInput = "datasets/GSE72326_wgcna_dataInput.rda"

#=====================================================================================
#
#  Code chunk 1
# 数据预处理
#=====================================================================================

library(tidyverse)
library(WGCNA)

## 1.读取数据
exp_matrix = read.csv(expfile, row.names = 1)
group = read.csv(group, row.names = 1)

exp = exp_matrix %>% 
  mutate(across(everything(), as.numeric)) %>% # 数据类型转换
  t() %>%  # 矩阵转置,适配后续分析
  as.data.frame()
identical(rownames(group), rownames(exp))
# 高变基因
var_res = apply(exp, 2, var)
per_res = quantile(var_res, probs = seq(0, 1, 0.25))
per_genes = exp[, which(var_res > per_res[4])] # 4:前25%，3：前50%，2：前70%
exp = data.matrix(per_genes)
## 2.数据检查
gsg = goodSamplesGenes(exp, verbose = 3)
if(gsg$allOK){
  print("gsg$allOK is True")
} else {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(exp)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(exp)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exp = exp[gsg$goodSamples, gsg$goodGenes]
}

# 样本聚类
sampleTree = hclust(dist(exp), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#切割末端聚类样本
abline(h = 76, col = "red") # 展示切割线，可适当调整
clust = cutreeStatic(sampleTree, cutHeight = 120, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = exp[keepSamples, ]

# 表型数据
table(group)
traitData = group %>% 
  mutate("Normal" = ifelse(group=="Normal", 1, 0),
         "Diseased" = ifelse(group=="Diseased", 1, 0)) %>% 
  select(-group) %>% 
  as.data.frame()

sameSample = intersect(rownames(datExpr), rownames(traitData))
datExpr = datExpr[sameSample, ]
datTraits = traitData[sameSample, ]

# 再聚类
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = dataInput)
## 
# 垃圾清理
rm(list = ls())
collectGarbage()
##


