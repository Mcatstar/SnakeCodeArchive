# SnakeCodeArchive

基于 GEO 数据库的转录组差异表达分析与功能挖掘流程，涵盖差异基因筛选、功能富集、WGCNA 共表达网络、STRING PPI 网络以及基于图卷积网络（GCN）的蛋白互作链路预测。

## 项目结构

```
SnakeCodeArchive/
├── RScripts/               # R 分析流程
├── PyScripts/              # Python 分析脚本
│   └── prolinks.py         # GCN 链路预测
├── datasets/               # 输入数据与中间文件（gitignored）
├── results/                # 最终结果（gitignored）
│   ├── ppi_links_pred.pth  # 训练好的 GCN 模型权重
│   └── predicted_interactions.csv  # Top 100 预测的新 PPI
├── plots/                  # 生成的可视化图表（gitignored）
├── demo.R                  # 示例脚本
└── SnakeCodeArchive.Rproj  # RStudio 项目文件
```

## 分析流程

### 1. 数据获取与预处理（R）
- `preprocessing.R` — 从 GEO 下载 GSE series，进行归一化、log2 转换、探针注释与去重

### 2. 差异表达分析（R）
- `limma.R` — 使用 limma 识别差异表达基因（DEGs），筛选标准：\|logFC\| > 0.5，adj.P.Val < 0.05

### 3. 差异基因可视化（R）
- `volcano.R` — 火山图，标注 Top 10 显著基因
- `venn.R` — 韦恩图比较差异基因重叠
- `heatmap.R` — 差异基因表达热图

### 4. 功能富集分析（R）
- `enrich.R` — GO（BP/MF/CC）和 KEGG 通路富集
- `enrichplots.R` — 富集结果气泡图

### 5. 基因集富集分析（R）
- `GSEA.R` — 基于 MSigDB Hallmark 基因集的 GSEA

### 6. 加权基因共表达网络（R）
- `WGCNA-DataInput.R` — 筛选高变异基因，构建表达矩阵与性状矩阵
- `WGCNA-networkConstr-auto.R` — 自动构建共表达网络
- `WGCNA-networkConstr-man.R` — 手动逐步构建（备选）
- `WGCNA-relateModsToExt.R` — 模块-性状相关性分析
- `WGCNA-Visualization.R` — TOM 热图与特征基因网络可视化

### 7. 蛋白互作网络（R）
- `PPI-Network.R` — 通过 STRINGdb 获取 PPI 子网络
- `PPI-Igraph.R` — PPI 网络注释与清洗（score ≥ 700）

### 8. GCN 链路预测（Python）
- `prolinks.py` — 构建 GCN 模型（2 层 GCN + 内积解码器），预测缺失的蛋白互作关系

## 使用

根据需求运行 `RScripts/` 中的 R 脚本和 Python 脚本，脚本正在不断完善

```r
# 在 R 中运行
source("RScripts/preprocessing.R")
source("RScripts/limma.R")
# ...
```

```bash
# Python 运行
python PyScripts/prolinks.py
```

## 依赖

- **R**: Bioconductor（GEOquery, limma, clusterProfiler, WGCNA, STRINGdb），ggplot2 等
- **Python**: PyTorch, pandas, numpy, scikit-learn
