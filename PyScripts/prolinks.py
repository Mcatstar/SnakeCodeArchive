#!.venv\Scripts\python
# -*- coding: utf-8 -*-
"""
基于 GCN 的 PPI 网络链接预测(链路预测)
输入:ppi_edges.csv, node_features.csv
输出:
  - 模型训练日志、评估指标
  - 预测的潜在相互作用列表(top_k 高分对)
"""

import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import train_test_split
import random
import warnings
warnings.filterwarnings('ignore')


# 1. 参数设置

SEED = 42
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
HIDDEN_DIM = 128
EMBEDDING_DIM = 64
EPOCHS = 200
LR = 0.01
TEST_SIZE = 0.2
VAL_SIZE = 0.1
NEG_SAMPLING_RATIO = 1.0   # 负采样比例（正:负 = 1:1）
TOP_K_PREDICT = 100        # 输出前 K 个潜在相互作用

# 设置随机种子
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)


# 2. 加载数据并构建图

print("Loading data...")
edges_df = pd.read_csv("datasets\GSE72326ppi_edges.csv")
nodes_df = pd.read_csv("datasets\GSE72326ppi_node_features.csv")
mapping = pd.read_csv("datasets\GSE72326ppi_node_mapping.csv")

# 使用 STRING ID 作为节点标识符
genes = nodes_df['string_id'].astype(str).tolist()
gene2idx = {gene: i for i, gene in enumerate(genes)}
idx2gene = {i: gene for gene, i in gene2idx.items()}   # 可选，用于输出时转换

# 构建边索引（无向图，双向边）
edge_list = []
for _, row in edges_df.iterrows():
    u = gene2idx[row['source']]
    v = gene2idx[row['target']]
    edge_list.append([u, v])
    edge_list.append([v, u])
edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()

# 节点特征（使用您需要的列）
feature_cols = ['logFC', 'pvalue', 'degree']
for col in feature_cols:
    if col not in nodes_df.columns:
        nodes_df[col] = 0.0
x = torch.tensor(nodes_df[feature_cols].values, dtype=torch.float)

# 标准化特征
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
x_np = scaler.fit_transform(x.numpy())
x = torch.tensor(x_np, dtype=torch.float)

data = Data(x=x, edge_index=edge_index).to(DEVICE)
print(f"Nodes: {data.num_nodes}, Edges: {data.num_edges // 2} (undirected)")


# 3. 划分训练/验证/测试边

def split_edges(edge_index, test_ratio=0.2, val_ratio=0.1):
    """将正边划分为训练、验证、测试集，并生成负边"""
    # 获取所有无向边（去重）
    edges = edge_index.t().cpu().numpy()
    # 只保留 u < v 的无向边
    undirected_edges = []
    for u, v in edges:
        if u < v:
            undirected_edges.append([u, v])
        else:
            undirected_edges.append([v, u])
    undirected_edges = np.array(undirected_edges)
    undirected_edges = np.unique(undirected_edges, axis=0)  # 去重
    n_edges = len(undirected_edges)
    
    # 随机划分
    indices = np.arange(n_edges)
    np.random.shuffle(indices)
    test_start = int(n_edges * (1 - test_ratio))
    val_start = int(n_edges * (1 - test_ratio - val_ratio))
    train_indices = indices[:val_start]
    val_indices = indices[val_start:test_start]
    test_indices = indices[test_start:]
    
    train_edges = undirected_edges[train_indices]
    val_edges = undirected_edges[val_indices]
    test_edges = undirected_edges[test_indices]
    
    # 构建训练、验证、测试的边索引（包含双向）
    def to_full_edges(edges_undirected):
        full = []
        for u, v in edges_undirected:
            full.append([u, v])
            full.append([v, u])
        return torch.tensor(full, dtype=torch.long).t().contiguous()
    
    train_edge_index = to_full_edges(train_edges)
    val_edge_index = to_full_edges(val_edges)
    test_edge_index = to_full_edges(test_edges)
    
    return train_edge_index, val_edge_index, test_edge_index, train_edges, val_edges, test_edges

# 划分正边
train_edge_index, val_edge_index, test_edge_index, train_edges, val_edges, test_edges = split_edges(
    data.edge_index, test_ratio=TEST_SIZE, val_ratio=VAL_SIZE
)
print(f"Train edges: {len(train_edges)}, Val edges: {len(val_edges)}, Test edges: {len(test_edges)}")


# 4. 负采样（在训练/验证/测试时动态采样）

def sample_negative_edges(positive_edges, num_nodes, num_neg_samples):
    """从非连接节点对中随机采样负边（无向）"""
    # 将正边转换为集合用于快速查找
    pos_set = set()
    for u, v in positive_edges:
        if u < v:
            pos_set.add((u, v))
        else:
            pos_set.add((v, u))
    negatives = []
    while len(negatives) < num_neg_samples:
        u = random.randint(0, num_nodes - 1)
        v = random.randint(0, num_nodes - 1)
        if u == v:
            continue
        if u > v:
            u, v = v, u
        if (u, v) not in pos_set:
            negatives.append([u, v])
            pos_set.add((u, v))  # 避免重复采样
    return torch.tensor(negatives, dtype=torch.long)


# 5. 定义图神经网络模型

class GCNEncoder(nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, out_channels)
    
    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        return x

class LinkPredictor(nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.encoder = GCNEncoder(in_channels, hidden_channels, out_channels)
        self.decoder = lambda z, edge_index: (z[edge_index[0]] * z[edge_index[1]]).sum(dim=-1)
    
    def forward(self, x, edge_index):
        z = self.encoder(x, edge_index)
        return z
    
    def decode(self, z, edge_index):
        return self.decoder(z, edge_index)

model = LinkPredictor(in_channels=data.x.size(1), hidden_channels=HIDDEN_DIM, out_channels=EMBEDDING_DIM).to(DEVICE)
optimizer = torch.optim.Adam(model.parameters(), lr=LR)


# 6. 训练函数

def train():
    model.train()
    optimizer.zero_grad()
    
    # 使用训练边进行编码
    z = model(data.x, train_edge_index)
    
    # 正边得分
    pos_edge_index = train_edge_index
    pos_score = model.decode(z, pos_edge_index)
    pos_loss = -F.logsigmoid(pos_score).mean()
    
    # 负采样（数量与正边相同）
    num_neg = pos_edge_index.size(1) // 2  # 因为正边包含双向，所以实际正边对数为 size(1)/2
    neg_edges = sample_negative_edges(train_edges, data.num_nodes, num_neg)
    neg_edge_index = torch.cat([neg_edges, neg_edges.flip(1)], dim=0).t().contiguous()  # 双向
    neg_edge_index = neg_edge_index.to(DEVICE)
    neg_score = model.decode(z, neg_edge_index)
    neg_loss = -F.logsigmoid(-neg_score).mean()
    
    loss = pos_loss + neg_loss
    loss.backward()
    optimizer.step()
    
    return loss.item()


# 7. 评估函数（AUC, AP）

@torch.no_grad()
def evaluate(positive_edges, negative_edges):
    model.eval()
    z = model(data.x, data.edge_index)  # 使用所有边编码，但实际预测时应该用全图
    
    # 正边得分：将 numpy 数组转为 tensor，构造双向边
    pos_edges_tensor = torch.tensor(positive_edges, dtype=torch.long, device=DEVICE)
    pos_edges_full = torch.cat([pos_edges_tensor, pos_edges_tensor.flip(1)], dim=0).t().contiguous()
    pos_scores = model.decode(z, pos_edges_full).cpu().numpy()
    
    # 负边得分
    neg_edges_tensor = torch.tensor(negative_edges, dtype=torch.long, device=DEVICE)
    neg_edges_full = torch.cat([neg_edges_tensor, neg_edges_tensor.flip(1)], dim=0).t().contiguous()
    neg_scores = model.decode(z, neg_edges_full).cpu().numpy()
    
    scores = np.concatenate([pos_scores, neg_scores])
    labels = np.concatenate([np.ones_like(pos_scores), np.zeros_like(neg_scores)])
    
    auc = roc_auc_score(labels, scores)
    ap = average_precision_score(labels, scores)
    return auc, ap

# 8. 训练循环

print("Training...")
best_val_auc = 0.0
for epoch in range(1, EPOCHS + 1):
    loss = train()
    # 验证集评估（采样负边）
    val_neg = sample_negative_edges(val_edges, data.num_nodes, len(val_edges))
    val_auc, val_ap = evaluate(val_edges, val_neg)
    if epoch % 20 == 0:
        print(f"Epoch {epoch:03d}, Loss: {loss:.4f}, Val AUC: {val_auc:.4f}, Val AP: {val_ap:.4f}")
    if val_auc > best_val_auc:
        best_val_auc = val_auc
        torch.save(model.state_dict(), "results/ppi_links_pred.pth")
        print(f"  -> New best model saved (AUC: {val_auc:.4f})")


# 9. 测试集评估

print("\nLoading best model for test...")
model.load_state_dict(torch.load("results/ppi_links_pred.pth"))   # 模型路径
test_neg = sample_negative_edges(test_edges, data.num_nodes, len(test_edges))
test_auc, test_ap = evaluate(test_edges, test_neg)
print(f"Test AUC: {test_auc:.4f}, Test AP: {test_ap:.4f}")


# 10. 预测潜在相互作用（所有未连接节点对）

print("\nPredicting missing links...")
model.eval()
z = model(data.x, data.edge_index).cpu().detach().numpy()  # 使用全图嵌入

# 构建所有可能的节点对（无向，u < v）
num_nodes = data.num_nodes
all_pairs = []
for u in range(num_nodes):
    for v in range(u+1, num_nodes):
        all_pairs.append([u, v])
all_pairs = np.array(all_pairs)

# 获取现有边集合（无向）
existing_edges = set()
edges_np = data.edge_index.cpu().numpy().T
for u, v in edges_np:
    if u < v:
        existing_edges.add((u, v))
    else:
        existing_edges.add((v, u))

# 筛选未连接的边
candidates = []
for u, v in all_pairs:
    if (u, v) not in existing_edges:
        candidates.append([u, v])
candidates = np.array(candidates)
print(f"Total candidate pairs: {len(candidates)}")

# 计算得分（内积）
scores = []
for u, v in candidates:
    score = np.dot(z[u], z[v])
    scores.append(score)
scores = np.array(scores)

# 取 top-k
top_k_idx = np.argsort(scores)[::-1][:TOP_K_PREDICT]
top_pairs = candidates[top_k_idx]
top_scores = scores[top_k_idx]

# 输出结果
results = []
for i in range(len(top_pairs)):
    u, v = top_pairs[i]
    gene_u = idx2gene[u]
    gene_v = idx2gene[v]
    results.append([gene_u, gene_v, top_scores[i]])

pred_df = pd.DataFrame(results, columns=['gene1', 'gene2', 'score'])

# 第一次合并：为 gene1 添加基因符号
pred_df = pred_df.merge(
    nodes_df[['string_id', 'gene']],
    left_on='gene1',
    right_on='string_id',
    how='left'
).rename(columns={'gene': 'gene1_symbol'}).drop('string_id', axis=1)

# 第二次合并：为 gene2 添加基因符号
pred_df = pred_df.merge(
    nodes_df[['string_id', 'gene']],
    left_on='gene2',
    right_on='string_id',
    how='left'
).rename(columns={'gene': 'gene2_symbol'}).drop('string_id', axis=1)

# 选择需要的列（保留原始 STRING ID 和基因符号）
pred_df = pred_df[['gene1', 'gene2', 'score', 'gene1_symbol', 'gene2_symbol']]
#保存
pred_df.to_csv('results/predicted_interactions.csv', index=False)
##
print(f"Top {TOP_K_PREDICT} predicted interactions saved to 'predicted_interactions.csv'.")
print("\nFirst 10 predictions:")
print(pred_df.head(10))
