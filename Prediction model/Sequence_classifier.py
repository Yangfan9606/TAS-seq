#!/usr/bin/env python
# encoding: utf-8
"""
@author:    Yangfan Zhou
@email:     yangfanzhou9606@gmail.com
@date:      2025-06-17
@version:   1.0
@license:   MIT
@brief:     使用TAS-seq修饰A前后10bp（共21bp）序列数据，用于预测模型训练。

"""
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, roc_auc_score, roc_curve
import os
import time
import json
import sys
from typing import List, Tuple, Dict  # 添加缺失的导入

start_time = time.time()

class RNASequenceProcessor:
    """RNA序列数据预处理类"""

    def __init__(self, window_size: int = 21):
        self.window_size = window_size
        self.base_to_idx = {'A': 0, 'U': 1, 'G': 2, 'C': 3, 'N': 4}  # N表示未知碱基
        self.padding_idx = 4  # 使用N作为padding

    def encode_sequence(self, sequence: str) -> torch.Tensor:
        """将RNA序列编码为one-hot向量"""
        # 转换为大写
        sequence = sequence.upper()

        # 转换为索引
        indices = [self.base_to_idx.get(base, self.base_to_idx['N']) for base in sequence]

        # 转换为one-hot编码
        one_hot = torch.zeros((len(sequence), 5))  # 5种碱基类型
        for i, idx in enumerate(indices):
            one_hot[i, idx] = 1

        return one_hot

    def load_data_from_file(self, file_path: str) -> Tuple[List[str], List[int]]:
        """从文件加载序列和标签数据（新格式）"""
        windows = []
        labels = []

        # 读取文件
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    sequence = parts[0].upper()
                    label = int(parts[1])  # 标签 (0或1)

                    # 验证序列长度是否为21bp
                    if len(sequence) != 21:
                        print(f"警告: 序列长度应为21bp, 实际为{len(sequence)}bp: '{sequence}'")

                    # 验证中心位置是否为A
                    if sequence[10].upper() != 'A':
                        print(f"警告: 序列中心位置应为'A', 实际为'{sequence[10]}': '{sequence}'")

                    windows.append(sequence)
                    labels.append(label)

        return windows, labels

class RNADataset(Dataset):
    """RNA序列数据集类"""

    def __init__(self, sequences: List[str], labels: List[int], processor: RNASequenceProcessor):
        self.sequences = sequences
        self.labels = labels
        self.processor = processor

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        label = self.labels[idx]

        # 编码序列
        encoded_seq = self.processor.encode_sequence(sequence)

        return encoded_seq, torch.tensor(label, dtype=torch.float32)

class RNAClassifier(nn.Module):
    """RNA序列分类模型 - 结合CNN和LSTM"""

    def __init__(self, input_size: int = 5, hidden_size: int = 128, num_filters: int = 64,
                 filter_sizes: List[int] = [3, 5, 7], dropout_rate: float = 0.3):
        super(RNAClassifier, self).__init__()

        self.input_size = input_size
        self.hidden_size = hidden_size
        self.num_filters = num_filters

        # 多尺度CNN层
        self.convs = nn.ModuleList([
            nn.Sequential(
                nn.Conv1d(input_size, num_filters, kernel_size=k, padding=k//2),
                nn.ReLU(),
                nn.BatchNorm1d(num_filters)
            ) for k in filter_sizes
        ])

        # LSTM层
        self.lstm = nn.LSTM(
            input_size=num_filters * len(filter_sizes),
            hidden_size=hidden_size,
            num_layers=2,
            dropout=dropout_rate,
            bidirectional=True,
            batch_first=True
        )

        # 注意力机制
        self.attention = nn.MultiheadAttention(
            embed_dim=hidden_size * 2,
            num_heads=8,
            dropout=dropout_rate,
            batch_first=True
        )

        # 分类层
        self.classifier = nn.Sequential(
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_size * 2, hidden_size),
            nn.ReLU(),
            nn.BatchNorm1d(hidden_size),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_size, 32),
            nn.ReLU(),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        # x shape: (batch_size, seq_len, input_size)
        batch_size, seq_len, _ = x.shape

        # 转换为CNN输入格式: (batch_size, input_size, seq_len)
        x_conv = x.transpose(1, 2)

        # 多尺度CNN特征提取
        conv_outputs = []
        for conv in self.convs:
            conv_out = conv(x_conv)  # (batch_size, num_filters, seq_len)
            conv_outputs.append(conv_out)

        # 拼接不同尺度的特征
        conv_features = torch.cat(conv_outputs, dim=1)  # (batch_size, num_filters*len(filter_sizes), seq_len)

        # 转换为LSTM输入格式: (batch_size, seq_len, features)
        conv_features = conv_features.transpose(1, 2)

        # LSTM特征提取
        lstm_out, (hidden, cell) = self.lstm(conv_features)

        # 自注意力机制
        attn_out, attn_weights = self.attention(lstm_out, lstm_out, lstm_out)

        # 全局平均池化
        pooled = torch.mean(attn_out, dim=1)  # (batch_size, hidden_size*2)

        # 分类
        output = self.classifier(pooled)

        return output.squeeze()

class ModelTrainer:
    """模型训练器"""

    def __init__(self, model: nn.Module, device: torch.device):
        self.model = model.to(device)
        self.device = device
        self.train_losses = []
        self.val_losses = []
        self.val_accuracies = []
        self.val_fpr = []  # 存储每个epoch的假阳性率
        self.val_tpr = []  # 存储每个epoch的真阳性率

    def train_epoch(self, train_loader: DataLoader, optimizer: optim.Optimizer,
                   criterion: nn.Module) -> float:
        """训练一个epoch"""
        self.model.train()
        total_loss = 0

        for sequences, labels in train_loader:
            sequences = sequences.to(self.device)
            labels = labels.to(self.device)

            optimizer.zero_grad()
            outputs = self.model(sequences)
            loss = criterion(outputs, labels)
            loss.backward()

            # 梯度裁剪
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)

            optimizer.step()
            total_loss += loss.item()

        return total_loss / len(train_loader)

    def validate(self, val_loader: DataLoader, criterion: nn.Module) -> Tuple[float, float, List, List]:
        """验证模型并计算ROC曲线数据"""
        self.model.eval()
        total_loss = 0
        all_predictions = []
        all_labels = []
        all_probabilities = []

        with torch.no_grad():
            for sequences, labels in val_loader:
                sequences = sequences.to(self.device)
                labels = labels.to(self.device)

                outputs = self.model(sequences)
                loss = criterion(outputs, labels)
                total_loss += loss.item()

                # 获取预测结果
                probabilities = torch.sigmoid(outputs)
                predictions = probabilities > 0.5

                all_predictions.extend(predictions.cpu().numpy())
                all_probabilities.extend(probabilities.cpu().numpy())
                all_labels.extend(labels.cpu().numpy())

        avg_loss = total_loss / len(val_loader)
        accuracy = accuracy_score(all_labels, all_predictions)

        # 计算ROC曲线
        fpr, tpr, _ = roc_curve(all_labels, all_probabilities)

        return avg_loss, accuracy, fpr, tpr

    def train(self, train_loader: DataLoader, val_loader: DataLoader,
              num_epochs: int = 100, learning_rate: float = 0.001,
              weight_decay: float = 1e-4, patience: int = 10) -> Dict:
        """完整训练流程"""

        # 处理类别不平衡
        criterion = nn.BCEWithLogitsLoss()
        optimizer = optim.Adam(self.model.parameters(), lr=learning_rate, weight_decay=weight_decay)
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', patience=5, factor=0.5)

        best_val_loss = float('inf')
        patience_counter = 0
        best_model_state = None

        print("开始训练...")
        for epoch in range(num_epochs):
            # 训练
            train_loss = self.train_epoch(train_loader, optimizer, criterion)

            # 验证
            val_loss, val_acc, fpr, tpr = self.validate(val_loader, criterion)

            # 学习率调度
            scheduler.step(val_loss)

            # 记录指标
            self.train_losses.append(train_loss)
            self.val_losses.append(val_loss)
            self.val_accuracies.append(val_acc)
            self.val_fpr.append(fpr.tolist())  # 存储为列表
            self.val_tpr.append(tpr.tolist())  # 存储为列表

            # 早停机制
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
                best_model_state = self.model.state_dict().copy()
            else:
                patience_counter += 1

            if epoch % 10 == 0:
                print(f'Epoch {epoch:3d}: Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.4f}')

            if patience_counter >= patience:
                print(f'早停于第 {epoch} 轮')
                break

        # 恢复最佳模型
        if best_model_state is not None:
            self.model.load_state_dict(best_model_state)

        return {
            'train_losses': self.train_losses,
            'val_losses': self.val_losses,
            'val_accuracies': self.val_accuracies,
            'val_fpr': self.val_fpr,
            'val_tpr': self.val_tpr,
            'best_val_loss': best_val_loss
        }

class ModelEvaluator:
    """模型评估器"""

    def __init__(self, model: nn.Module, device: torch.device):
        self.model = model.to(device)
        self.device = device

    def evaluate(self, test_loader: DataLoader) -> Dict:
        """在测试集上评估模型并计算ROC曲线数据"""
        self.model.eval()
        all_predictions = []
        all_probabilities = []
        all_labels = []
        all_sequences = []  # 保存原始序列

        with torch.no_grad():
            for sequences, labels in test_loader:
                # 保存原始序列
                batch_sequences = sequences  # 直接保存序列

                sequences = sequences.to(self.device)
                labels = labels.to(self.device)

                outputs = self.model(sequences)
                probabilities = torch.sigmoid(outputs)
                predictions = probabilities > 0.5

                all_predictions.extend(predictions.cpu().numpy())
                all_probabilities.extend(probabilities.cpu().numpy())
                all_labels.extend(labels.cpu().numpy())
                all_sequences.extend(batch_sequences)  # 保存原始序列

        # 计算评估指标
        accuracy = accuracy_score(all_labels, all_predictions)
        precision, recall, f1, _ = precision_recall_fscore_support(all_labels, all_predictions, average='binary')
        auc_roc = roc_auc_score(all_labels, all_probabilities)

        # 计算ROC曲线
        fpr, tpr, thresholds = roc_curve(all_labels, all_probabilities)

        return {
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'auc_roc': auc_roc,
            'predictions': all_predictions,
            'probabilities': all_probabilities,
            'labels': all_labels,
            'sequences': all_sequences,
            'fpr': fpr.tolist(),
            'tpr': tpr.tolist(),
            'thresholds': thresholds.tolist()
        }

def save_training_history(history: Dict, filename: str = "training_history.csv"):
    """保存训练历史到CSV文件（包含每个epoch的ROC数据）"""
    # 创建基础数据框
    df = pd.DataFrame({
        'epoch': range(len(history['train_losses'])),
        'train_loss': history['train_losses'],
        'val_loss': history['val_losses'],
        'val_accuracy': history['val_accuracies']
    })

    # 添加文件头说明
    header = [
        "# 训练历史记录",
        f"# 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "# 列说明:",
        "# epoch: 训练轮次",
        "# train_loss: 训练集损失值",
        "# val_loss: 验证集损失值",
        "# val_accuracy: 验证集准确率"
    ]

    # 保存文件
    with open(filename, 'w') as f:
        f.write('\n'.join(header) + '\n')
        df.to_csv(f, index=False)

    print(f"训练历史已保存到 {filename}")

    # 保存ROC曲线数据（每个epoch的FPR和TPR）
    roc_data = {
        'epochs': list(range(len(history['val_fpr']))),
        'fpr': history['val_fpr'],
        'tpr': history['val_tpr']
    }

    roc_filename = "validation_roc_data.json"
    with open(roc_filename, 'w') as f:
        json.dump(roc_data, f, indent=4)

    print(f"验证集ROC曲线数据已保存到 {roc_filename}")

def save_test_results(results: Dict, filename: str = "test_results.csv"):
    """保存测试结果到CSV文件"""
    # 创建结果DataFrame
    result_df = pd.DataFrame({
        'sequence': results['sequences'],
        'true_label': results['labels'],
        'prediction': results['predictions'],
        'probability': results['probabilities']
    })

    # 添加评估指标
    metrics = {
        'accuracy': results['accuracy'],
        'precision': results['precision'],
        'recall': results['recall'],
        'f1_score': results['f1_score'],
        'auc_roc': results['auc_roc']
    }

    # 添加文件头说明
    header = [
        "# 模型测试结果",
        f"# 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "# 列说明:",
        "# sequence: 输入的21bp RNA序列",
        "# true_label: 真实标签 (0=未修饰, 1=修饰)",
        "# prediction: 模型预测标签 (0=未修饰, 1=修饰)",
        "# probability: 修饰概率 (0-1)",
        "#",
        "# 评估指标汇总:",
        *[f"# {k}: {v:.4f}" for k, v in metrics.items()]
    ]

    # 保存文件
    with open(filename, 'w') as f:
        f.write('\n'.join(header) + '\n')
        result_df.to_csv(f, index=False)

    print(f"测试结果已保存到 {filename}")

    # 保存ROC曲线数据
    roc_data = {
        'fpr': results['fpr'],
        'tpr': results['tpr'],
        'thresholds': results['thresholds'],
        'auc': results['auc_roc']
    }

    roc_filename = "test_roc_data.json"
    with open(roc_filename, 'w') as f:
        json.dump(roc_data, f, indent=4)

    print(f"测试集ROC曲线数据已保存到 {roc_filename}")

    # 单独保存评估指标到JSON文件
    metrics_filename = "evaluation_metrics.json"
    with open(metrics_filename, 'w') as f:
        json.dump(metrics, f, indent=4)

    print(f"评估指标已保存到 {metrics_filename}")

def save_model_config(model, processor, filename: str = "model_config.json"):
    """保存模型配置信息"""
    # 获取卷积核大小
    filter_sizes = []
    for conv in model.convs:
        # 获取卷积层的kernel_size
        for module in conv:
            if isinstance(module, nn.Conv1d):
                filter_sizes.append(module.kernel_size[0])

    config = {
        'model_params': {
            'input_size': model.input_size,
            'hidden_size': model.hidden_size,
            'num_filters': model.num_filters,
            'filter_sizes': filter_sizes,
            'dropout_rate': model.lstm.dropout
        },
        'processor_params': {
            'window_size': processor.window_size,
            'base_to_idx': processor.base_to_idx
        },
        'save_time': time.strftime('%Y-%m-%d %H:%M:%S')
    }

    with open(filename, 'w') as f:
        json.dump(config, f, indent=4)

    print(f"模型配置已保存到 {filename}")

def main(input_file: str):
    """主函数 - 完整的训练和评估流程"""

    # 设置随机种子
    torch.manual_seed(42)
    np.random.seed(42)

    # 设置设备
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'使用设备: {device}')
    print(f'输入文件: {input_file}')

    # 1. 数据预处理
    print("1. 数据预处理...")
    processor = RNASequenceProcessor(window_size=21)

    # 从文件加载数据
    windows, labels = processor.load_data_from_file(input_file)

    print(f"总样本数: {len(windows)}")
    print(f"正样本数: {sum(labels)}")
    print(f"负样本数: {len(labels) - sum(labels)}")

    # 2. 数据划分
    print("2. 数据划分...")
    X_temp, X_test, y_temp, y_test = train_test_split(
        windows, labels, test_size=0.2, stratify=labels, random_state=42
    )
    X_train, X_val, y_train, y_val = train_test_split(
        X_temp, y_temp, test_size=0.25, stratify=y_temp, random_state=42
    )

    print(f"训练集: {len(X_train)}, 验证集: {len(X_val)}, 测试集: {len(X_test)}")

    # 3. 创建数据集和数据加载器
    train_dataset = RNADataset(X_train, y_train, processor)
    val_dataset = RNADataset(X_val, y_val, processor)
    test_dataset = RNADataset(X_test, y_test, processor)

    batch_size = 32
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    # 4. 创建模型
    print("3. 创建模型...")
    model = RNAClassifier(
        input_size=5,
        hidden_size=128,
        num_filters=64,
        filter_sizes=[3, 5, 7],
        dropout_rate=0.3
    )

    print(f"模型参数数量: {sum(p.numel() for p in model.parameters()):,}")

    # 5. 训练模型
    print("4. 训练模型...")
    trainer = ModelTrainer(model, device)
    training_history = trainer.train(
        train_loader=train_loader,
        val_loader=val_loader,
        num_epochs=50,
        learning_rate=0.001,
        patience=10
    )

    # 6. 评估模型
    print("5. 评估模型...")
    evaluator = ModelEvaluator(model, device)
    test_results = evaluator.evaluate(test_loader)

    print("\n=== 测试集结果 ===")
    print(f"准确率: {test_results['accuracy']:.4f}")
    print(f"精确率: {test_results['precision']:.4f}")
    print(f"召回率: {test_results['recall']:.4f}")
    print(f"F1分数: {test_results['f1_score']:.4f}")
    print(f"AUC-ROC: {test_results['auc_roc']:.4f}")

    # 7. 保存结果
    print("6. 保存结果...")
    save_training_history(training_history, "training_history.csv")
    save_test_results(test_results, "test_results.csv")
    save_model_config(model, processor, "model_config.json")

    # 8. 保存模型
    model_save_path = "rna_classifier_model.pth"
    torch.save({
        'model_state_dict': model.state_dict(),
        'processor_config': {
            'window_size': processor.window_size,
            'base_to_idx': processor.base_to_idx
        },
        'model_config': {
            'input_size': 5,
            'hidden_size': 128,
            'num_filters': 64,
            'filter_sizes': [3, 5, 7],
            'dropout_rate': 0.3
        }
    }, model_save_path)

    print(f"模型已保存到 {model_save_path}")
    end_time = time.time()
    run_time = end_time - start_time
    print(f"程序运行时间为: {run_time} 秒")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python script.py <输入文件路径>")
        print("输入文件格式: 每行包含21bp序列和标签，以空格分隔")
        print("示例: ATCAATATCTATAACTGCATT 1")
        sys.exit(1)

    input_file = sys.argv[1]
    main(input_file)
