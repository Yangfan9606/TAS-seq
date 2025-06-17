#!/usr/bin/env python
# encoding: utf-8
"""
@author:    Yangfan Zhou
@email:     yangfanzhou9606@gmail.com
@date:      2025-06-17
@version:   1.0
@license:   MIT
@brief:     使用21bp序列数据，用于预测11位的A是否为TadA-8e修饰位点。

"""
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import sys
import os
import json
import argparse
from typing import List, Tuple, Dict, Union

class RNASequenceProcessor:
    """RNA序列数据预处理类"""

    def __init__(self, window_size: int = 21):
        self.window_size = window_size
        self.base_to_idx = {'A': 0, 'U': 1, 'G': 2, 'C': 3, 'N': 4}  # N表示未知碱基
        self.padding_idx = 4  # 使用N作为padding

    def encode_sequence(self, sequence: str) -> torch.Tensor:
        """将RNA序列编码为one-hot向量"""
        # 转换为大写并将T替换为U（如果输入的是DNA序列）
        sequence = sequence.upper().replace('T', 'U')

        # 转换为索引
        indices = [self.base_to_idx.get(base, self.base_to_idx['N']) for base in sequence]

        # 转换为one-hot编码
        one_hot = torch.zeros((len(sequence), 5))  # 5种碱基类型
        for i, idx in enumerate(indices):
            one_hot[i, idx] = 1

        return one_hot

    def validate_sequence(self, sequence: str) -> Tuple[bool, str]:
        """验证序列格式"""
        sequence = sequence.upper().strip()

        # 检查长度
        if len(sequence) != 21:
            return False, f"序列长度应为21bp，实际为{len(sequence)}bp"

        # 检查碱基类型
        valid_bases = set('AUCGT')  # 允许DNA和RNA碱基
        invalid_bases = set(sequence) - valid_bases
        if invalid_bases:
            return False, f"包含无效碱基: {', '.join(invalid_bases)}"

        # 转换为RNA格式并检查中心位置
        rna_sequence = sequence.replace('T', 'U')
        if rna_sequence[10] != 'A':
            return False, f"序列中心位置应为'A'，实际为'{rna_sequence[10]}'"

        return True, ""


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


class RNAPredictor:
    """RNA修饰预测器"""

    def __init__(self, model_path: str):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model = None
        self.processor = None
        self.load_model(model_path)

    def load_model(self, model_path: str):
        """加载训练好的模型"""
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"模型文件不存在: {model_path}")

        print(f"加载模型: {model_path}")
        checkpoint = torch.load(model_path, map_location=self.device)

        # 获取模型配置
        model_config = checkpoint['model_config']
        processor_config = checkpoint['processor_config']

        # 创建处理器
        self.processor = RNASequenceProcessor(window_size=processor_config['window_size'])
        self.processor.base_to_idx = processor_config['base_to_idx']

        # 创建模型
        self.model = RNAClassifier(
            input_size=model_config['input_size'],
            hidden_size=model_config['hidden_size'],
            num_filters=model_config['num_filters'],
            filter_sizes=model_config['filter_sizes'],
            dropout_rate=model_config['dropout_rate']
        )

        # 加载模型权重
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.to(self.device)
        self.model.eval()

        print(f"模型加载成功! 使用设备: {self.device}")

    def predict_single(self, sequence: str) -> Dict:
        """预测单个序列"""
        # 验证序列
        is_valid, error_msg = self.processor.validate_sequence(sequence)
        if not is_valid:
            return {
                'sequence': sequence,
                'error': error_msg,
                'is_modified': None,
                'probability': None,
                'confidence': None
            }

        # 预处理序列
        sequence = sequence.upper().replace('T', 'U')  # 转换为RNA
        encoded_seq = self.processor.encode_sequence(sequence)

        # 添加batch维度
        input_tensor = encoded_seq.unsqueeze(0).to(self.device)

        # 预测
        with torch.no_grad():
            output = self.model(input_tensor)
            probability = torch.sigmoid(output).item()
            prediction = probability > 0.5
            confidence = max(probability, 1 - probability)  # 置信度

        return {
            'sequence': sequence,
            'error': None,
            'is_modified': bool(prediction),
            'probability': float(probability),
            'confidence': float(confidence)
        }

    def predict_batch(self, sequences: List[str]) -> List[Dict]:
        """预测多个序列"""
        results = []

        print(f"开始预测 {len(sequences)} 个序列...")

        for i, seq in enumerate(sequences):
            if i % 100 == 0 and i > 0:
                print(f"已处理 {i}/{len(sequences)} 个序列")

            result = self.predict_single(seq)
            results.append(result)

        print("预测完成!")
        return results

    def predict_from_file(self, file_path: str) -> List[Dict]:
        """从文件读取序列并预测"""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"输入文件不存在: {file_path}")

        sequences = []
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith('#'):  # 跳过空行和注释行
                    # 如果行包含多个字段，只取第一个字段作为序列
                    sequence = line.split()[0]
                    sequences.append(sequence)

        print(f"从文件 {file_path} 读取到 {len(sequences)} 个序列")
        return self.predict_batch(sequences)


def save_results(results: List[Dict], output_file: str):
    """保存预测结果到文件"""
    # 创建DataFrame
    df_data = []
    for result in results:
        df_data.append({
            'sequence': result['sequence'],
            'is_modified': result['is_modified'],
            'probability': result['probability'],
            'confidence': result['confidence'],
            'error': result['error']
        })

    df = pd.DataFrame(df_data)

    # 添加文件头说明
    header = [
        "# RNA修饰预测结果",
        f"# 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "# 列说明:",
        "# sequence: 输入的21bp RNA序列",
        "# is_modified: 预测结果 (True=修饰, False=未修饰)",
        "# probability: 修饰概率 (0-1)",
        "# confidence: 预测置信度 (0.5-1.0)",
        "# error: 错误信息 (如果有)",
        "#"
    ]

    # 保存文件
    with open(output_file, 'w') as f:
        f.write('\n'.join(header) + '\n')
        df.to_csv(f, index=False)

    print(f"结果已保存到 {output_file}")


def print_result(result: Dict):
    """打印单个预测结果"""
    print("\n" + "="*50)
    print(f"序列: {result['sequence']}")

    if result['error']:
        print(f"错误: {result['error']}")
    else:
        modification_status = "修饰" if result['is_modified'] else "未修饰"
        print(f"预测结果: {modification_status}")
        print(f"修饰概率: {result['probability']:.4f}")
        print(f"置信度: {result['confidence']:.4f}")


def print_summary(results: List[Dict]):
    """打印预测结果摘要"""
    total = len(results)
    valid_results = [r for r in results if r['error'] is None]
    error_count = total - len(valid_results)

    if valid_results:
        modified_count = sum(1 for r in valid_results if r['is_modified'])
        unmodified_count = len(valid_results) - modified_count
        avg_probability = np.mean([r['probability'] for r in valid_results])
        avg_confidence = np.mean([r['confidence'] for r in valid_results])
    else:
        modified_count = unmodified_count = 0
        avg_probability = avg_confidence = 0

    print("\n" + "="*50)
    print("预测结果摘要:")
    print(f"总序列数: {total}")
    print(f"有效预测: {len(valid_results)}")
    print(f"错误序列: {error_count}")
    print(f"预测为修饰: {modified_count}")
    print(f"预测为未修饰: {unmodified_count}")
    print(f"平均修饰概率: {avg_probability:.4f}")
    print(f"平均置信度: {avg_confidence:.4f}")


def main():
    parser = argparse.ArgumentParser(description='RNA修饰预测工具')
    parser.add_argument('--model', '-m', type=str, default='rna_classifier_model.pth',
                        help='模型文件路径 (默认: rna_classifier_model.pth)')
    parser.add_argument('--input', '-i', type=str,
                        help='输入文件路径 (每行一个21bp序列)')
    parser.add_argument('--sequence', '-s', type=str,
                        help='单个21bp序列')
    parser.add_argument('--output', '-o', type=str,
                        help='输出文件路径 (可选)')

    args = parser.parse_args()

    # 检查输入参数
    if not args.input and not args.sequence:
        print("错误: 必须提供 --input 或 --sequence 参数")
        parser.print_help()
        sys.exit(1)

    if args.input and args.sequence:
        print("错误: --input 和 --sequence 参数不能同时使用")
        sys.exit(1)

    try:
        # 创建预测器
        predictor = RNAPredictor(args.model)

        # 执行预测
        if args.sequence:
            # 单序列预测
            print(f"预测序列: {args.sequence}")
            result = predictor.predict_single(args.sequence)
            print_result(result)

            # 保存结果（如果指定了输出文件）
            if args.output:
                save_results([result], args.output)

        else:
            # 文件预测
            results = predictor.predict_from_file(args.input)

            # 打印摘要
            print_summary(results)

            # 保存结果
            output_file = args.output or f"prediction_results_{time.strftime('%Y%m%d_%H%M%S')}.csv"
            save_results(results, output_file)

            # 如果结果不多，也打印详细结果
            if len(results) <= 10:
                print("\n详细结果:")
                for result in results:
                    print_result(result)

    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # 如果没有命令行参数，提供交互式使用方式
    if len(sys.argv) == 1:
        print("RNA修饰预测工具")
        print("="*50)
        print("使用方法:")
        print("1. 预测单个序列:")
        print("   python rna_prediction.py -s ATCAATATCTATAACTGCATT")
        print()
        print("2. 预测文件中的序列:")
        print("   python rna_prediction.py -i sequences.txt")
        print()
        print("3. 指定输出文件:")
        print("   python rna_prediction.py -i sequences.txt -o results.csv")
        print()
        print("4. 指定模型文件:")
        print("   python rna_prediction.py -m my_model.pth -s ATCAATATCTATAACTGCATT")
        print()
        print("输入文件格式:")
        print("- 每行一个21bp序列")
        print("- 支持DNA和RNA格式 (T会自动转换为U)")
        print("- 中心位置必须为A")
        print("- 支持注释行 (以#开头)")
        print()

        # 交互式输入
        try:
            sequence = input("请输入21bp序列进行预测 (或按Enter退出): ").strip()
            if sequence:
                predictor = RNAPredictor('rna_classifier_model.pth')
                result = predictor.predict_single(sequence)
                print_result(result)
        except KeyboardInterrupt:
            print("\n退出程序")
        except Exception as e:
            print(f"错误: {e}")
    else:
        main()
