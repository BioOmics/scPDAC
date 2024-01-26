#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Haoyu Chao (haoyuchao@zju.edu.cn)
'''
===================
一、脚本使用方法
===================
1. 先使用pancancer账户切换到gpuserver节点:
ssh gpuserver
2. 再使用绝对路径调用conda激活环境:
source /public/workspace202011/pancancer/chaohy/software/miniconda3/bin/activate gpu
3. 调用脚本处理demo数据
python /public/workspace202011/pancancer/chaohy/scripts/auto_nmf.py \
    -i /public/workspace202011/pancancer/chaohy/NMF_count/demo/demo.count.txt \
    -W /public/workspace202011/pancancer/chaohy/NMF_count/demo/demo.w.txt \
    -H /public/workspace202011/pancancer/chaohy/NMF_count/demo/demo.h.txt \
    --rank_min 1 --rank_max 30 -n 3000
4. 查看脚本帮助文档
python /public/workspace202011/pancancer/chaohy/scripts/auto_nmf.py -h
===================
二、其他有帮助的命令
===================
1. 查看GPU使用情况, 每秒刷新一次. ctrl+c退出
watch -n 1 nvidia-smi
2. 想回到node5节点, 需要先切换到server, 再切到node5。:)
ssh server
ssh node5
'''
print('Loading modules...')
import argparse
import torch
import pandas as pd
import matplotlib.pyplot as plt

# 定义非负矩阵分解函数
def nmf(data, rank, num_iterations, tol=1e-5, seed=520):

    torch.manual_seed(seed) # set random seed

    # Convert data to tensor
    data_tensor = torch.tensor(data.values, dtype=torch.float32).cuda() # fp32

    # Random initialization of factor matrices
    m, n = data_tensor.shape
    W = torch.abs(torch.randn((m, rank), device='cuda'))
    H = torch.abs(torch.randn((rank, n), device='cuda'))
    # print('W shape: {}, H shape: {}'.format(W.shape, H.shape))

    # Iterative update of factor matrices
    prev_recon_error = float('inf') # Set previous reconstruction error to infinity
    for i in range(num_iterations):
        # Update H
        numerator = torch.matmul(W.t(), data_tensor)
        denominator = torch.matmul(torch.matmul(W.t(), W), H)
        denominator = torch.where(denominator < 1.3e-7, torch.tensor(1.3e-7, device='cuda'), denominator) # Protect the denominator
        H *= numerator / denominator

        # Update W
        numerator = torch.matmul(data_tensor, H.t())
        denominator = torch.matmul(torch.matmul(W, H), H.t())
        denominator = torch.where(denominator < 1.3e-7, torch.tensor(1.3e-7, device='cuda'), denominator) # Protect the denominator
        W *= numerator / denominator

        # Calculate reconstruction error and check for convergence
        if i % 10 == 0:
            recon_error = torch.norm(data_tensor - torch.matmul(W, H), p='fro')
            rel_change = torch.abs(recon_error - prev_recon_error) / prev_recon_error
            prev_recon_error = recon_error

            # print('iterations: {}, Reconstruction error: {:.2f}, Relative change: {:.2e}'.format(i, recon_error, rel_change))
            if rel_change <= tol:
                print('Converged after {} iterations, Recon_error: {:.2f}, Tol: {:.2e}'.format(i, recon_error, rel_change))
                break

        if i == num_iterations - 1:
            print('Reached maximum number of iterations, Recon_error: {:.2f}, Tol: {:.2e}'.format(recon_error, rel_change))
    
    # Convert results to pandas DataFrame and add row and column names
    W_df = pd.DataFrame(W.cpu().numpy(), index=data.index)
    H_df = pd.DataFrame(H.cpu().numpy(), columns=data.columns)

    return W_df, H_df, recon_error


def find_optimal_rank(data, rank_min, rank_max, num_iterations):
    best_rank = None
    similarity_values = []
    recon_error_values = []
    score_values = []
    
    rank_values = range(rank_min, rank_max + 1)

    for rank in range(rank_min, rank_max + 1):
        print('Evaluating rank: {}'.format(rank))

        W, H, recon_error = nmf(data, rank, num_iterations)
        recon_error_values.append(recon_error.item())

        # Convert data to tensor
        data_tensor = torch.tensor(data.values, dtype=torch.float32).cuda() # fp32
        W = torch.tensor(W.values, dtype=torch.float32).cuda() # fp32
        H = torch.tensor(H.values, dtype=torch.float32).cuda() # fp32

        # Calculate similarity between data and reconstruction
        similarity = torch.cosine_similarity(data_tensor, torch.matmul(W, H), dim=1) # 按细胞所在的行计算余弦相似度
        similarity_avg = torch.mean(similarity)  # 使用平均值作为衡量相似性的指标
        similarity_values.append(similarity_avg.item())

        print('Similarity: {}'.format(similarity_avg))

    # 将重构误差归一化到0-1之间
    scaled_recon_error_values = scale_values(recon_error_values) 

    # Calculate the combined score as weighted average
    alpha = 0.3 # 重构误差的权重
    beta = 0.7 # 相似性的权重
    score_values = [alpha * (1 - x) + beta * y for x, y in zip(scaled_recon_error_values, similarity_values)]
    for i in range(len(score_values)):
        if i > 0:
            if abs(score_values[i] - score_values[i-1]) < 0.001:
                best_rank = i + rank_min
                best_similarity = similarity_values[i]
                break

    if best_rank is None:
        best_rank = score_values.index(max(score_values))
        best_similarity = similarity_values[best_rank]
    
    # 保存结果图像为PDF
    plot_results(rank_values, similarity_values, scaled_recon_error_values, score_values, 'result_plot.pdf', best_rank)

    return best_rank, best_similarity

def plot_results(rank_values, similarity_values, scaled_recon_error_values, score_values, output_file, best_rank):
    fig, ax1 = plt.subplots(figsize=(len(rank_values) * 0.4, 4))

    color = 'tab:green'
    ax1.set_xlabel('Rank')
    ax1.set_ylabel('Similarity', color=color)
    ax1.plot(rank_values, similarity_values, color=color, marker='o')
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('Error', color=color)
    ax2.plot(rank_values, scaled_recon_error_values, color=color, marker='o')
    ax2.tick_params(axis='y', labelcolor=color)

    color='tab:red'
    ax3 = ax1.twinx()  # instantiate a third axes that shares the same x-axis
    ax3.spines['right'].set_position(('outward', 60))  # move the third axis to the right
    ax3.set_ylabel('Score', color=color)
    ax3.plot(rank_values, score_values, color=color, marker='o')
    ax3.tick_params(axis='y', labelcolor=color)  # Set label color to red

    plt.xticks(range(min(rank_values), max(rank_values) + 1), range(min(rank_values), max(rank_values) + 1))

    plt.axvline(x=best_rank, color='purple', linestyle='--')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(output_file, format='pdf')
    plt.close()

def scale_values(values):
    min_val = min(values)
    max_val = max(values)
    scaled_values = [(val - min_val) / (max_val - min_val) for val in values]
    return scaled_values

if __name__ == '__main__':
    # 定义参数
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, help='输入文件路径')
    parser.add_argument('-W', '--output_w', type=str, required=True, help='输出W的文件路径')
    parser.add_argument('-H', '--output_h', type=str, required=True, help='输出H的文件路径')
    parser.add_argument('--rank_min', type=int, default=2, help='矩阵分解的最小秩')
    parser.add_argument('--rank_max', type=int, default=30, help='矩阵分解的最大秩')
    parser.add_argument('-n', '--num_iterations', type=int, default=3000, help='迭代次数')
    args = parser.parse_args()

    # 读取数据
    print('Reading data...')
    data = pd.read_csv(args.input, sep='\t', index_col=0)

    # 寻找最优的rank值
    best_rank, best_similarity = find_optimal_rank(data, args.rank_min, args.rank_max, args.num_iterations)

    # 进行非负矩阵分解
    print('=====================================================')
    print('Performing NMF with Best rank: {}, Best similarity: {}'.format(best_rank, best_similarity))
    W, H, _ = nmf(data, best_rank, args.num_iterations)

    # 保存结果
    print('Saving results...')
    W.to_csv(args.output_w, sep='\t')
    H.to_csv(args.output_h, sep='\t')
