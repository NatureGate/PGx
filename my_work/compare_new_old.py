import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import glob

# 设置文件夹路径
folder_path = './pgx_results'
# print(os.getcwd())
# # 收集所有 .xls 和 .report.tsv 文件
# xls_files = glob.glob(os.path.join(folder_path, '*.xls'))
# report_files = glob.glob(os.path.join(folder_path, '*.report.tsv'))
# report_list = []
# # 遍历每个 .report.tsv 文件
# count = 0
# for report_path in report_files:
#     count += 1
#     print(f"Processing {count}")
#     # 提取文件名（不带扩展名）
#     base_name_full = os.path.splitext(os.path.basename(report_path))[0]
#     main_base = base_name_full.split('.')[0]  # 去除 .drug 部分
#     report_name = f"{main_base}.report.tsv"
#     report_path = os.path.join(folder_path, report_name)

    
#     # 检查对应的 .report.tsv 文件是否存在
#     if os.path.exists(report_path):
#         #print(f"Processing pair:  {report_path}")
        
#         # 读取 .xls 文件
#         # df_xls = pd.read_table(xls_path, sep='\t')
#         df_report = pd.read_table(report_path, sep='\t',skiprows=1)
#         gene_mask = ~df_report['Gene'].isin(['CYP2D6', 'HLA-A', 'HLA-B'])
#         diplotype_mask = ~(df_report['Source Diplotype'].isin(['Unknown/Unknown']))
#         filtered_df = df_report[gene_mask & diplotype_mask]
#         report_list.append(filtered_df)

# # print(f"Total reports processed: {len(report_list)}")
# # 合并所有报告的结果
# combined_report = pd.concat(report_list, ignore_index=True)
# # 计算出每个基因的Source Diplotype 不包含 ' AND ' 或者 ' OR ' 的概率
gene_column = 'gene'
# source_diplotype_column = 'Source Diplotype'
# # === 计算 total_counts ===
# total_counts = combined_report.groupby(gene_column).size().reset_index(name='total_counts')

# # === 计算概率（不含 ' AND ' 或 ' OR ' 的记录比例）===
# mask = ~combined_report[source_diplotype_column].str.contains(' AND | OR ', na=False, case=False)
# valid_counts = combined_report.groupby(gene_column)[source_diplotype_column].apply(
#     lambda x: ((~x.str.contains(' AND ', case=False)) & (~x.str.contains(' OR ', case=False))).sum()
# ).reset_index(name='valid_counts')

# # 合并数据框并计算概率
# merged_df = pd.merge(total_counts, valid_counts, on=gene_column)
# merged_df['probability'] = (merged_df['valid_counts'] / merged_df['total_counts'] * 100).round(2)

# 按 total_counts 排序
# merged_sorted = merged_df.sort_values('total_counts', ascending=False)
# merged_sorted.to_excel(os.path.join(folder_path, 'merged_sorted.xlsx'), index=False)

old_sorted = pd.read_excel(os.path.join(folder_path, 'compare_old_new.xlsx'))
old_sorted = old_sorted.sort_values(by=['old','new'], ascending=False)
# === 绘制分组柱状图（双纵轴）===
plt.figure(figsize=(12, 6))
ax1 = plt.gca()  # 左纵轴（数量）
ax2 = ax1.twinx()  # 右纵轴（概率）

# 生成x轴的位置
x = np.arange(len(old_sorted))
width = 0.4  # 柱子宽度

# 绘制左轴（total_counts）
bars1 = ax1.bar(
    x - width/2,        # 左侧柱子偏移
    old_sorted['old'],
    width,
    color='skyblue',
    label='Total Counts'
)

# 绘制右轴（probability）
bars2 = ax2.bar(
    x + width/2,        # 右侧柱子偏移
    old_sorted['new'],
    width,
    color='coral',
    label=''
)

# 设置坐标轴标签和标题
ax1.set_xlabel('gene', fontsize=12)
ax1.set_ylabel('Counts', color='darkblue', fontsize=12)
ax1.tick_params(axis='y', labelcolor='darkblue')

ax2.set_ylabel('', color='darkred', fontsize=12)
ax2.tick_params(axis='y', labelcolor='darkred')
#ax2.set_ylim(0, 110)  # 设置概率轴范围（0-110%）

# 设置x轴标签（基因名称）
ax1.set_xticks(x)
ax1.set_xticklabels(old_sorted[gene_column], rotation=45, ha='right')

# 添加图例
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
combined_lines = lines1 + lines2
combined_labels = labels1 + labels2

# 创建统一的图例
plt.legend(
    combined_lines,
    combined_labels,
    loc='upper right',
    bbox_to_anchor=(0.99, 1.0),  # 右外侧对齐
    frameon=True,
    fontsize=9,                # 统一字体大小
    borderpad=0.4,             # 统一边距
    labelspacing=0.5           # 统一标签间距
)

# 显示数值标签（可选）
def add_labels(bars, ax, is_percent=False):
    for bar in bars:
        height = bar.get_height()
        if is_percent:
            ax.annotate(f'{height:.1f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),  # 3点垂直偏移
                        textcoords="offset points",
                        ha='center', va='bottom',
                        fontsize=6)
        else:
            ax.annotate(f'{int(height)}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom',
                        fontsize=6)

add_labels(bars1, ax1)
add_labels(bars2, ax2)

# 保存并显示
plt.tight_layout()
output_path = os.path.join(folder_path, 'old_new_bar_plot.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"图表已保存至：{output_path}")
plt.show()