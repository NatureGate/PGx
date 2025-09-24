import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import glob

# 设置文件夹路径
folder_path = './pgx_results'
print(os.getcwd())
# 收集所有 .xls 和 .report.tsv 文件
xls_files = glob.glob(os.path.join(folder_path, '*.xls'))
# report_files = glob.glob(os.path.join(folder_path, '*.report.tsv'))
report_list = []
# 遍历每个 .xls 文件
count = 0
for xls_path in xls_files:
    count += 1
    print(f"Processing {count}")
    # 提取文件名（不带扩展名）
    # base_name_full = os.path.splitext(os.path.basename(xls_path))[0]
    # main_base = base_name_full.split('.')[0]  # 去除 .drug 部分
    # report_name = f"{main_base}.report.tsv"
    # report_path = os.path.join(folder_path, report_name)

    
    # 检查对应的 .report.tsv 文件是否存在
    if os.path.exists(xls_path):
        #print(f"Processing pair:  {report_path}")
        
        # 读取 .xls 文件
        df_drugxls = pd.read_table(xls_path, sep='\t')
        gene_mask = ~df_drugxls['gene'].isin(['CYP2D6', 'HLA-A', 'HLA-B'])
        diplotype_mask = ~(df_drugxls['diplotype'].isin(['*1/*1','.']))
        filtered_df = df_drugxls[gene_mask & diplotype_mask]
        report_list.append(filtered_df)

# print(f"Total reports processed: {len(report_list)}")
# 合并所有报告的结果
combined_report = pd.concat(report_list, ignore_index=True)
# 计算出每个基因的Source Diplotype 不包含 ' AND ' 或者 ' OR ' 的概率
gene_column = 'gene'
source_diplotype_column = 'diplotype'
combined_report = combined_report.groupby(['sample','gene']).first()
# === 计算 total_counts ===
total_counts = combined_report.groupby(gene_column).size().reset_index(name='total_counts')
sorted_data = total_counts.sort_values('total_counts', ascending=False)
sorted_data.to_excel(os.path.join(folder_path, 'old_sorted_total_counts.xlsx'), index=False)
# === 计算概率（不含 ' AND ' 或 ' OR ' 的记录比例）===
#mask = ~combined_report[source_diplotype_column].str.contains(' AND | OR ', na=False, case=False)
# valid_counts = combined_report.groupby(gene_column)[source_diplotype_column].apply(
#     lambda x: ((~x.str.contains(' AND ', case=False)) & (~x.str.contains(' OR ', case=False))).sum()
# ).reset_index(name='valid_counts')

# 合并数据框并计算概率
# merged_df = pd.merge(total_counts, valid_counts, on=gene_column)
# merged_df['probability'] = (merged_df['valid_counts'] / merged_df['total_counts'] * 100).round(2)

# 按 total_counts 排序
#merged_sorted = merged_df.sort_values('total_counts', ascending=False)

# === 绘制分组柱状图（双纵轴）===
plt.figure(figsize=(12, 6))
ax1 = plt.gca()  # 左纵轴（数量）


# 生成x轴的位置
x = np.arange(len(sorted_data))
width = 0.4  # 柱子宽度

# 绘制左轴（total_counts）
bars1 = ax1.bar(
    x - width/2,        # 左侧柱子偏移
    sorted_data['total_counts'],
    width,
    color='skyblue',
    label='Total Counts'
)



# 设置坐标轴标签和标题
ax1.set_xlabel('Gene', fontsize=12)
ax1.set_ylabel('Counts', color='darkblue', fontsize=12)
ax1.tick_params(axis='y', labelcolor='darkblue')



# 设置x轴标签（基因名称）
ax1.set_xticks(x)
ax1.set_xticklabels(sorted_data[gene_column], rotation=45, ha='right')

# 添加图例
lines1, labels1 = ax1.get_legend_handles_labels()

combined_lines = lines1 
combined_labels = labels1 

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


# 保存并显示
plt.tight_layout()
output_path = os.path.join(folder_path, 'grouped_bar_plot_old.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"图表已保存至：{output_path}")
plt.show()