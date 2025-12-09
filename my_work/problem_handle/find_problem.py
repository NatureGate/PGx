import pandas as pd


# 读取 old_new_genotypes.csv
df = pd.read_csv('problem_handle/old_new_genotypes.csv')

# 找出所有列中 isna 或值为 'Unknown/Unknown' 的行
problem_mask = (df['MT-RNR1_x'] == 'Reference') & (df['MT-RNR1_y'] != '.')
mask = df.isna() | (df == 'Unknown/Unknown')|problem_mask
# 找到MT-RNR1_x为Reference但是MT-RNR1_y不是.的行
# 获取存在问题的样本索引
problem_idx = mask.any(axis=1)

# 提取对应的 sample_id
problem_samples = df.loc[problem_idx, 'sample_id'].dropna().unique()

# 输出结果
print("存在缺失值或 Unknown/Unknown 的样本 ID：")
for sid in problem_samples:
    print(sid)
    
# 读取 PGx-1112.csv
pgx_df = pd.read_csv('problem_handle/PGx-1112.csv')

# 根据 problem_samples 过滤
filtered_df = pgx_df[pgx_df['EntityID'].isin(problem_samples)]

# 输出新的文件
filtered_df.to_csv('problem_handle/PGx-1125_problem_samples.csv', index=False)
print("已生成包含 problem_samples 的新文件：PGx-1125_problem_samples.csv")


