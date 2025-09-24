import pandas as pd
import re

# # 读文件
# df = pd.read_csv('var_transcript_data/combin_all.csv', dtype=str)

# # 1) 只保留 REF 和 ALT 都是单碱基 A/C/T/G 的行
# valid_bases = {'A', 'C', 'T', 'G'}
# single_df = df[df['REF'].isin(valid_bases) & df['ALT'].isin(valid_bases)]
# multi_df = df[~(df['REF'].isin(valid_bases) & df['ALT'].isin(valid_bases))]
# # 2) 定义映射关系
# MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# def _should_drop(row):
#     ref, alt = row['REF'], row['ALT']
#     pat1 = f'{ref}>{alt}'
#     pat2 = f'{MAP[ref]}>{MAP[alt]}'
#     name = str(row['Name'])
#     # 如果 Name 中既没有 pat1 也没有 pat2，则返回 True（表示不符合）
#     return not (re.search(pat1, name, re.I) or re.search(pat2, name, re.I))

# # 3) 筛选出“不符合”条件的记录
# df_good = single_df[~single_df.apply(_should_drop, axis=1)]

# # 查看结果
# print(df_good.shape)
# # 可选：保存
# df_good.to_csv('var_transcript_data/good_records.csv', index=False)

# pd.concat([df_good, multi_df]).to_csv('var_transcript_data/filter_combine_all.csv', index=False)


import pandas as pd

# 读取 CSV
df = pd.read_csv('var_transcript_data/filter_combine_all.csv', dtype=str)

# 清洗 ALT 列：去掉引号和空格
df['ALT'] = df['ALT'].str.replace('"', '').str.replace(' ', '')

def split_alt_by_group(group):
    # 拆分所有 ALT
    alt_list = group['ALT'].str.split(',').explode().tolist()
    # 只保留前 n 个（n = 组内行数）
    n = len(group)
    alt_list = alt_list[:n]
    # 赋值回 ALT 列
    group['ALT'] = alt_list
    return group

# 应用分组处理
df_final = df.groupby(['CHROM', 'POS', 'ID'], group_keys=False).apply(split_alt_by_group)

# 排序（可选）
# df_final = df_final.sort_values(['CHROM', 'POS']).reset_index(drop=True)

# 输出结果
print(df_final.head())
# 可选：保存
df_final.to_csv('var_transcript_data/cleaned.csv', index=False)