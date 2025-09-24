import pandas as pd
# df1 = pd.read_excel('pgx_bgdata/final_guideline_positions_all_filter.xlsx')[['CHROM','POS','ID','REF','ALT','Gene','Chemical']].astype(str)
# df2 = pd.read_excel('pgx_bgdata/missed_with_var_filter.xlsx')[['CHROM','POS','ID','REF','ALT','Gene','Chemical']].astype(str)
# data = pd.concat([df1,df2])
# print(data.shape)

# ann_data = pd.read_csv('pgx_bgdata/merge_pos_var_ann.csv', dtype=str)
# results_data = pd.merge(ann_data, data, on=['POS'], how='left')
# results_data.to_excel('pgx_bgdata/results.xlsx', index=False)
# # 3. 核心逻辑：Name 列是否包含 Gene 内容
# def fix_name(row):
#     name, gene = row['Name'], row['Gene']
#     if pd.isna(gene) or gene == '':
#         return name                       # 没有 Gene 信息，保持原样
#     if str(gene).strip() in str(name):
#         return name                       # 已包含，不动
#     # 不包含，按 ":" 拆分取第 0 段，再拼 "(Gene)"
#     prefix = str(name).split(':')[0]
#     position = str(name).split(':')[1] if ':' in str(name) else ''
#     return f'{prefix}({gene}):{position}'

# results_data['Name'] = results_data.apply(fix_name, axis=1)

# results_data.to_csv('pgx_bgdata/pgx_data.csv', index=False)
# data = pd.read_csv('pgx_bgdata/pgx_data.csv', dtype=str)
# mask = (data['Gene'].notna() & (data['Gene'] != '')&(~data['Chemical'].str.contains('/'))&(~data['Chemical'].str.contains(' and ')))
# data[mask].to_csv('pgx_bgdata/pgx_data_bg_results.csv', index=False)

data = pd.read_csv('pgx_bgdata/pgx_data_bg_results.csv', dtype=str)
mask = (~data['Chemical'].str.contains(' and '))&(~data['Chemical'].str.contains('/'))
clean_data = data[mask]
clean_data = clean_data[clean_data['ALT_x'].str.contains(',')]
clean_data.to_csv('pgx_bgdata/pgx_data_bg_results_clean.csv', index=False)