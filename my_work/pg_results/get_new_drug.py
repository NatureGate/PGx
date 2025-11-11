import os
import glob
import pandas as pd
# 设置 pandas 打印最大行数为无限制，方便调试查看完整数据
pd.set_option('display.max_rows', None)


# 读取pg-results目录下所有以E-B开头的文件
# eb_files = os.path.join('pg-results', 'E-B24169320585.drug.tsv')
eb_df = pd.read_table('pg_results/E-B24169320585.drug.tsv', sep=',')
# 读取pg-results目录下的drug_detail.txt文件
# drug_detail_path = os.path.join('pg-results', 'drug_detail.txt')
drug_detail_df = pd.read_table('pg_results/drug_detail.txt', sep='\t')
# 不考虑大小写，找到eb_df 的name_en与drug_detail_df 的drug_EN相同的数据
common_drugs = pd.merge(
    eb_df['name_en'].str.lower().to_frame('name_en'),
    drug_detail_df['drugEN'].str.lower().to_frame('drugEN'),
    left_on='name_en',
    right_on='drugEN',
    how='inner'
)['name_en'].drop_duplicates()

# 找到在eb_df中但不在drug_detail_df中的数据
new_drugs = eb_df[~eb_df['name_en'].str.lower().isin(drug_detail_df['drugEN'].str.lower())]['name_en']
print('新的药物：')
# print(new_drugs)
print('新的药物数量：', len(new_drugs))

# 读取pg_results目录下的药物核对.xlsx文件
check_df = pd.read_excel('F:\code\PGx\my_work\pg_results\药物核对.xlsx')
new_drugs_in_check = set(check_df['name'])-set(new_drugs)
print('新的药物在核对文件中不存在的：')
print(new_drugs_in_check)
print('新的药物在核对文件中不存在的数量：', len(new_drugs_in_check))
# 找到在common_drugs中但不在eb_df中的数据
old_drugs = drug_detail_df[~drug_detail_df['drugEN'].str.lower().isin(eb_df['name_en'].str.lower())]
# print('旧的药物：')
# print(old_drugs)
# print('旧的药物数量：', len(old_drugs))

# 获取old_drugs的drug_id和drugCN
# old_drugs = pd.merge(
#     old_drugs,
#     drug_detail_df,
#     left_on='drugEN',
#     right_on='drugEN',
#     how='left'
# )
# print('旧的药物：')
# print(old_drugs)
# print('旧的药物数量：', len(old_drugs))
# old_drugs.to_excel('pg_results/old_drugs.xlsx', index=False)


# 找到同时在eb_df和drug_detail_df中的数据
# common_drugs = pd.merge(
#     eb_df['name_en'].str.lower().to_frame('name_en'),
#     drug_detail_df['drugEN'].str.lower().to_frame('drugEN'),
#     left_on='name_en',
#     right_on='drugEN',
#     how='inner'
# )['name_en'].drop_duplicates()
# print('共同的药物：')
# print(common_drugs)
# print('共同的药物数量：', len(common_drugs))

