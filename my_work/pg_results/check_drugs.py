import pandas as pd
df = pd.read_excel('pg_results/药物核对-1107-zw.xlsx')
print(df)
eb_df = pd.read_table('pg_results/E-B24169320585.drug.tsv', sep=',')
print(eb_df)
eb_df['name_en'] = eb_df['name_en'].str.capitalize()
df['drugEN'] = df['drugEN'].str.capitalize()
df.to_excel('pg_results/drug_details_1107_zw.xlsx', index=False)
drug_details = pd.merge(
    eb_df,
    df,
    left_on='name_en',
    right_on='drugEN',
    how='left'
)

# drug_details.to_excel('pg_results/drug_details.xlsx', index=False)
# 找出drug与drugEN不相同的数据
diff_mask = drug_details['drug'] != drug_details['drugCN']
diff_drug = drug_details[diff_mask]
print('drug与drugEN不相同的数据：')
print(diff_drug[['name_en', 'drug', 'drugCN']])
print('drug与drugCN不相同的数量：', len(diff_drug))
mask = drug_details['drugEN'].isna()
mask_drug = drug_details[mask]
print('新的药物：')
print(mask_drug['name_en'], mask_drug['drug'])
print('新的药物数量：', len(mask_drug))

