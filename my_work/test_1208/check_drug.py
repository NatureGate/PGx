import pandas as pd
df = pd.read_excel('test_1208/药物背景汇总-新-1203.xlsx',sheet_name='药物核对')
data = pd.read_table('test_1208/2901227132C.drug.tsv',sep=',')
merge_data = pd.merge(data,df,left_on='drug',right_on='药物',how='left')
mask = merge_data['药物'].isna()
print(merge_data[mask]['drug'].tolist())