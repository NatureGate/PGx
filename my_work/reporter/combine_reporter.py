import pandas as pd
df = pd.read_excel('reporter/guideline.xlsx')
df.head()
df2 = pd.read_excel('reporter/merged_short_1208.xlsx')
merged_data = pd.merge(df2,df,on=['Text','Related Chemicals'],how='left')
merged_data.to_excel('reporter/merged_1209.xlsx', index=False)
