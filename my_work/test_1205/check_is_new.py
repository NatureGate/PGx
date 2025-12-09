import pandas as pd
import numpy as np
df1 = pd.read_table('test_1205/E-B25997180293.drug.xls')
print(df1.head())
drugs = df1['drug'].tolist()
df2 = pd.read_excel('test_1205/drug_names.xlsx')
df2['is_new'] = np.where(df2['name_zh'].isin(drugs), 'No', 'Yes')
print(df2.head())
df2.to_excel('test_1205/drug_names.xlsx', index=False)
