import pandas as pd

df = pd.read_csv('reporter/merged_drugs.csv')
data = df.loc[df['id'].isna(), :]
data[['name', 'atc_code']].to_csv('reporter/additional_drugs.csv', index=False)
