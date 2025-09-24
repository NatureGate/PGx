import pandas as pd
df1 = pd.read_csv('var_transcript_data/fill_na_transcript_919.csv',sep=',')
df2 = pd.read_csv('var_transcript_data/merged_all.csv',sep=',')
df3 = pd.read_csv('var_transcript_data/na_transcript.csv',sep=',')

mask1 = df1['Name'].notna()
mask2 = df2['Name'].notna()
mask3 = df3['Name'].notna()

df1 = df1[mask1]
df2 = df2[mask2]    
df3 = df3[mask3]

big_table = pd.concat([df1,df2,df3]).reset_index(drop=True)[['CHROM','POS','ID','REF','ALT','Name']]

big_table = big_table[big_table['Name'].str.contains(':')]
big_table.to_csv('var_transcript_data/combin_all.csv', index=False)
