import pandas as pd 
df = pd.read_table('new_position.txt')
df['INFO'] = df['INFO'].str.replace('PX=', '')
# old_positon = pd.read_table('drug_pots.bed')
df['chr'] = df['CHROM']+':'+df['POS'].astype(str)+':'+df['REF']+'>'+df['ALT']
df.to_csv('position.tsv', sep='\t', index=False)
# df_new = pd.merge(df, old_positon, left_on='chr', right_on='chr', how='left')
# mask = df_new['rs'].isna()
# df_new = df_new[mask]
# df_new.to_excel('new_position.xlsx', index=False)
## new_position.xlsx是新增的位点
# new_positions = pd.read_excel('new_position.xlsx')
# gene_type_table = pd.read_excel('gene_var/all_variants.xlsx')
# new_positions = pd.merge(new_positions, gene_type_table, left_on='chr', right_on='Variant', how='left')
# new_positions.to_excel('new_position_with_gene_type.xlsx', index=False)
# relationships = pd.read_table('relationships.tsv')
# print(relationships.shape)
# mask = ((relationships['Entity1_type']=='Variant') & (relationships['Entity2_type']=='Chemical')) | \
#        ((relationships['Entity1_type']=='Chemical') & (relationships['Entity2_type']=='Variant')|
#         (relationships['Entity1_type']=='Haplotype') & (relationships['Entity2_type']=='Chemical')|
#        (relationships['Entity1_type']=='Chemical') & (relationships['Entity2_type']=='Haplotype'))
# relationships = relationships[mask]  
# print(relationships.shape)    
import pandas as pd

# 假设原始 DataFrame 叫 df
# df = pd.read_csv('xxx.tsv', sep='\t')

# 1. 先拆出 Entity2 是 Chemical 的行
# var_side = (relationships[relationships['Entity2_type'] == 'Chemical']
#             .rename(columns={'Entity1_id':   'variant_id',
#                              'Entity1_name': 'variant_name',
#                              'Entity2_id':   'chemical_id',
#                              'Entity2_name': 'chemical_name'}))

# 2. 再拆出 Entity1 是 variant 的行（也就是 chemical 在左边）
# chem_side = (relationships[relationships['Entity1_type'] == 'Chemical']
#              .rename(columns={'Entity2_id':   'variant_id',
#                               'Entity2_name': 'variant_name',
#                               'Entity1_id':   'chemical_id',
#                               'Entity1_name': 'chemical_name'}))

# 3. 纵向合并，得到最终表
# new_relationships = (pd.concat([var_side, chem_side], ignore_index=True)
#          .loc[:, ['variant_id', 'variant_name',
#                   'chemical_id', 'chemical_name']])
# new_relationships.to_excel('new_relationships.xlsx', index=False)
# data = pd.merge(new_relationships, new_positions, left_on='variant_name', right_on='ID', how='right')

# data.to_excel('final_variants.xlsx', index=False)
# mask = data['variant_id'].isna()
# data_new = data[~mask]
## relationship联合drug
# drugs = pd.read_table('reporter/all_drugs.txt')
# drugs_merged_variants = pd.merge(drugs, new_relationships, left_on='drug', right_on='chemical_name', how='left')
# drugs_merged_variants.to_excel('drugs_merged_variants.xlsx', index=False)
# ## 将单倍型改成相关的rsID
# rsmask = drugs_merged_variants['variant_name'].str.startswith('rs')
# rsdrugs_merged_variants = drugs_merged_variants[rsmask]
# hadrugs_merged_variants = drugs_merged_variants[~rsmask]
# hadrugs_merged_variants['variant_name'] = hadrugs_merged_variants['variant_name'].map(gene_type_table.set_index('Haplotype')['rsID']).fillna(hadrugs_merged_variants['variant_name'])

# var_drug_ann = pd.read_table('var_drug_ann.tsv')[['Variant/Haplotypes', 'Gene']]
# drugs_merged_variants_gene = pd.merge(drugs_merged_variants, var_drug_ann, left_on='variant_name', right_on='Variant/Haplotypes', how='left').sort_values(by=['Gene'])
# mask = drugs_merged_variants_gene['Gene'].isin(['ABCG2',
# 'CACNA1S',
# 'CYP2B6',
# 'CYP2C19',
# 'CYP2C9',
# 'CYP2D6',
# 'CYP3A4',
# 'CYP3A5',
# 'CYP4F2',
# 'CFTR',
# 'DPYD',
# 'G6PD',
# 'HLA-A',
# 'HLA-B',
# 'IFNL3',
# 'MT-RNR1',
# 'NUDT15',
# 'RYR1',
# 'SLCO1B1',
# 'TPMT',
# 'UGT1A1',
# 'VKORC1',
# ]) | drugs_merged_variants_gene['Gene'].isna()
# drugs_merged_variants_gene= drugs_merged_variants_gene[mask]

# drugs_merged_variants_gene.to_excel('drugs_merged_variants_gene.xlsx', index=False)
# results = pd.merge(new_positions, drugs_merged_variants_gene, left_on='ID', right_on='variant_name', how='left')
