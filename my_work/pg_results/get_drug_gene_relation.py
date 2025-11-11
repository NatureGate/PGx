import pandas as pd 
# drugs = pd.read_excel('pg_results/old_drugs.xlsx')
# drug_gene_relation = pd.read_table('pg_results/drug_pots.bed', sep='\t')[['drug_id', 'gene']]
# merged_data = pd.merge(drugs,drug_gene_relation, left_on='id', right_on='drug_id', how='left')
# merged_data.to_excel('pg_results/drug_gene_relation_with_old_drugs.xlsx', index=False)
# merged_data = merged_data.drop_duplicates()
# merged_data.to_excel('pg_results/drug_gene_relation_with_old_drugs.xlsx', index=False)
relations  = pd.read_excel('pg_results/drug_gene_relation_with_old_drugs.xlsx')

cpic = pd.read_excel('pg_results/cpic_gene-drug_pairs.xlsx')
# 合并两个数据框，基于gene和drug name字段
relations['drugEN'] = relations['drugEN'].str.lower()
cpic['Drug'] = cpic['Drug'].str.lower()
merged = pd.merge(
    relations,
    cpic,
    left_on=['drugEN'],
    right_on=['Drug'],
    how='left'
)
merged.to_excel('pg_results/drug_gene_relation_with_cpic.xlsx', index=False)