import json
import pandas as pd

# 这里假定guideline.json是你的json文件``
# with open('reporter/prescribing_guidance.json','r',encoding='utf-8') as f:
#     data = json.load(f)

# f = open('guideline_positions_new.tsv', 'w', encoding='utf-8')
# f.write("Chemical\tGene\tAllele\tGene_number\n")
# for guideline in data['guidelines']:
#     guideline_entity = guideline['guideline']
#     related_chemicals = guideline_entity['relatedChemicals']
#     related_genes = guideline_entity['relatedGenes']
#     related_alleles = guideline_entity['relatedAlleles']
#     for chemical in related_chemicals:
#         chemical_name = chemical['name']
#         if chemical_name == 'quetiapine':
#             print('debug')
#         for gene in related_genes:
#             gene_name = gene['symbol']
            
#             if len(related_alleles)>0:
#                 for allele in related_alleles:
#                     allele_name = allele['symbol']
#                     f.write(f"{chemical_name}\t{gene_name}\t{allele_name}\t{len(related_genes)}\n")
#             else:
#                 f.write(f"{chemical_name}\t{gene_name}\tall alleles\t{len(related_genes)}\n")
df = pd.read_table('guideline_positions_new.tsv', sep='\t')
print(set(df['Gene']))
grouped_df = df.groupby(['Chemical', 'Gene','Allele']).first().reset_index()
mask = grouped_df[(grouped_df['Gene_number'] > 1)&(grouped_df['Allele'] != 'all alleles')].index
multi_gene_df = grouped_df.loc[mask]
# multi_gene_df.to_excel('multi_gene_df.xlsx', index=False)
import pandas as pd

# 1. 先建好“基因 → 合法 allele 集合”
legal = {
    'RYR1':    {'rs1800559', 'rs772226819'},   # 对 RYR1 来说这两个是 **非法**
    'CACNA1S': {'rs1800559', 'rs772226819'}    # 对 CACNA1S 来说这两个是 **合法**
}

# 2. 定义一个小函数：只处理 RYR1 和 CACNA1S，其余基因直接放行
def fix_allele(row):
    g, a = row['Gene'], row['Allele']
    if g not in legal:          # 非目标基因，一律放行
        return a
    if g == 'RYR1':             # RYR1 不允许出现在 legal[g] 里
        return a if a not in legal[g] else pd.NA
    if g == 'CACNA1S':          # CACNA1S 只允许出现在 legal[g] 里
        return a if a in legal[g] else pd.NA
    return a

# 3. 应用函数并一次性丢掉非法行（或保留 NA 后续再填）
multi_gene_df['Allele'] = multi_gene_df.apply(fix_allele, axis=1)
clean = multi_gene_df.dropna(subset=['Allele'])          # 去掉被标记为 NA 的非法行
# 手动修改
# clean.to_excel('multi_gene_df.xlsx', index=False)
all_mask = grouped_df[grouped_df['Allele'] == 'all alleles'].index
all_allele_df = grouped_df.loc[all_mask]

all_var = pd.read_excel('gene_var/all_variants.xlsx')
all_merged = pd.merge(all_allele_df, all_var,left_on='Gene', right_on='Gene', how='left').groupby(['Chemical', 'Gene','Variant']).first().reset_index()
all_merged.to_excel('all_merged.xlsx', index=False)

################ 处理multi_gene_df.xlsx ################
data = pd.read_excel('multi_gene_df.xlsx')
exp = all_var.assign(rsID=all_var['rsID']).explode('rsID')
print(exp)
mapper = exp.rename(columns={'Haplotype': 'old_rs', 'rsID': 'new_rs'})
multi_allele_table = (data
       .merge(mapper, left_on='Allele', right_on='old_rs', how='left')
       .assign(rsid=lambda x: x['new_rs'].fillna(x['Allele']))  # 没匹配到的保留原值
       .drop(columns=['old_rs', 'new_rs'])
      )
multi_allele_table['Allele'] = multi_allele_table['rsid']
multi_allele_table = multi_allele_table[['Chemical', 'Gene_x', 'Allele', 'Gene_number']] 
multi_allele_table.to_excel('multi_gene_fixed.xlsx', index=False)

pd.merge(multi_allele_table, all_var, left_on=['Gene_x','Allele'], right_on=['Gene','rsID'], how='left').groupby(['Chemical', 'Gene','Variant']).first().reset_index().to_excel('multi_gene_with_var.xlsx', index=False)

# #################处理单个的情况###############
mask = grouped_df[(grouped_df['Gene_number'] == 1)&(grouped_df['Allele'] != 'all alleles')].index
single_gene_df = grouped_df.loc[mask]
single_gene_df.to_excel('single_gene_df.xlsx', index=False)
exp = all_var.assign(rsID=all_var['rsID']).explode('rsID')
mapper = exp.rename(columns={'Haplotype': 'old_rs', 'rsID': 'new_rs'})
single_gene_table = (single_gene_df
       .merge(mapper, left_on='Allele', right_on='old_rs', how='left')
       .assign(rsid=lambda x: x['new_rs'].fillna(x['Allele']))  # 没匹配到的保留原值
       .drop(columns=['old_rs', 'new_rs'])
      )
single_gene_table['Allele'] = single_gene_table['rsid']
single_gene_table = single_gene_table[['Chemical', 'Gene_x', 'Allele', 'Gene_number']]

pd.merge(single_gene_table, all_var, left_on=['Gene_x','Allele'], right_on=['Gene','rsID'], how='left').groupby(['Chemical', 'Gene','Variant']).first().reset_index().to_excel('single_gene_with_var.xlsx', index=False)

####################################################
all_drug_gene_var = pd.read_excel('all_drug_gene_var.xlsx')
positions = pd.read_csv('position.tsv', sep='\t')
final = pd.merge(all_drug_gene_var, positions, left_on=['Variant'], right_on=['chr'], how='right')
final.to_excel('final_guideline_positions.xlsx', index=False)

# rsfinal = pd.merge(all_drug_gene_var, positions, left_on=['rsID'], right_on=['ID'], how='right')
# rsfinal.to_excel('final_guideline_positions_rs.xlsx', index=False)
miss_mask = final['Chemical'].isna()& (final['ID']!='.')
miss_match_col = final[miss_mask][positions.columns]
rs_final = pd.merge(all_drug_gene_var, miss_match_col, left_on=['rsID'], right_on=['ID'], how='right')
rs_drop_mask = rs_final['Chemical'].isna()
rs_final = rs_final[~rs_drop_mask]
final_drop_mask = final['Chemical'].isna()& (final['ID']!='.')
final = final[~final_drop_mask]

final_result = pd.concat([final, rs_final], ignore_index=True)
mask = (final_result['INFO']!='CYP3A4')&(final_result['Chemical'].isna())
missed_final_result = final_result[mask]
missed_final_result[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'PharmCAT', 'chr']].to_excel('missed_final_result.xlsx', index=False)
pd.merge(all_allele_df, missed_final_result, left_on=['Gene'], right_on=['INFO'], how='right').to_excel('missed_with_var.xlsx', index=False)
