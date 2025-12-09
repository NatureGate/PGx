import pandas as pd

# 拆分REF/ALT中的逗号，并生成CHROM:POS:REF>ALT列
def split_ref_alt(row):
    chrom = row['CHROM']
    pos = row['POS']
    ref = row['REF']
    alt = row['ALT']
    # 如果ALT中有逗号，拆分成多个ALT
    alts = alt.split(',') if pd.notna(alt) else []
    results = []
    for a in alts:
        results.append(f"{chrom}:{pos}:{ref}>{a}")
    return results


# df = pd.read_table('site_1206/2901227132C.drug.tsv',sep=',')[['drug', 'gene']]
# vcf_df = pd.read_table('site_1206/pharmcat_positions.vcf',sep='\t').astype(str)
# print(df.head())
# print(vcf_df.head())

# merged_df = pd.merge(df, vcf_df, left_on='gene', right_on='INFO', how='left')
# 应用拆分函数并展开为多行
# expanded = merged_df.apply(lambda row: pd.Series(split_ref_alt(row)), axis=1).stack().reset_index(level=1, drop=True)
# expanded.name = 'variant'
# merged_df = merged_df.join(expanded)[['drug', 'gene', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'variant']]
# # merged_df.to_csv('site_1206/2901227132C.drug.variant.tsv', sep='\t', index=False)
# old_drug_pots = pd.read_table('site_1206/drug_pots_full.bed').astype(str)
# data = pd.merge(merged_df,old_drug_pots,left_on='variant',right_on='chr',how='left')
# data = data.groupby(['variant','drug'],as_index=False).first()
# # data.to_csv('site_1206/site_results.csv')
# data_add_mask = data['id'].isna()
# need_add_data = data[data_add_mask]
# old_data = data[~data_add_mask]
# add = need_add_data.groupby(['variant','drug'],as_index=False)[['variant','drug']].first()
# all_var = pd.read_excel('site_1206/all_variants_1206.xlsx')
# add_var = pd.merge(add,all_var,left_on='variant',right_on='Variant',how='left')
# add_var = add_var.groupby(['variant','drug'],as_index=False).first()
# add_var.to_excel('site_1206/add_variants_1206.xlsx', index=False)
# 
########################################################
# data = pd.read_excel('site_1206/add_variants_1207.xlsx')
# bg_data = pd.read_table('site_1206/pgx_data_bg_results.csv',sep=',').astype(str)
# need_add_indexes = data[data['Gene'].isna()].index
# print(need_add_indexes)
# for need_add_index in need_add_indexes:
    
#     variant = data.iloc[need_add_index,0]
#     print(variant)
#     pos = variant.split(':')[1]
#     var = variant.split(':')[2]
#     ref = var.split('>')[0]
#     alt = var.split('>')[1] 
#     is_contains_var = (bg_data['REF_x']==ref) & (bg_data['ALT_x']==alt)
#     is_startswith_NM = bg_data['Name'].str.startswith('NM_')
#     is_contains_position = bg_data['POS'] == pos
#     mask = is_contains_position & is_startswith_NM & is_contains_var
#     bg_mask_data = bg_data[mask]
    
#     if mask.sum() > 0:
#         data.iloc[need_add_index,3] = bg_data[mask]['ID_x'].values[0]
#         name = bg_data[mask]['Name'].values[0]
#         # NM_000771.4(CYP2C9):c.226G>A(p.Val76Met)
#         lname = name.split(':')[0]
#         rname = name.split(':')[1]
#         data.iloc[need_add_index,2] = lname[lname.index('(')+1:-2]
#         data.iloc[need_add_index,-1] = rname[:rname.index('(')]
#         data.iloc[need_add_index,-2] = rname[rname.index('(')+1:-1]
#         print(data.iloc[need_add_index,:])
# data.to_excel('site_1206/add_variants_1207_2.xlsx', index=False)

# 继续处理
data = pd.read_excel('site_1206/add_variants_1207_2.xlsx')
multi_var_indexes = data[data['GeneChange'].str.contains(';')].index
for multi_var_index in multi_var_indexes:
    gene_changes = data.iloc[multi_var_index,-1].split(';')
    pro_change = data.iloc[multi_var_index,-2].split(';')
    var = data.iloc[multi_var_index,0].split(':')[2]
    ref_map = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    trans_table = str.maketrans(ref_map)
    el_var = var.translate(trans_table)
    for idx,gene_change in enumerate(gene_changes):
        if var in gene_change:
            data.iloc[multi_var_index,-1] = gene_changes[idx]
            data.iloc[multi_var_index,-2] = pro_change[idx]
        elif el_var in gene_change:
            data.iloc[multi_var_index,-1] = gene_changes[idx]
            data.iloc[multi_var_index,-2] = pro_change[idx]

data.to_excel('site_1206/add_variants_1207_3.xlsx', index=False)
