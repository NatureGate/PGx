import json
import os
import sys
import argparse
import pandas as pd
import numpy as np

COLUMNS = [
    "sample", "drug", "gatkscore", "depth", "readsratio", "gene", "diplotype",
    "cHGVS", "pHGVS", "zyg", "guide", "effect", "advice", "suggest", "ref_guide",
    "drugid", "rsID", "location"
]

origin_col = ['sample','drug','QUAL','depth','readsratio','gene','diplotype','cHGVS','pHGVS','zyg','guide','implications','phenotype','recommendation','ref_guide','drug_id','rsid','location']

DRUG_ANNOTATION = {
    'CPIC':'CPIC Guideline Annotation',
    'DPWG':'DPWG Guideline Annotation',
    'FDA':'FDA Label Annotation',
    'FDA-PGx':'FDA PGx Association',
}

OUTCALL_GENE={'CYP2D6','HLA-A','HLA-B','MT-RNR1'}

DRUG_NAME_JSON = '/home/stereonote/data/drugs_union.json'
# DRUG_NAME_JSON = 'drugs_union.json'
def drug_name_en2zh(df):

    # 读取 JSON 文件
    with open(DRUG_NAME_JSON, 'r', encoding='utf-8') as f:
        drug_dict = json.load(f)

    # 替换药物名称为中文，空值保持不变（或设为空字符串）
    df['drug'] = (
        df['drug']
        .map(lambda x: drug_dict.get(x, x) if pd.notnull(x) else x)  # 先保留 NaN
        .fillna('')                                                   # 如需空字符串
    )
    return df
def reverse_string(s):
    # 将字符串按斜杠分割为两部分
    parts = s.split('/')
    if len(parts) != 2:
        return "输入格式错误"
    left_part, right_part = parts[0], parts[1]
    return f"{right_part}/{left_part}"


def prepare_variant_msg(file1, file2):
    print(f'file1, file2:{file1}, {file2}')
    df = pd.read_csv(file1,sep='\t',dtype=str)
    df['depth'] = df.iloc[:, -1].str.split(':').str[1]
    split_depth = df['depth'].str.split(',', expand=True).astype(float)
    split_depth.fillna(0, inplace=True)
# 计算 readsratio（保留两位小数）
    denominator = split_depth[0] + split_depth[1]
    df['readsratio'] = np.where(
        denominator == 0,  # 分母为0时返回0
        0,
        np.round(split_depth[1] / denominator, 2)  # 保留两位小数
    )
    df['readsratio'] = df['readsratio'].astype(str)
    df['depth'] = df['depth'].astype(str)

# 添加 location 列
    df['location'] = (
        df['CHROM'] + ':' + 
        df['POS'].astype(str) + ':' + 
        df['REF'] + '>' + 
        df['ALT']
    )
    prepharmcat_data = pd.read_csv(file2, sep='\t',dtype=str)
    prepharmcat_data['PX'] = prepharmcat_data['INFO'].str.extract(r'PX=([^;]+)')
    prepharmcat_data['cHGVS'] = prepharmcat_data['INFO'].str.split('|').str[9]
    prepharmcat_data['pHGVS'] = prepharmcat_data['INFO'].str.split('|').str[10]
    df['gene'] = prepharmcat_data['PX']
    df['rsid'] = prepharmcat_data['ID']
    df['cHGVS'] = prepharmcat_data['cHGVS']
    df['pHGVS'] = prepharmcat_data['pHGVS']
    grouped_df = df.groupby('gene', as_index=False).agg({
        'cHGVS': lambda x: ','.join(set(x)),  
        'pHGVS': lambda x: ','.join(set(x)),  
        'rsid': lambda x: ','.join(set(x)),  
        'location': lambda x: ','.join(set(x)),  
        'cHGVS': lambda x: ','.join(set(x)),
        'pHGVS': lambda x: ','.join(set(x)),
        'QUAL': lambda x: ','.join(set(x)),
        'depth': lambda x: ';'.join(set(x)),
        'readsratio': lambda x: ';'.join(set(x)),
 
    })
    print(f"grouped_df:{grouped_df}")
    grouped_df.to_csv('grouped_df.csv',index=False)
    return grouped_df

def get_gene(df:pd.DataFrame, outcall_file):
    
    
      
    columns = df.columns.tolist()
    if os.path.getsize(outcall_file) != 0:
        print("outcall_file is not empty ")
        outcallgene = set(pd.read_table(outcall_file, sep='\t', header=None)[0])
        print(outcallgene-set(df['gene']))
        df = df[~df['gene'].isin(outcallgene-set(df['gene']))]
        len(outcallgene)
        outcall_df = pd.DataFrame(columns=columns, index=range(len(outcallgene-set(df['gene']))))
        outcall_df['gene'] = list(outcallgene-set(df['gene']))  
        df = pd.concat([df,outcall_df]).reset_index(drop=True)
        print(f"outcall_df:{outcall_df}['gene']")
    
    return df

def add_drugs(df,related_drugs, gene,diplotype,annotation_type):
    
    for drugs in related_drugs:
        # if drugs['name'].__contains__('/'):
        #     # for drug in drugs['name'].split('/'):
        #     #     new_row=[drug,drugs['id'], gene,diplotype,annotation_type]
        #     #     df.loc[len(df)] = new_row  
        #     # # new_row=[drugs['name'], drugs['id'], gene,diplotype,annotation_type]
        #     # # df.loc[len(df)] = new_row
        #     continue
        # else:
        new_row=[drugs['name'], drugs['id'], gene,diplotype,annotation_type]
        df.loc[len(df)] = new_row  
    return df

def get_all_guideline(df:pd.DataFrame,data,annotation_type,diplotype,gene_table):
    _guide = data['drugs'][annotation_type]
    df = df[(df['drug'].notna()) & (df['drug'] != '')]
    all_drugs = df['drug'].unique()
    all_guidelines = []
    for drug_name in all_drugs:
        drug_name = drug_name.strip()
        if drug_name not in _guide:
            print(f"Warning: {drug_name} not found in CPIC Guideline Annotation")
            continue
        guidelines = _guide[drug_name]['guidelines']
        drug_id = _guide[drug_name]['id']
        for guideline_item in guidelines:
            annotations = guideline_item['annotations'] if 'annotations' in guideline_item.keys() else []
            if len(annotations) == 0:
                continue
            all_guidelines = get_annotation_guideline(all_guidelines, drug_name, annotations,annotation_type,drug_id,diplotype,gene_table)
    # print(f"all_guidelines:{all_guidelines}")

    return all_guidelines

def get_annotation_guideline(all_guidelines, drug_name, annotations,annotation_type,drug_id,diplotype,gene_table):
    for annotation in annotations:
        implications = ';'.join(annotation['implications']) if 'implications' in annotation.keys() else ''
        drug_recommendation = annotation['drugRecommendation'] if 'drugRecommendation' in annotation.keys() else 'Not Found'
        for genotype in annotation['genotypes']:
                    # diplotypes = list(filter(lambda x: x['allele1']['gene'] == gene, genotype['diplotypes']))
            diplotypes = genotype['diplotypes']
            if len(diplotypes) == 0:
                
                continue
            all_guidelines = add_guide_line(all_guidelines, drug_name, implications, drug_recommendation, diplotypes,annotation_type,drug_id,gene_table)
    return all_guidelines

# 在这里排查gene与type是否匹配
# 如果gene与type不匹配，则不添加到all_guidelines中
def add_guide_line(all_guidelines, drug_name, implications, drug_recommendation, diplotypes,annotation_type,drug_id,gene_table):
    for diplotype in diplotypes:
        di_gene = diplotype['gene']
        if di_gene=='MT-RNR1':
            print(f"di_gene:{di_gene},diplotype:{diplotype}")
        gene_type = diplotype['label']
        filtered_gene_table = gene_table[~gene_table['Source Diplotype'].str.contains(' OR | AND ', na=False)]
        mask = (gene_table['Gene'] == di_gene) & gene_table['Source Diplotype'].str.contains('OR|AND', case=False, na=False)
        if di_gene not in OUTCALL_GENE and gene_table[(gene_table['Gene'] == di_gene)].empty:
            continue
        if di_gene not in OUTCALL_GENE and filtered_gene_table[(filtered_gene_table['Gene'] == di_gene)&(filtered_gene_table['Source Diplotype'] == gene_type)].empty:
            if gene_type.__contains__(' (heterozygous)'):
                gene_type = gene_type[:gene_type.index(' (heterozygous)')]
                var_type = filtered_gene_table[(filtered_gene_table['Gene'] == di_gene)]['Source Diplotype'].values[0]
                if var_type.__contains__('Reference/'+gene_type):
                    # print(f'gene,type,{di_gene},{gene_type}')
                    gene_type = var_type
                else:
                    continue
            elif len(gene_table[mask]) > 0 :
                print(f'gene:{di_gene},type:{gene_type}')
            else:
                print(f"Warning: {di_gene} with type {gene_type} not found in gene table")
                continue
            #print(f"di_gene,type:{di_gene},{type}")
            # continue
        
        di_implications = implications
        di_recommendation = drug_recommendation
        
        all_guidelines.append([
            drug_name, di_gene, gene_type, di_implications, di_recommendation, annotation_type, drug_id
        ])
    return all_guidelines

def get_report(report_json_file, var_gene,gene_table):
    columns = ['drug', 'drug_id','gene','diplotype','ref_guide']
    df1 = pd.DataFrame(columns=columns)
    df2 = pd.DataFrame(columns=columns)
    df3 = pd.DataFrame(columns=columns)
    df4 = pd.DataFrame(columns=columns)
    with open(report_json_file, 'r',encoding='utf-8') as file:
        data = json.load(file)
        cpic_all_guides:list = get_annotation(var_gene, df1, data,DRUG_ANNOTATION['CPIC'],gene_table)
        dpwg_all_guides:list = get_annotation(var_gene, df2, data,DRUG_ANNOTATION['DPWG'],gene_table)
        fda_all_guides:list = get_annotation(var_gene, df3, data,DRUG_ANNOTATION['FDA'],gene_table)
        fda_pgx_all_guides:list = get_annotation(var_gene, df4, data,DRUG_ANNOTATION['FDA-PGx'],gene_table)
        all_guides = (
            cpic_all_guides + 
            dpwg_all_guides + 
            fda_all_guides + 
            fda_pgx_all_guides
        )
        df = pd.DataFrame(all_guides, columns=['drug', 'gene', 'diplotype', 'implications', 'recommendation','ref_guide','drug_id'])
        df.to_csv('all_guidelines.csv', index=False, encoding='utf-8')
        return df

def get_annotation(diplotype_gene_set, df, data, annotation_type,gene_table):

    for gene in diplotype_gene_set:
        
        if gene in data['genes']['CPIC'].keys():
            
                
            related_drugs = data['genes']['CPIC'][gene]['relatedDrugs']
            # diplotype = data['genes']['CPIC'][gene]['sourceDiplotypes'][0]['allele1']['name'] \
            #                 + '/' + data['genes']['CPIC'][gene]['sourceDiplotypes'][0]['allele2']['name']
            diplotype = data['genes']['CPIC'][gene]['sourceDiplotypes'][0]['label']
            df = add_drugs(df,related_drugs, gene,diplotype,annotation_type)
        else:
            print(f"Warning: {gene} not found in CPIC Guideline Annotation")
            new_row = ['', '', gene, '', annotation_type]
            df.loc[len(df)] = new_row
    
    results = get_all_guideline(df,data,annotation_type,diplotype,gene_table)
    
    return results

def get_phenotype(report_json_file, diplotype_gene_set, results_table):
    results_table['phenotype'] = ''
    with open(report_json_file, 'r',encoding='utf-8') as file:
        data = json.load(file)
        for gene in diplotype_gene_set:
            if gene in data['genes']['CPIC']:
                #gene_all_diplotypes = data['genes']['CPIC'][gene]['sourceDiplotypes']
                gene_all_diplotypes = data['genes']['CPIC'][gene]['sourceDiplotypes']
                phenotype = ';'.join(set(map(lambda x: x['phenotypes'][0] if len(x['phenotypes']) > 0 else ''   , gene_all_diplotypes)))
                #phenotype = ';'.join(gene_all_diplotypes)
                results_table.loc[results_table['gene'] == gene, 'phenotype'] = phenotype
    return results_table        

def add_filter_gene(results_table,filtered_df,sample_id):
    
    filtered_df = filtered_df.rename(columns={
        'QUAL': 'gatkscore',
        'depth': 'depth',
        'readsratio': 'readsratio',
        'gene': 'gene',
        'cHGVS': 'cHGVS',
        'pHGVS': 'pHGVS',
        'rsid': 'rsID',
        'location': 'location'
    })
    filtered_df['drug'] = ''
    filtered_df['diplotype'] = ''
    filtered_df['zyg'] = ''  
    filtered_df['guide'] = ''   
    filtered_df['effect'] = ''  
    filtered_df['advice'] = ''
    filtered_df['suggest'] = ''
    filtered_df['ref_guide'] = ''
    filtered_df['drugid'] = ''
    filtered_df['sample'] = sample_id
    filtered_df = filtered_df[COLUMNS]
    results_table = pd.concat([results_table, filtered_df], ignore_index=True)      
    return results_table  

# 这里添加的是有CPIC用药指南，但是没有检测出变异的基因。
def add_reference_gene_report(results_table,var_gene):
    
    data = pd.read_csv('/home/stereonote/data/merged_output.csv',sep='\t',dtype=str)
    # data = pd.read_csv('merged_output.csv',sep='\t',dtype=str)
    default_gene = set(data['gene'])
    results_gene = var_gene
    print(f"results_gene:{results_gene}")
    default_gene = default_gene - results_gene
    print(f"default_gene:{default_gene}")
    print(f"data:{data.columns}")
    print(f"results_table:{results_table.columns}")
    print(f"results_table['sample']:{results_table['sample']}")
    data['sample'] = results_table['sample'].iloc[0]
    results_table.astype(str)
    data = data[results_table.columns]
    results_table.astype(str)
    data = data[data['gene'].isin(default_gene)]
    results_table = pd.concat([results_table, data], ignore_index=True)
    print(f"results_table:{results_table}")
    results_table = results_table.drop_duplicates(subset=['gene', 'drug', 'diplotype'], keep='first')
    results_table = results_table.reset_index(drop=True)
    print(f"results_table after drop duplicates:{results_table}")
    results_table.to_csv('results_table.csv', index=False, encoding='utf-8')
    return results_table



def filter_single_effect_genes(unknown_type_gene_table: pd.DataFrame,
                               gene_col: str = 'gene',
                               effect_col: str = 'effect',
                               group_cols: list = None) -> pd.DataFrame:
    """
    筛选gene diplotype不明的数据

    参数
    ----
    df : pd.DataFrame
        输入的 DataFrame
    gene_col : str, 默认 'gene'
        gene 列名
    effect_col : str, 默认 'effect'
        effect 列名
    group_cols : list[str], 默认 ['gene', 'drug']
        用于分组的列名列表

    返回
    ----
    pd.DataFrame
        已分组并去重后的结果
    """
    if group_cols is None:
        group_cols = ['gene', 'drug']
    # 过滤 effect 不含分号的行
    mask = unknown_type_gene_table[effect_col].astype(str).str.contains(';', na=False)
    unknown_type_gene_table = unknown_type_gene_table[~mask]

    # 按指定列分组，每组保留第一行
    result = unknown_type_gene_table.groupby(group_cols, as_index=False).first()

    return result

def split_drug(df):
    # 2. 拆分 drug 列中含 "/" 的行
    has_slash = df['drug'].str.contains('/', na=False)
    df_to_split = df[has_slash].copy()
    df_intact   = df[~has_slash].copy()
    # 3. 真正拆分
    #    先按 / 拆成多列，再 stack 成 Series，再还原成 DataFrame
    split_series = (
        df_to_split['drug']
        .str.split('/', expand=True)   # 拆成多列
        .stack()                       # 变成长 Series
        .reset_index(level=1, drop=True)
        .rename('drug')                # 保持列名一致
    )

    # 4. 把拆分后的 drug 拼回原 DataFrame
    df_split = (
        df_to_split
        .drop(columns=['drug'])
        .join(split_series)
        .reset_index(drop=True)
    )

    # 5. 合并结果
    result = pd.concat([df_intact, df_split], ignore_index=True)    

    return result

def main():
    # describe the tool
    parser = argparse.ArgumentParser(description='Extract results from PharmCAT JSONs into a TSV-formatted file.')
    # list input arguments
    parser.add_argument("-report_json_file", type=str,help="Pharmcat report result file",default='data/E-B24925859205.report.json')
    parser.add_argument("-result_tsv_file", type=str,help="Pharmcat result tsv file",default='data/E-B24925859205.report.tsv')
    parser.add_argument("-outcall_file", type=str,help="pharcat outcall file",default='data/E-B24925859205.outcall.tsv')
    parser.add_argument("-intersect_file", type=str,help="pharmcat pipeline inputfile",default='data/intersect.csv')
    parser.add_argument("-prepharmcat_file", type=str,help="",default='data/annotated.csv')
    parser.add_argument("-sample_id", type=str,help="",default='Test')
    args = parser.parse_args()
    print(f"args:{args}")
    
    df = prepare_variant_msg(args.intersect_file, args.prepharmcat_file)
    vcf_gene_set = set(df['gene'])
    vcf_df = df
    print(f'vcf_gene_set:{vcf_gene_set}')
    gene_table = pd.read_table(args.result_tsv_file,sep='\t',skiprows=1)
    unknown_type_gene = list(gene_table[gene_table['Source Diplotype'].str.contains(' OR | AND ', na=False)]['Gene'])
    unknown_type_gene = set(unknown_type_gene)-OUTCALL_GENE
    print(f'unknown_type_gene:{unknown_type_gene}')
    var_gene = set(gene_table['Gene'])|OUTCALL_GENE
    df = get_gene(df, args.outcall_file)
    
    df.to_csv('variant_msg.csv', index=False, encoding='utf-8')
    diplotype_gene_set = set(df['gene'])
    print(f'diplotype_gene_set:{diplotype_gene_set}')

    drug_gene_table = get_report(args.report_json_file, var_gene,gene_table)
    results_table = drug_gene_table.merge(df,on='gene', how='left')
    results_table = results_table.groupby(['gene','drug','diplotype','ref_guide'], as_index=False).first()
    results_table.to_csv('results_table.csv', index=False, encoding='utf-8')
    mask = (results_table['diplotype'] == 'Unknown/Unknown') | (results_table['diplotype'] == 'Not Found')
    results_table = results_table[~mask]
    results_table = get_phenotype(args.report_json_file,diplotype_gene_set,results_table)
    # 创建mask筛选出包含斜杠的行
    mask = results_table['diplotype'].str.contains('/')

    # 对符合条件的行应用lambda函数计算zyg列
    results_table.loc[mask, 'zyg'] = results_table.loc[mask, 'diplotype'].apply(
        lambda x: 'Hom' if x.split('/')[0] == x.split('/')[1] else 'Het'
    )
    # results_table['zyg'] = results_table['diplotype'].apply(lambda x: 'Hom' if x.split('/')[0] == x.split('/')[1] else 'Het')
    results_table['sample']= args.sample_id  
    results_table['guide'] = '' 
    results_table = results_table[origin_col]
    results_table = results_table.rename(columns={
        'sample': 'sample',
        'drug': 'drug',
        'QUAL': 'gatkscore',
        'depth': 'depth',
        'readsratio': 'readsratio',
        'gene': 'gene',
        'diplotype': 'diplotype',
        'cHGVS': 'cHGVS',
        'pHGVS': 'pHGVS',
        'zyg': 'zyg',
        'phenotype':'effect',
        'guide': 'guide',
        'implications': 'advice',
        'recommendation': 'suggest',
        # 'classification': 'suggest',
        'ref_guide': 'ref_guide',
        'drug_id': 'drugid',
        'rsid': 'rsID',
        'location': 'location'
    })
#     result_df = results_table.groupby(['gene', 'drug', 'diplotype']).apply(lambda group: group[group['ref_guide'].str.startswith('CPIC')]
# ).reset_index(drop=True)
    
    results_table = add_reference_gene_report(results_table,var_gene)
    results_table = results_table.fillna('.').replace('Not Found', '.')
    results_table['sample'] = args.sample_id
    results_table = split_drug(results_table)
    selected_genes = results_table.groupby('gene')['diplotype'].nunique()[lambda x: x > 1].index
    unknown_type_gene_table = results_table[results_table['gene'].isin(selected_genes)]
    known_type_gene_table = results_table[~results_table['gene'].isin(selected_genes)]
    known_type_gene_table = drug_name_en2zh(known_type_gene_table)
    if len(unknown_type_gene_table)>0:
        unknown_type_gene_table = drug_name_en2zh(unknown_type_gene_table)
        unknown_type_gene_table.to_csv(f'{args.sample_id}_unknown.pgx.tsv', index=False, encoding='utf-8', sep='\t')
        add_diplotype_table = filter_single_effect_genes(unknown_type_gene_table)
        known_type_gene_table = pd.concat([known_type_gene_table, add_diplotype_table], ignore_index=True) if 'add_diplotype_table' in locals() else known_type_gene_table
    known_type_gene_table.to_csv(f'{args.sample_id}.pgx.tsv', index=False, encoding='utf-8', sep='\t')
    # results_table = add_filter_gene(results_table,filtered_df,args.sample_id)
    # results_table.to_csv(f'{args.sample_id}.pgx.tsv', index=False, encoding='utf-8',sep='\t')



if __name__ == "__main__":
    main()
