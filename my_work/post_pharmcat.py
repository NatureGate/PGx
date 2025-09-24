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

OUTCALL_GENE={'CYP2D6','HLA-A','HLA-B'}


def reverse_string(s):
    # 将字符串按斜杠分割为两部分
    parts = s.split('/')
    if len(parts) != 2:
        return "输入格式错误"
    left_part, right_part = parts[0], parts[1]
    return f"{right_part}/{left_part}"


def prepare_variant_msg(file1, file2):
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

def get_gene(df:pd.DataFrame, report_json_file, result_tsv_file, outcall_file):
    
    
      
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
        if drugs['name'].__contains__('/'):
            for drug in drugs['name'].split('/'):
                new_row=[drug,drugs['id'], gene,diplotype,annotation_type]
                df.loc[len(df)] = new_row  
        else:
            new_row=[drugs['name'], drugs['id'], gene,diplotype,annotation_type]
            df.loc[len(df)] = new_row  
    return df

def find_guideline(guidelines, df, i, drug_name, guideline, gene, diplotype):
    if drug_name == 'allopurinol' and gene=='HLA-B':
        print('allopurinol')
    if len(guidelines) == 0:
        raise ValueError(f"No guidelines found for drug {drug_name} in {guideline} for gene {gene} with diplotype {diplotype}") 
    for guideline_item in guidelines:
        annotations = guideline_item['annotations'] if 'annotations' in guideline_item.keys() else []
        if len(annotations) == 0:
            print(f"No annotations found for drug {drug_name} in {guideline} for gene {gene} with diplotype {diplotype}")
            continue
        for annotation in annotations:
            for genotype in annotation['genotypes']:
                # diplotypes = list(filter(lambda x: x['allele1']['gene'] == gene, genotype['diplotypes']))
                diplotypes = genotype['diplotypes']
                flag = False
                is_in_df = False
                di_gene = ''
                if len(diplotypes) == 0:
                    print(f"No diplotypes found for drug {drug_name} in {guideline} for gene {gene} with diplotype {diplotype}") 
                    continue
                # 这里循环必须走完，因为存在多个基因控制一个性状，目前策略是只要一个基因二倍型对不上，这个药物指南就不适用
                for diplotype in diplotypes:
                    di_gene = diplotype['gene']
                    type = diplotype['allele1']['name']+'/'+diplotype['allele2']['name']
                    r_type = reverse_string(type)
                    #保证在diplotypes里面有一个diplotype是当前df[i]的diplotype
                    if (type == df.iloc[i,3]) and (df.iloc[i,2] == diplotype['gene']):
                        is_in_df = True
                    if (r_type == df.iloc[i,3]) and (df.iloc[i,2] == diplotype['gene']):
                        is_in_df = True
                    if len(df[(df['gene'] == di_gene)&((df['diplotype']==type)|(df['diplotype']==r_type))])!=0:
                        print(f"Warning: {drug_name}  for gene {di_gene} with diplotype {diplotype['label']} found in the input data")
                        flag = True
                    else: 
                        flag = False
                        break
                        print(f"Warning: {drug_name}  for gene {di_gene} with diplotype {diplotype['label']} not found in the input data")   
                if (flag and is_in_df):
                    # 判断diplotypes中的基因二倍型是否在df中
                    df.iloc[i, 6] = annotation['drugRecommendation'] 
                    df.iloc[i, 7] = annotation['classification']
                    if len(annotation['implications']) == 0:
                        df.iloc[i, 5] = 'Not Found'
                    elif len(annotation['implications'])==1:
                        df.iloc[i, 5] = annotation['implications'][0]
                    else:
                        implication_filter = list(filter(lambda x:x.startswith(gene),annotation['implications']))
                        #df.iloc[i, 5] = implication_filter[0]
                        df.iloc[i, 5] = '\n'.join(annotation['implications'])
                else:
                    df.iloc[i, 5] = 'Not Found'
                    df.iloc[i, 6] = 'Not Found'
                    df.iloc[i, 7] = 'Not Found'       
                    
    return df

def get_all_guideline(df:pd.DataFrame,data,annotation_type):
    _guide = data['drugs'][annotation_type]
    for i in range(len(df)):
        drug_name=df.iloc[i, 0]
        guideline = df.iloc[i, 1]
        gene = df.iloc[i, 2]
        diplotype = df.iloc[i, 3]
        ref_guide = df.iloc[i, 4]
        
        if _guide.__contains__(drug_name):
            try:
                find_guideline(_guide[drug_name]['guidelines'], df, i, drug_name, guideline, gene, diplotype)
            except Exception as e:
                df.iloc[i, 5] = 'Not Found'
                df.iloc[i, 6] = 'Not Found'
                df.iloc[i, 7] = 'Not Found' 
        else:
            print(f"Warning: {drug_name} not found in CPIC Guideline Annotation")
            df.iloc[i, 5] = 'Not Found'
            df.iloc[i, 6] = 'Not Found'
            df.iloc[i, 7] = 'Not Found'
def get_report(report_json_file, diplotype_gene_set,unknown_type_gene=None):    
    columns = ['drug', 'drug_id','gene','diplotype','ref_guide']
    df1 = pd.DataFrame(columns=columns)
    df2 = pd.DataFrame(columns=columns)
    df3 = pd.DataFrame(columns=columns)
    df4 = pd.DataFrame(columns=columns)
    with open(report_json_file, 'r',encoding='utf-8') as file:
        data = json.load(file)
        df_cpic = get_annotation(diplotype_gene_set, df1, data,DRUG_ANNOTATION['CPIC'],unknown_type_gene)
        df_cpic.to_csv('cpic.csv', index=False, encoding='utf-8')
        df_dpwg = get_annotation(diplotype_gene_set, df2, data,DRUG_ANNOTATION['DPWG'],unknown_type_gene)
        df_dpwg.to_csv('dpwg.csv', index=False, encoding='utf-8')
        df_fda = get_annotation(diplotype_gene_set, df3, data,DRUG_ANNOTATION['FDA'],unknown_type_gene)
        df_fda.to_csv('fda.csv', index=False, encoding='utf-8')
        df_fda_pgx = get_annotation(diplotype_gene_set, df4, data,DRUG_ANNOTATION['FDA-PGx'],unknown_type_gene)
        df_fda_pgx.to_csv('fda_pgx.csv', index=False, encoding='utf-8')
        df = pd.concat([df_cpic, df_dpwg, df_fda, df_fda_pgx]).reset_index(drop=True)
        df.to_csv('all_guidelines.csv', index=False, encoding='utf-8')
        return df

def get_annotation(diplotype_gene_set, df, data, annotation_type,unknown_type_gene):
    
    for gene in diplotype_gene_set:
        if gene in unknown_type_gene:
            diplotype = 'Unknown/Unknown'
            continue
        if gene in data['genes']['CPIC'].keys():
            related_drugs = data['genes']['CPIC'][gene]['relatedDrugs']
            diplotype = data['genes']['CPIC'][gene]['sourceDiplotypes'][0]['allele1']['name'] \
                            + '/' + data['genes']['CPIC'][gene]['sourceDiplotypes'][0]['allele2']['name']
            
            df = add_drugs(df,related_drugs, gene,diplotype,annotation_type)
        else:
            print(f"Warning: {gene} not found in CPIC Guideline Annotation")
            new_row = ['', '', gene, '', annotation_type]
            df.loc[len(df)] = new_row
    df['implications'] = ''
    df['recommendation'] = ''
    df['classification'] = ''
    results = get_all_guideline(df,data,annotation_type)
    print(f"results:{results}")
    return df

def get_phenotype(report_json_file, diplotype_gene_set, results_table):
    results_table['phenotype'] = ''
    with open(report_json_file, 'r',encoding='utf-8') as file:
        data = json.load(file)
        for gene in diplotype_gene_set:
            if gene in data['genes']['CPIC']:
                phenotype = ';'.join(data['genes']['CPIC'][gene]['sourceDiplotypes'][0]['phenotypes'])
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
def main():
    # describe the tool
    parser = argparse.ArgumentParser(description='Extract results from PharmCAT JSONs into a TSV-formatted file.')
    # list input arguments
    parser.add_argument("-report_json_file", type=str,help="Pharmcat report result file",default='data/E-B22444326276.report.json')
    parser.add_argument("-result_tsv_file", type=str,help="Pharmcat result tsv file",default='data/E-B22444326276.report.tsv')
    parser.add_argument("-outcall_file", type=str,help="pharcat outcall file",default='data/outcall.tsv')
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
    #condition = (
    #gene_table['Source Diplotype'].str.contains(' OR | AND ') &
    #~gene_table['Gene'].isin(OUTCALL_GENE)
#)
    #gene_table.loc[condition, 'Source Diplotype'] = 'Unknown/Unknown'
    df = get_gene(df, args.report_json_file, args.result_tsv_file, args.outcall_file)
    
    df.to_csv('variant_msg.csv', index=False, encoding='utf-8')
    diplotype_gene_set = set(df['gene'])
    print(f'diplotype_gene_set:{diplotype_gene_set}')

    drug_gene_table = get_report(args.report_json_file, diplotype_gene_set,unknown_type_gene)
    results_table = drug_gene_table.merge(df,on='gene', how='left')
    results_table.to_csv('results_table.csv', index=False, encoding='utf-8')
    cols_to_check = ['implications', 'recommendation', 'classification']
    mask = (results_table[cols_to_check] == 'Not Found').all(axis=1)
    results_table = results_table[~mask]
    results_table = get_phenotype(args.report_json_file,diplotype_gene_set,results_table)
    results_table['zyg'] = results_table['diplotype'].apply(lambda x: 'Hom' if x.split('/')[0] == x.split('/')[1] else 'Het')
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
        'implications': 'advice',
        'recommendation': 'suggest',
        # 'classification': 'suggest',
        'ref_guide': 'ref_guide',
        'drug_id': 'drugid',
        'rsid': 'rsID',
        'location': 'location'
    })
    result_df = results_table.groupby(['gene', 'drug', 'diplotype']).apply(lambda group: group[group['ref_guide'].str.startswith('CPIC')]
).reset_index(drop=True)
    #results_table = add_filter_gene(results_table,filtered_df,args.sample_id)
    results_table.to_csv(f'{args.sample_id}.pgx.tsv', index=False, encoding='utf-8',sep='\t')
    result_df.to_csv(f'{args.sample_id}.pgx.result.tsv', index=False, encoding='utf-8',sep='\t')
    # print(f"drug_gene_table:{drug_gene_table}")
    # print(drug_gene_table)
if __name__ == "__main__":
    main()