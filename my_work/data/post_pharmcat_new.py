import json
import pandas as pd
import argparse
import numpy as np
from handle_results import *



def get_relations(annotation_type,data):
    gene_drug_relation = []
    CPIC_gene = data[annotation_type]
    for gene in CPIC_gene.keys():
        gene_name = gene
        gene_data = CPIC_gene[gene]
        sourceDiplotypes = gene_data['sourceDiplotypes']
        if len(sourceDiplotypes) > 1:
            continue
        genotype = sourceDiplotypes[0]['label']
        drugs = gene_data['relatedDrugs']
        for drug in drugs:
            drug_name = drug['name']
            gene_drug_relation.append({
                'gene': gene_name,
                'genotype': genotype,
                'drug': drug_name,
                'annotation_type': annotation_type,
                'phenotypes': ','.join(sourceDiplotypes[0]['phenotypes'])
            })
    return gene_drug_relation


def read_gene_data(report_file):
    
    # 读取JSON文件
    with open(report_file, 'r', encoding='utf-8') as file:
        data = json.load(file)
    data = data['genes']
    CPIC_gene_drug_relation = get_relations('CPIC',data)
    DPWG_gene_drug_relation = get_relations('DPWG',data)
    CPIC_gene_drug_relation.extend(DPWG_gene_drug_relation)
    gene_drug_relation_df = pd.DataFrame(CPIC_gene_drug_relation)
# 按gene，genotype，drug进行分组，遇到重复时，保留annotation_type为CPIC的数据
    gene_drug_relation_df = gene_drug_relation_df.sort_values(by=['annotation_type'], ascending=[False])
    gene_drug_relation_df = gene_drug_relation_df.groupby(['gene', 'genotype', 'drug']).first().reset_index()
    # gene_drug_relation_df.to_csv('data/test_gene_drug_relation.csv',index=False)
    return gene_drug_relation_df


def get_drug_guideline_relations(drugs,anno_type):
    drug_drug_relation = []
    empty_annotations = []
    anntype_drugs = drugs[anno_type]
    for drug_name in anntype_drugs.keys():
        drug_data = anntype_drugs[drug_name]
        guidelines = drug_data['guidelines']
        for guideline in guidelines:
            annotations = guideline['annotations']
            if len(annotations) == 0:
                drug_drug_relation.append({
                            'drug': drug_name,
                            'drugRecommendation': '.',
                            'genotype': '.',
                            'gene': '.',
                            'annotation_type': anno_type,
                            # 'phenotypes': ','.join(sourceDiplotypes[0]['phenotypes'])
                        })
                
                continue
            for annotation in annotations:
                drugRecommendation = annotation['drugRecommendation']
                genotypes = annotation['genotypes']
                for genotype in genotypes:
                    diplotypes = genotype['diplotypes']
                    for diplotype in diplotypes:
                        gene = diplotype['gene']
                        genotype = diplotype['label']
                        drug_drug_relation.append({
                            'drug': drug_name,
                            'drugRecommendation': drugRecommendation,
                            'genotype': genotype,
                            'gene': gene,
                            'annotation_type': anno_type,
                            # 'phenotypes': ','.join(sourceDiplotypes[0]['phenotypes'])
                        })
    return drug_drug_relation
def read_drugs_data(report_file):
    # 读取JSON文件
    with open(report_file, 'r', encoding='utf-8') as file:
        data = json.load(file)
    drugs = data['drugs']
    CPIC_drug_guideline_relation = get_drug_guideline_relations(drugs,'CPIC Guideline Annotation')
    DPWG_drug_guideline_relation = get_drug_guideline_relations(drugs,'DPWG Guideline Annotation')
    FDA_LABEL_drug_guideline_relation = get_drug_guideline_relations(drugs,'FDA Label Annotation')
    FDA_PGX_drug_guideline_relation = get_drug_guideline_relations(drugs,'FDA PGx Association')
    
    all_drug_guidelines = CPIC_drug_guideline_relation+(DPWG_drug_guideline_relation)+(FDA_LABEL_drug_guideline_relation)+(FDA_PGX_drug_guideline_relation)
    drug_guideline_relation_df = pd.DataFrame(all_drug_guidelines)
    # drug_guideline_relation_df.to_csv('data/test_drug_guideline_relation.csv',index=False)
    return drug_guideline_relation_df                    



def merge_gene_drug_guideline(reprot_file):
    # 调用函数处理数据
    gene_drug_relation_df = read_gene_data(reprot_file)
    drug_guideline_relation_df = read_drugs_data(reprot_file)
    merged_df = pd.merge(gene_drug_relation_df, drug_guideline_relation_df, on=['drug'], how='left')
    # 提取gene_y!='.'并且gene_x!=gene_y的数据
    condition1 = merged_df['gene_y'] != '.'
    condition2 = merged_df['gene_x'] == merged_df['gene_y']
    part1 = merged_df[condition1 & condition2]
    
    # 提取gene_y=='.'的数据
    part2 = merged_df[merged_df['gene_y'] == '.']
    
    # 合并两个数据
    merged_df = pd.concat([part1, part2])
    
    # 筛选出 gene_x 和 drug 组合下，删除 gene_y 为 '.' 且存在非 '.' 的数据
    # 先按 gene_x 和 drug 分组，判断每组中 gene_y 是否有非 '.' 的值
    group_condition = merged_df.groupby(['gene_x', 'drug'])['gene_y'].transform(lambda x: (x != '.').any())
    # 筛选出不需要删除的数据：要么 gene_y 不为 '.'，要么该组中没有非 '.' 的 gene_y 值
    merged_df = merged_df[~((merged_df['gene_y'] == '.') & group_condition)]
    # merged_df.to_csv('data/test_1.csv',index=False)
    # 定义 annotation_type_y 的优先级
    priority_map = {
        'CPIC Guideline Annotation': 0,
        'DPWG Guideline Annotation': 1,
        'FDA Label Annotation': 2,
        'FDA PGx Association': 3
    }
    # 添加优先级列
    merged_df['priority'] = merged_df['annotation_type_y'].map(priority_map)
    # 按优先级排序
    merged_df = merged_df.sort_values(by=['priority'])
    # 按 gene_x 和 drug 分组，每组保留第一条数据
    merged_df = merged_df.groupby(['gene_x', 'drug']).first().reset_index()
    # 删除优先级列
    merged_df = merged_df.drop(columns=['priority'])
    # merged_df.to_csv('data/test_merged_gene_drug_guideline.csv',index=False)
    return merged_df



def prepare_variant_msg(intersect_file, prepharmcat_file):
    print(f'file1, file2:{intersect_file}, {prepharmcat_file}')
    df = pd.read_csv(intersect_file,sep='\t',dtype=str)
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
    prepharmcat_data = pd.read_csv(prepharmcat_file, sep='\t',dtype=str)
    prepharmcat_data['PX'] = prepharmcat_data['INFO'].str.extract(r'PX=([^;]+)')
    prepharmcat_data['cHGVS'] = 'test'
    prepharmcat_data['pHGVS'] = 'test'
    # prepharmcat_data['cHGVS'] = prepharmcat_data['INFO'].str.split('|').str[9]
    # prepharmcat_data['pHGVS'] = prepharmcat_data['INFO'].str.split('|').str[10]
    df['gene'] = prepharmcat_data['PX']
    df['rsid'] = prepharmcat_data['ID']
    df['cHGVS'] = prepharmcat_data['cHGVS']
    df['pHGVS'] = prepharmcat_data['pHGVS']
    grouped_df = df.groupby('gene', as_index=False).agg({
        'cHGVS': lambda x: ','.join(set(x)),  
        'pHGVS': lambda x: ','.join(set(x)),  
        'rsid': lambda x: ','.join(set(x)),  
        'location': lambda x: ','.join(set(x)),  
        
        'QUAL': lambda x: ','.join(set(x)),
        'depth': lambda x: ';'.join(set(x)),
        'readsratio': lambda x: ';'.join(set(x)),
 
    })
    print(f"grouped_df:{grouped_df}")
    # grouped_df.to_csv('grouped_df.csv',index=False)
    return grouped_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract results from PharmCAT JSONs into a TSV-formatted file.')
    # list input arguments
    parser.add_argument("-report_json_file", type=str,help="Pharmcat report result file",default='data/test_new.report.json')
    parser.add_argument("-result_tsv_file", type=str,help="Pharmcat result tsv file",default='data/test_new.report.tsv')
    parser.add_argument("-outcall_file", type=str,help="pharcat outcall file",default='data/test_new_outcall.tsv')
    parser.add_argument("-intersect_file", type=str,help="pharmcat pipeline inputfile",default='data/intersect.csv')
    parser.add_argument("-prepharmcat_file", type=str,help="",default='data/annotated.csv')
    parser.add_argument("-sample_id", type=str,help="",default='Test')
    args = parser.parse_args()
    print(f"args:{args}")
    var_df = prepare_variant_msg(args.intersect_file, args.prepharmcat_file)
    guideline_df = merge_gene_drug_guideline(args.report_json_file)
    guideline_df.rename(columns={'gene_x': 'gene'}, inplace=True)
    merged_df = pd.merge(var_df, guideline_df, on=['gene'], how='right')
    merged_df['drugid'] = ''
    merged_df = merged_df[['gene','drug','genotype_x','rsid','cHGVS','pHGVS','QUAL','depth','readsratio','phenotypes','drugRecommendation','annotation_type_y','location','drugid']]
    # 创建mask筛选出包含斜杠的行
    merged_df.rename(columns={'genotype_x': 'diplotype'
                              
                              }, inplace=True)
    mask = merged_df['diplotype'].str.contains('/')

    # 对符合条件的行应用lambda函数计算zyg列
    merged_df.loc[mask, 'zyg'] = merged_df.loc[mask, 'diplotype'].apply(
        lambda x: 'Hom' if x.split('/')[0] == x.split('/')[1] else 'Het'
    )
    # results_table['zyg'] = results_table['diplotype'].apply(lambda x: 'Hom' if x.split('/')[0] == x.split('/')[1] else 'Het')
    merged_df['sample']= args.sample_id  
    merged_df = merged_df.rename(columns={
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
        'phenotypes':'effect',
        'guide': 'guide',
        'implications': 'advice',
        'drugRecommendation': 'suggest',
        # 'classification': 'suggest',
        'ref_guide': 'ref_guide',
        'drug_id': 'drugid',
        'rsid': 'rsID',
        'location': 'location'
    })
    known_type_gene_table = get_short_guide(merged_df)
    known_type_gene_table.to_csv(f'{args.sample_id}_drug.tsv',index=False)
    