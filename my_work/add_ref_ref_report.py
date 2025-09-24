import os
import json
import pandas as pd
from glob import glob

# 定义文件夹路径
folder_path = 'phenotype/cpic'

# 获取文件夹下所有JSON文件
json_files = glob(os.path.join(folder_path, '*.json'))

# 创建数据列表
data = []

for file in json_files:
    # 读取JSON文件内容
    with open(file, 'r', encoding='utf-8') as f:
        content = json.load(f)
    
    # 提取gene信息
    gene = content.get('gene')
    
    # 获取diplotype数组的第一个元素
    diplotype_list = content.get('diplotypes', [])
    if diplotype_list:
        default_diplotype = diplotype_list[0]
        if gene == 'G6PD':
            default_diplotype = diplotype_list[20]
        if gene == 'IFNL3' or gene == 'CFTR':
            continue
        
        diplotype_val = default_diplotype.get('diplotype')
        lookupkey_val = default_diplotype.get('lookupkey')
    else:
        diplotype_val = None
        lookupkey_val = None
    
    # 将数据添加到列表
    data.append({
        'gene': gene,
        'diplotype': diplotype_val,
        'lookupkey': lookupkey_val
    })

# 转换为DataFrame
phe_df = pd.DataFrame(data)


# 显示结果（可选）
print(phe_df)

# 保存为CSV文件（可选）
# df.to_csv('output.csv', index=False)
import json
import pandas as pd

# 输入JSON文件路径
input_json_file = 'reporter/prescribing_guidance.json'

# 读取JSON文件
with open(input_json_file, 'r', encoding='utf-8') as json_file:
    data = json.load(json_file)

# 初始化一个空的列表，用于存储提取的数据
rows = []

# 遍历JSON数据中的每个条目
for guideline in data['guidelines']:
    ref_guide = guideline['guideline'].get('source', '')
    for recommendation in guideline['recommendations']:
        relatedChemicals = recommendation.get('relatedChemicals', '')
        for chemical in relatedChemicals:
            drug = chemical.get('name', '')
            drug_id = chemical.get('id', '')
            html_content = recommendation['text']['html']
            if html_content.__contains__('xN is not a valid input for CYP2D6 copy number'):
                continue
            implications = '; '.join(recommendation['implications'])
            
            lookup_key = recommendation['lookupKey']
            # 提取基因名称
            genes = []
            values = []
            for key in lookup_key.keys():
                genes.append(key)
                value = lookup_key[key]
                if not isinstance(value, str):
                    value = json.dumps(value)
                values.append(value)
            # print(genes, values)
            lookupkeys = ','.join(values)
            
            # 如果有多个基因，用逗号分隔
            gene_list = ','.join(genes)
            ref_guide = ref_guide
            # 将提取的数据添加到列表中
            rows.append([drug, drug_id, html_content, implications, lookupkeys, gene_list, ref_guide])

# 创建DataFrame
df = pd.DataFrame(rows, columns=['drug', 'drugid', 'suggest', 'advice', 'lookupkey', 'gene','ref_guide'])
df = df[['drug','drugid','gene', 'suggest', 'advice', 'lookupkey','ref_guide']]
df['sample'] = '.'  # 添加一个示例样本ID列
df['gatkscore'] = '.'  # 添加一个QUAL列
df['depth'] = '.'  # 添加一个depth列
df['readsratio'] = '.'  # 添加一个readsratio列  
df['cHGVS'] = '.'  # 添加一个cHGVS列
df['pHGVS'] = '.'  # 添加一个pHGVS列
df['zyg'] = '.'  # 添加一个zyg列
df['rsID'] = '.'  # 添加一个rsid列
df['location'] = '.'  # 添加一个location列
df['effect'] = '.' 
df['guide'] = df['lookupkey']  # 添加一个guide列
# 显示DataFrame
#print(df)

# 如果需要，可以将DataFrame保存为CSV文件
#phe_df.to_csv('output.csv', sep='\t', index=False, encoding='utf-8')
merge_df = pd.merge(phe_df, df, on=['gene','lookupkey'], how='left')
merge_df.to_csv('merged_output.csv', sep='\t', index=False, encoding='utf-8')