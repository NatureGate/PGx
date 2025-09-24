import json
import csv
import os
import pandas as pd

# 输入文件路径
input_file = 'reporter/prescribing_guidance.json'
output_file = 'reporter/recommendations_summary.xlsx'

# 检查输入文件是否存在
if not os.path.isfile(input_file):
    raise FileNotFoundError(f"输入文件 {input_file} 未找到，请确保文件路径正确。")

# 读取 JSON 文件
with open(input_file, 'r', encoding='utf-8') as f:
    data = json.load(f)

# 准备输出数据
records = []

# 遍历所有 guideline
for guideline in data.get('guidelines', []):
    guideline_name = guideline.get('guideline', {}).get('name', 'Unknown Guideline')
    related_genes = ', '.join([g.get('symbol', '') for g in guideline.get('guideline', {}).get('relatedGenes', [])])
    # 遍历该 guideline 下的所有 recommendations
    for rec in guideline.get('recommendations', []):
        rec_id = rec.get('id', '')
        name = rec.get('name', '')
        population = rec.get('population', '')
        # classification = rec.get('classification', {}).get('term', '')
        
        # 获取相关药物列表
        related_chemicals = ', '.join([c.get('name', '') for c in rec.get('relatedChemicals', [])])
        
        # 获取文本内容（去除 HTML 标签）
        text_html = rec.get('text', {}).get('html', '')
        text_content = text_html.strip()
        if text_html.startswith('<p>xN is not a valid input for CYP2D6'):
            text_content = '.'
        
        # 获取基因型含义
        implications = '; '.join(rec.get('implications', []))
        
        # 获取 lookupKey
        lookup_keys = '; '.join([f"{k}={v}" for k, v in rec.get('lookupKey', {}).items()])
        
        # 添加到输出数据
        records.append({
            'related_genes': related_genes,
            'Guideline Name': guideline_name,
            #'Recommendation ID': rec_id,
            'Name': name,
            #'Population': population,
            #'Classification': classification,
            'Related Chemicals': related_chemicals,
            'Text': text_content,
            'Implications': implications,
            'Lookup Key': lookup_keys
        })

# 生成 DataFrame 并写出 CSV
df = pd.DataFrame(records, columns=[
    'related_genes',
    'Guideline Name',
    'Name',
    'Related Chemicals',
    'Text',
    'Implications',
    'Lookup Key'
])

df.to_excel(output_file, index=False, encoding='utf-8')
print(f"成功提取 {len(df)} 条 recommendation 信息，已保存至 {output_file}")