import json
import csv
import os
import pandas as pd

# 输入文件路径
input_file = 'reporter/prescribing_guidance.json'
output_file = 'reporter/drug_genotype_variant.xlsx'

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
    guideline_entity = guideline.get('guideline', {})
    related_chemicals = guideline_entity.get('relatedChemicals', [])
    related_genes = guideline_entity.get('relatedGenes', [])
    related_alleles = guideline_entity.get('relatedAlleles', [])
    for chemical in related_chemicals:
        chemical_name = chemical.get('name', '')
        chemical_id = chemical.get('id', '')
        for gene in related_genes:
            gene_name = gene.get('symbol', '')
            gene_id = gene.get('id', '')
            if len(related_alleles) == 0:
                records.append({
                    'Guideline Name': guideline_name,
                    'Gene': gene_name,
                    'Gene ID': gene_id,
                    'Allele': '',
                    'Allele ID': '',
                    'Chemical': chemical_name,
                    'Chemical ID': chemical_id,
                    'Allele Symbol': ''
                })
            else:
                for allele in related_alleles:
                    allele_name = allele.get('name', '')
                    allele_id = allele.get('id', '')
                    allele_symbol = allele.get('symbol', '')
                    records.append({
                        'Guideline Name': guideline_name,
                        'Gene': gene_name,
                        'Gene ID': gene_id,
                        'Allele': allele_name,
                        'Allele ID': allele_id,
                        'Chemical': chemical_name,
                        'Chemical ID': chemical_id,
                        'Allele Symbol': allele_symbol
                    })

# 生成 DataFrame 并写出 CSV
df = pd.DataFrame(records, columns=[
    'Guideline Name',
    'Gene',
    'Gene ID',
    'Allele',
    'Allele ID',
    'Chemical',
    'Chemical ID',
    'Allele Symbol'
])

df.to_excel(output_file, index=False, encoding='utf-8')
print(f"成功提取 {len(df)} 条 recommendation 信息，已保存至 {output_file}")