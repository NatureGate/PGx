import json
import csv
import os
import pandas as pd

# 输入文件路径
input_file = 'reporter/prescribing_guidance.json'
# output_file = 'reporter/recommendations_summary.csv'

# 检查输入文件是否存在
if not os.path.isfile(input_file):
    raise FileNotFoundError(f"输入文件 {input_file} 未找到，请确保文件路径正确。")

# 读取 JSON 文件
with open(input_file, 'r', encoding='utf-8') as f:
    data = json.load(f)

# 准备输出数据
output_data = []
all_chemicals = set()
single_chemical_set = set()
# 遍历所有 guideline
for guideline in data.get('guidelines', []):
    guideline_name = guideline.get('guideline', {}).get('name', 'Unknown Guideline')
    related_chemicals = set([c.get('name', '') for c in guideline['guideline'].get('relatedChemicals', [])])
    
    all_chemicals.update(related_chemicals)

    # 遍历该 guideline 下的所有 recommendations
    for rec in guideline.get('recommendations', []):

        # 获取相关药物列表
        related_chemicals = set([c.get('name', '') for c in rec.get('relatedChemicals', [])])
        
        all_chemicals.update(related_chemicals)
print(len(all_chemicals))
with open('reporter/all_drugs.txt', 'w', encoding='utf-8') as f: 
    f.write('drug\n')
    for chem in all_chemicals:
        f.write(chem + '\n')
# for chem in all_chemicals:
#     if ' / ' not in chem and ' and ' not in chem :
#         single_chemical_set.add(chem)
#     elif ' and ' in chem:
#         components = chem.split(' and ')
#         single_chemical_set.update(components)
#     elif ' / ' in chem:
#         components = chem.split(' / ')
#         single_chemical_set.update(components)

# print(len(single_chemical_set))



# 将drugEN变为小写
drug_detail = pd.read_table('reporter/drug_detail.txt')
old_drug = drug_detail['drugEN'].str.lower().tolist()
print(single_chemical_set - set(old_drug))
add_drug = single_chemical_set - set(old_drug)
# 将add_drug 按字母顺序排序
add_drug = sorted(add_drug)
# 将add_drug 写出文件，一行一个，守字母大写
with open('reporter/additional_drug.txt', 'w', encoding='utf-8') as f:
    f.write('drug_en\n')
    for drug in add_drug:
        f.write(drug.capitalize() + '\n')

add_drug = pd.read_table('reporter/additional_drug.txt')
new_drug = pd.read_excel('reporter/新增药物范围.xlsx')
pd.merge(add_drug, new_drug, left_on='drug_en', right_on='drug_en', how='outer').to_excel('reporter/additional_drugs.xlsx', index=False, encoding='utf-8')


# old_drug.extend(new_drug)
# for single_chemical in single_chemical_set:
#     if single_chemical not in old_drug:
#         # print(single_chemical)
#         pass
        # old_drug.append(single_chemical)
# print(all_chemicals)

# import pandas as pd 
# all_drugs = pd.read_table('drugs.tsv')
# mask = all_drugs['PharmGKB Accession Id'].isin(all_chemicals)
# all_drugs = all_drugs[mask]
# all_drugs = all_drugs[['PharmGKB Accession Id', 'Name', 'ATC Identifiers']]
# all_drugs = all_drugs.rename(columns={'PharmGKB Accession Id': 'id', 'Name': 'name', 'ATC Identifiers': 'atc_code'})
# all_drugs.to_csv('reporter/drugs.csv', index=False, encoding='utf-8')

# new_drugs = pd.read_csv('reporter/drugs.csv')
# old_drugs = pd.read_table('reporter/drug_detail.txt')
# merged_drugs = pd.merge(new_drugs, old_drugs, on='atc_code', how='outer')
# merged_drugs.to_csv('reporter/merged_drugs.csv', index=False, encoding='utf-8')
# print(merged_drugs.head())