import json
# f_w = open('gene_drug.csv', 'w', encoding='utf-8')
# f_w.write('drug\tgene\n')
# with open('guideline.json', 'r', encoding='utf-8') as f:
#     gene_drugs = json.load(f)
    
        
#         # 遍历每个基因药物指南
#     for guide_line in gene_drugs['guidelines']:
#         #print(guide_line['guideline']['relatedChemicals'])
#         relatedChemicals = guide_line.get('guideline', ).get('relatedChemicals', [])
#         realrelatedChemicalNames = ','.join([chemical['name'] for chemical in relatedChemicals])
#         relatedGenes = guide_line.get('guideline', '').get('relatedGenes', [])
#         relatedGeneNames = ','.join([gene['symbol'] for gene in relatedGenes])
#         for relatedChemical in relatedChemicals:
#             name = relatedChemical['name']
#             f_w.write(f'{name}\t{relatedGeneNames}\n')
            
import pandas as pd
df = pd.read_csv('gene_drug.csv', sep='\t')
df = df.drop_duplicates()
df.to_csv('gene_drug_relation.csv', index=False, encoding='utf-8')