#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从phenotyper流程结果中筛选：
sourceDiplotypes 长度 != 1 的基因；

"""

from typing import Any


import json
from pathlib import Path
import argparse

# --------------- 配置 ---------------
# INPUT_FILE  = 'problem_handle/1900617694C.preprocessed.new.phenotype.json'   # 原始文件
                     # 结果保存路径
# --------------- 配置 ---------------

def main(input_file, outcall_file):
    data = json.loads(Path(input_file).read_text(encoding='utf-8'))


    for gene_symbol, gene_data in data['geneReports']['CPIC'].items():
        multi_phenotype = []  
        zero_phenotype_gene = set()
        dts = gene_data.get('sourceDiplotypes', [])
        gene = gene_data.get('geneSymbol', '')
        variants = gene_data.get('variants', '')
        # sourceDiplotypes 长度 ≠ 1
        if len(dts) > 1:
            phenotype_all = []
            for diplotype in dts:
                phenotypes = diplotype.get('phenotypes', [])
                genotype = diplotype.get('label', '')
                if len(phenotypes) > 0:
                    for phenotype in phenotypes:
                        if phenotype in phenotype_all:
                            continue
                        multi_phenotype.append({                    
                            'gene': gene,
                            'genotype': genotype,
                            'phenotype': phenotype,
                        })
                        phenotype_all.append(phenotype)
            if len(phenotype_all) == 1:
                # 写出gene\tgenotype\tphenotype
                # 将 gene、genotype 追加写入新文件（若文件不存在则自动创建）
                with open(outcall_file, 'a', encoding='utf-8') as f:
                    f.write(f'{gene}\t{multi_phenotype[0]["genotype"]}\n')
                print(f'{gene}\t{multi_phenotype[0]["genotype"]}\t{multi_phenotype[0]["phenotype"]}')
                continue

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='输入pharmcat的phenotyper结果JSON文件路径')
    parser.add_argument('--phenotype_result', required=True, help='输入的phenotype结果JSON文件路径')
    parser.add_argument('--outcall_file', required=True, help='输出的outcall文件路径')
    args = parser.parse_args()
    input_file = args.phenotype_result
    outcall_file = args.outcall_file
    main(input_file, outcall_file)
