#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 prescribing_guidance.json 中提取所有不重复的药物名称
"""

import json
from pathlib import Path

def main():
    # 1. 读入原始文件
    src = Path('prescribing_guidance.json')
    with src.open(encoding='utf-8') as f:
        data = json.load(f)

    # 2. 用 set 去重
    drug_names = set()

    # 3. 遍历所有 guideline
    for g in data.get('guidelines', []):
        # guideline 层的药物
        for chem in g.get('guideline', {}).get('relatedChemicals', []):
            drug_names.add(chem.get('name'))

        # recommendations 层的药物
        for rec in g.get('recommendations', []):
            for chem in rec.get('relatedChemicals', []):
                drug_names.add(chem.get('name'))

    # 4. 排序输出
    drug_list = sorted(drug_names)

    # 5. 保存到文件
    dst = Path('extracted_drug_names.json')
    with dst.open('w', encoding='utf-8') as f_out:
        json.dump(drug_list, f_out, ensure_ascii=False, indent=2)

    print(f'提取完成，共 {len(drug_list)} 种不重复药物：{drug_list}')

if __name__ == '__main__':
    main()