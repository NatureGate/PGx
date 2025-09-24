import pandas as pd
import numpy as np
import json
def prepare_gene_msg(file1, file2):
    df = pd.read_csv(file1,sep='\t')
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

# 添加 location 列
    df['location'] = (
        df['CHROM'] + ':' + 
        df['POS'].astype(str) + ':' + 
        df['REF'] + '>' + 
        df['ALT']
    )

    print(df.head())

    prepharmcat_data = pd.read_csv(file2, sep='\t')
    prepharmcat_data['PX'] = prepharmcat_data['INFO'].str.extract(r'PX=([^;]+)')
    df['gene'] = prepharmcat_data['PX']
    df['rsid'] = prepharmcat_data['ID']

    print(set(df['gene']))
    return df


