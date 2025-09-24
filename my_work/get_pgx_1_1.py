import json
import pandas as pd
from glob import glob
import os
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
        first_diplotype = diplotype_list[0]
        diplotype_val = first_diplotype.get('diplotype')
        lookupkey_val = first_diplotype.get('lookupkey')
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
df = pd.DataFrame(data)

# 显示结果（可选）
print(df)