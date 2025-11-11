import pandas as pd
from pathlib import Path

data_dir = Path(r'./pgx_drugs_result/')          # 改成你的目录
files  = data_dir.glob('*.tsv')     # 找到所有 tsv

# 一次性读完再合并
df = pd.concat(
        [pd.read_table(f, sep=',') for f in files],
        ignore_index=True           # 重新编号行索引
     )
print(df.shape)

mask = df['drug_id'].isna()