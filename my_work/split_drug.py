import pandas as pd
import csv
# 1. 读入数据
with open('drug_results.txt', 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    rows = list(reader)          # 这里已经把每行正确地读成 dict

# 2. 转成 DataFrame
df = pd.DataFrame(rows)

def split_drug(df):
    # 2. 拆分 drug 列中含 "/" 的行
    has_slash = df['drug'].str.contains('/', na=False)
    df_to_split = df[has_slash].copy()
    df_intact   = df[~has_slash].copy()

    # 3. 真正拆分
    #    先按 / 拆成多列，再 stack 成 Series，再还原成 DataFrame
    split_series = (
        df_to_split['drug']
        .str.split('/', expand=True)   # 拆成多列
        .stack()                       # 变成长 Series
        .reset_index(level=1, drop=True)
        .rename('drug')                # 保持列名一致
    )

    # 4. 把拆分后的 drug 拼回原 DataFrame
    df_split = (
        df_to_split
        .drop(columns=['drug'])
        .join(split_series)
        .reset_index(drop=True)
    )

    # 5. 合并结果
    result = pd.concat([df_intact, df_split], ignore_index=True)

    # 6. 若需要，可重新排序列
    result = result[df.columns]

# 7. 查看结果
print(result)