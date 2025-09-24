import pandas as pd
import re
from pathlib import Path
import warnings

# 方法 1：完全忽略所有 UserWarning
warnings.filterwarnings("ignore", category=UserWarning)

# 方法 2：仅忽略特定消息的警告（可选）
# warnings.filterwarnings("ignore", message=".*some.*text.*", category=UserWarning)
folder = Path(".")
records = []

def abbreviate_phenotype(s: pd.Series) -> pd.Series:
    """将 'Intermediate Metabolizer' 转为 'IM'"""
    return (
        s.str.split()
        .apply(lambda words: ''.join(word[0].upper() for word in words))
    )

for file_path in folder.glob("diplotype_phenotype/*_Diplotype_Phenotype_Table.xlsx"):
    gene = file_path.stem.replace("_Diplotype_Phenotype_Table", "")
    df = pd.read_excel(file_path)
    df.columns = ['diplotype', 'activity score', 'phenotype', 'phenotype abbreviation']
    gene_ = gene+' '
    df['phenotype'] = df['phenotype'].str.replace(gene, "")
    df['phenotype abbreviation'] = df['phenotype'].str.replace(gene_, "")
    df['phenotype abbreviation'] = abbreviate_phenotype(df['phenotype'])
    df['gene'] = gene
    records.append(df)
    print(f"Finished processing {gene}")
# 7. 写出 Excel
pd.concat(records, ignore_index=True).to_excel("diplotype_phenotype/diplotype_phenotype.xlsx", index=False)
print("结果已保存为 diplotype_phenotype/diplotype_phenotype.xlsx")