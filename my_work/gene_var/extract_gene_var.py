import pandas as pd
import re
from pathlib import Path

folder = Path("gene_var")
records = []

for file_path in folder.glob("*_allele_definition_table.xlsx"):
    if file_path.name in {"MT-RNR1_allele_definition_table.xlsx"}:
        print(file_path.name, "跳过")
    gene = file_path.stem.replace("_allele_definition_table", "")
    df = pd.read_excel(file_path, header=None)

    # 1. 染色体号（第1列第4行）
    chrom_line = str(df.iloc[3, 0]).strip()
    chrom_match = re.search(r'chromosome\s+([0-9]{1,2}|X|Y|MT|M)', chrom_line, flags=re.I)
    if chrom_match:
        
        chrom = chrom_match.group(1).upper()
    else:
        chrom = 'chrZ'

    # 2. rsID（第6行，从第2列开始）
    rs_row = df.iloc[5, 1:]      # 跳过第一列

    # 3. 基因组位置行（第4行，从第2列开始）
    pos_row = df.iloc[3, 1:]     # 跳过第一列

    # 4. 等位基因名称列（第9行开始，第0列）
    allele_col = df.iloc[8:, 0]
    # 5. 蛋白质变化行（第3行，从第2列开始）
    protein_change_row = df.iloc[2, 1:]  # 从第2列开始
    
    #6. 基因变化行（第四行，从第二列开始）
    gene_change_row = df.iloc[3, 1:]  # 从第2列开始
    # 5. 遍历每一列（从第2列开始）
    for col_idx in range(1, df.shape[1]):
        rs_id = str(rs_row.iloc[col_idx - 1]).strip()
        pos_str = str(pos_row.iloc[col_idx - 1]).strip()
        protein_change = str(protein_change_row.iloc[col_idx - 1]).strip()
        gene_change = str(gene_change_row.iloc[col_idx - 1]).strip()
        if pd.isna(rs_id) or pd.isna(pos_str):
            continue

        # 解析 g.12345A>G
        
        m = re.match(r'g\.(\d+)([ACGT])>([ACGT])', pos_str)
        if not m:
            if not pos_str.startswith("g."):
                variant = f"chr{chrom}:{pos_str}"
            else:
                idx = pos_str.index("g.")
                pos_str = pos_str[idx + 2:]
                variant = f"chr{chrom}:{pos_str}"
        else:
            pos, ref, alt = m.groups()
            variant = f"chr{chrom}:{pos}:{ref}>{alt}"

        
        # 6. 遍历该列从第9行开始
        for row_offset in range(8, len(df)):
            allele = str(df.iloc[row_offset, 0]).strip()
            base = str(df.iloc[row_offset, col_idx]).strip()
            if not (base.__contains__("A") or base.__contains__("C") or base.__contains__("G") or base.__contains__("T")):
                continue
            records.append({
                "Gene": gene,
                "rsID": rs_id,
                "Variant": variant,
                "Ref": "hg38",
                "Haplotype": f"{gene}{allele}",
                "ProteinChange": protein_change,
                "GeneChange": gene_change
            })
            #print(f"Processed {gene} - {rs_id} - {variant} - {allele}:{base}")
    print(f"Finished processing {gene}")
# 7. 写出 Excel
out_df = pd.DataFrame(records, columns=["Gene", "rsID", "Variant", "Ref", "Haplotype", "ProteinChange", "GeneChange"])
out_df.to_excel("gene_var/all_variants.xlsx", index=False)
print("结果已保存为 gene_var/all_variants.xlsx")