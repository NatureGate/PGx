import pandas as pd
import argparse


GENE = ['ABCG2','CACNA1S','CYP2B6','CYP2C19','CYP2C9','CYP2D6','CYP3A5','CYP4F2','CFTR','DPYD','G6PD','HLA-A','HLA-B','IFNL3','MT-RNR1','NUDT15','RYR1','SLCO1B1','TPMT','UGT1A1','VKORC1']
def load_drug_data(file_path):
    df = pd.read_table(file_path,sep='\t')
    return df


# 从merged_vcf提取rsID
def extract_rsid(ids_file):
    with open(ids_file, 'r') as f:
        rs_ids = f.readlines()
    results_ids = [x.strip() for x in rs_ids if x.strip()]
    return results_ids

def get_var_msg(df:pd.DataFrame,rs_ids,gene):
    filter_df = df[(~df['Gene'].isin(gene)) & (df['Variant/Haplotypes'].isin(rs_ids))]
    return filter_df

if __name__ == "__main__":
    # Example usage
    parser = argparse.ArgumentParser()
    # Add arguments to the parser
    # 这里可以添加更多的参数，例如输入文件路径等
    parser.add_argument('--valid_ids_path', type=str, help='Path to the valid IDs file', default='valid_ids.txt')
    parser.add_argument('--sample_id', type=str, help='Sample ID', default='Test')

    args = parser.parse_args()

    valid_ids_path = args.valid_ids_path
    sample_id = args.sample_id
    # Load drug data 这个是PharmGKB的变异数据
    drug_data = load_drug_data('/home/stereonote/data/var_drug_ann.tsv')
    print(drug_data.head())
    
    # Extract rsID from VCF file
    rs_ids = extract_rsid(valid_ids_path)
    
    # Example of getting variant message
    gene = GENE
    results = get_var_msg(drug_data, rs_ids, gene)
    results.to_csv(f'{sample_id}_other.pgx.tsv', sep='\t', index=False)
    print(f"Filtered results for genes {gene}:")
    print(results.head())
    print(f"Number of filtered results: {results.shape[0]}")

