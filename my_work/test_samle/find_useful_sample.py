import os
import pandas as pd

def find_useful_samples(base_dir):
    """
    在 base_dir 目录下，根据 allsite.csv 和 bam.csv 中列出的样本，
    检查是否存在对应的 vcf、vcf.idx、bam、bam.bai 四个文件，
    将满足条件的样本整理成 5 列表格返回。
    """
    allsite_csv = os.path.join(base_dir, "allsite.csv")
    bam_csv     = os.path.join(base_dir, "bam.csv")

    # 读取两个 csv，假设第一列是样本 ID
    allsite_df = pd.read_csv(allsite_csv)
    bam_df     = pd.read_csv(bam_csv)
    # 取交集，确保样本在两个表中都出现
    # 从 allsite.csv 中提取以 .allsite.vcf 结尾的文件名，并截取样本 ID
    vcf_samples = allsite_df.loc[
        allsite_df['文件名称'].str.endswith('.allsite.vcf.gz'), '文件名称'
    ].str.replace('.allsite.vcf.gz', '', regex=False)
    # 从 allsite.csv 中提取以 .allsite.vcf.idx 结尾的文件名，并截取样本 ID
    vcf_index_samples = allsite_df.loc[
        allsite_df['文件名称'].str.endswith('.allsite.vcf.gz.tbi'), '文件名称'
    ].str.replace('.allsite.vcf.gz.tbi', '', regex=False)   
    
    
    # 从 bam.csv 中提取以 .bam 结尾的文件名，并截取样本 ID
    bam_samples = bam_df.loc[
        bam_df['文件名称'].str.endswith('.bam'), '文件名称'
    ].str.replace('.bqsr.bam', '', regex=False)
    # 从 bam.csv 中提取以 .bam.bai 结尾的文件名，并截取样本 ID
    bam_index_samples = bam_df.loc[
        bam_df['文件名称'].str.endswith('.bam.bai'), '文件名称'
    ].str.replace('.bqsr.bam.bai', '', regex=False)
    # 取交集，确保样本在两个表中都出现
    common_samples = set(vcf_samples) & set(vcf_index_samples) & set(bam_samples) & set(bam_index_samples)
    
    useful_records = []
    # 在 allsite_df 与 bam_df 中分别建立“样本ID -> 完整文件名”的映射
    

    useful_records = []
    for sample_id in common_samples:
        # 构造四个必备文件路径
        bam_path = bam_df.loc[
            bam_df['文件名称']==f"{sample_id}.bqsr.bam", '文件路径'
        ].values[0]+"/"+f"{sample_id}.bqsr.bam"
        bam_bai_path = bam_df.loc[
            bam_df['文件名称']==f"{sample_id}.bqsr.bam.bai", '文件路径'
        ].values[0]+"/"+f"{sample_id}.bqsr.bam.bai"
        vcf_path = allsite_df.loc[
            allsite_df['文件名称']==f"{sample_id}.allsite.vcf.gz", '文件路径'
        ].values[0]+"/"+f"{sample_id}.allsite.vcf.gz"
        vcf_idx_path = allsite_df.loc[
            allsite_df['文件名称']==f"{sample_id}.allsite.vcf.gz.tbi", '文件路径'
        ].values[0]+"/"+f"{sample_id}.allsite.vcf.gz.tbi"
        
        # 四个文件同时存在才保留
        useful_records.append({
            "sample_id": sample_id,
            "vcf_path": vcf_path,
            "vcf_idx_path": vcf_idx_path,
            "bam_path": bam_path,
            "bam_bai_path": bam_bai_path
        })

    # 生成结果 DataFrame 并保存
    result_df = pd.DataFrame(useful_records, columns=["sample_id", "vcf_path", "vcf_idx_path", "bam_path", "bam_bai_path"])
    output_csv = os.path.join(base_dir, "useful_samples.csv")
    result_df.to_csv(output_csv, index=False)
    print(f"可用样本表格已保存至: {output_csv}")
    return result_df

if __name__ == "__main__":
    # 默认脚本所在目录即为 base_dir
    base_dir = os.path.dirname(os.path.abspath(__file__))
    find_useful_samples(base_dir)
    
    
