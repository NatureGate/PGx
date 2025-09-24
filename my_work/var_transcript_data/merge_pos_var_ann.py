import pandas as pd 
# df = pd.read_csv('var_transcript_data/filter_combine_all.csv', dtype=str,sep=',')
# df2 = pd.read_table('var_transcript_data/new_position.txt', dtype=str)
# merged_df = pd.merge(df2, df, on=['CHROM', 'POS', 'ID'], how='left')
# na_data = merged_df[merged_df['Name'].isna()]
# na_data.to_csv('var_transcript_data/na_anno_data.csv', index=False)
{
    "chr4":"NC_000004.12",
    "chr1":"NC_000001.11",
    "chr7":"NC_000007.14",
    "chr19":"NC_000019.10",
    "chr10":"NC_000010.11",
    "chr22":"NC_000022.11",
    "chr7":"NC_000007.14",
    "chr19":"NC_000019.10",
    "chr1":"NC_000001.11",
    "chrX":"NC_000023.11",
    "chr19":"NC_000019.10",
    "chr8":"NC_000008.11",
    "chr13":"NC_000013.11",
    "chr19":"NC_000019.10",
    "chr12":'NC_000012.12',
    "chr6":'NC_000006.12',
    "chr2":"NC_000002.12",
    "chr16":"NC_000016.10"
}
all_data = pd.read_csv('var_transcript_data/filter_combine_all.csv', dtype=str)
na_data = pd.read_csv('var_transcript_data/na_anno_data.csv', dtype=str)[all_data.columns]
pd.concat([all_data[~all_data['Name'].isna()], na_data]).to_csv('var_transcript_data/merge_pos_var_ann.csv', index=False)   