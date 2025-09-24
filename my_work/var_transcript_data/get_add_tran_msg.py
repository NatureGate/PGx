import pandas as pd
from get_transcript_by_worm import get_transcript_by_worm
import json
transcripts_genome = ['NM_004827','NM_000069','NM_000492','NM_000767','NM_000771','NM_000769','NM_000106','NM_017460','NM_001360016','NM_000110','NM_001082','NM_000777','NM_172139','NM_018283','NM_000540','NM_006446','NM_000367','NM_000463','NM_024006']

# data = pd.read_csv('var_transcript_data/na_transcript.csv')
# mask = data['Name'].isna()
# na_transcript = data[mask]  

# def get_trans_details(na_transcript:pd.DataFrame):
#     for i in range(na_transcript.shape[0]):
#         try:
            
#             rs_id = na_transcript.iloc[i]['ID']
#             gene = na_transcript.iloc[i,7].split('=')[1]
#             ALT = na_transcript.iloc[i]['ALT']
#             REF = na_transcript.iloc[i]['REF']
#             if str(na_transcript.iloc[i]['POS'])=='18132163':
#                 print(na_transcript.iloc[i])
#             pos = na_transcript.iloc[i]['CHROM'][3:] + ':' + str(na_transcript.iloc[i]['POS'])
#             change = na_transcript.iloc[i]['REF'] + '>' + na_transcript.iloc[i]['ALT']
#             query = f"{rs_id}" if rs_id !='.' else f"?term={pos}"
#             print(f'没有找到转录本信息，尝试从NCBI爬取: {rs_id} {pos}')
#         # 4. 这里可以调用爬虫脚本
#             df = get_transcript_by_worm(query)
#             if df is not None:
#                 if ALT in df.columns:
#                     transcript = df.loc[df[ALT].str.split('.').str[0].isin(transcripts_genome), ALT].iat[0]
#                     trans_nm = transcript.split(':')[0]
#                     trans_change = transcript.split(':')[1]
#                     print('找到转录本信息:', transcript)
#                     protein_change = df.loc[df[ALT].str.startswith('NP_'), ALT].iat[0]
#                     protein_change = protein_change.split(':')[1]
#                     print('找到蛋白质变化信息:', protein_change)
#                     na_transcript.iloc[i,11] = f'{trans_nm}({gene}):{trans_change}({protein_change})'
#                 else:
#                     msk = df.iloc[:,-1].str.split('.').str[0].isin(transcripts_genome)
#                     transcript = df.loc[msk].iat[0,-1]
#                     trans_nm = transcript.split(':')[0]
#                     trans_change = transcript.split(':')[1]
#                     print('找到转录本信息:', transcript)
#                     p_msk = df.iloc[:,-1].str.startswith('NP_')
#                     protein_change = df.loc[p_msk].iat[0, -1]
#                     protein_change = protein_change.split(':')[1]
#                     print('找到蛋白质变化信息:', protein_change)
#                     na_transcript.iloc[i,11] = f'{trans_nm}({gene}):{trans_change}({protein_change})'
#         except Exception as e:
#             print('没有找到转录本信息', e)

#     na_transcript.to_csv('var_transcript_data/fill_na_transcript.csv', index=False)

# 运行了这个之后，手动NCBI搜索查找转录本信息核对
# get_trans_details(na_transcript)  
####实在是找不到的从alleles里面找
####################################################
data = pd.read_csv('var_transcript_data/fill_na_transcript.csv')
for i in range(data.shape[0]):
    # print(data.iloc[i,11])
    if pd.isna(data.iloc[i,11]):
        try:
            gene = data.iloc[i,7].split('=')[1]
            pos = str(data.iloc[i]['POS'])
            with open(f'alleles/{gene}_translation.json') as f:
                pos_data = json.load(f)
                gene = pos_data['gene']
                chromosome = pos_data['chromosome']
                refSeqChromosomeId = pos_data['refSeqChromosomeId']
                for var in pos_data['variants']:
                    if str(var['position']) == pos:
                        chromosomeHgvsName = var['chromosomeHgvsName']
                        data.iloc[i,11] = f"{refSeqChromosomeId}:{chromosomeHgvsName}"
        except Exception as e:
            print(f'没有找到转录本信息{gene}:{pos}', e)
data.to_csv('var_transcript_data/fill_na_transcript_919.csv', index=False)