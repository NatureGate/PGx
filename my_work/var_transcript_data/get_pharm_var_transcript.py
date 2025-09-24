import pandas as pd
from pathlib import Path
from get_transcript_by_worm import get_transcript_by_worm
folder = Path(r'var_transcript_data')   # 1. 把路径换成你的

# 2. 想读哪些扩展名就写哪些
exts = [ '*_NM.txt']
transcripts_genome = ['NM_004827','NM_000069','NM_000492','NM_000767','NM_000771','NM_000769','NM_000106','NM_017460','NM_001360016','NM_000110','NM_001082','NM_000777','NM_172139','NM_018283','NM_000540','NM_006446','NM_000367','NM_000463','NM_024006']
frames = {}
for ext in exts:
    for file in folder.glob(ext):
        print('reading', file.name)
        df = pd.read_table(file)
        frames[file.stem] = df          # 用文件名（无扩展名）当 key

# 3. 如果还想拼成一张大表
big_table = pd.concat(frames).reset_index(level=0).astype(str)
target = pd.read_table('var_transcript_data/new_position.txt').astype(str)
print(big_table.columns)
# .to_csv('var_transcript_data/merged_all.csv', index=False)
all_transcript = pd.merge(target, big_table, left_on='POS', right_on='GRCh38Location', how='left')
mask = all_transcript['Name'].isna()
na_transcript = all_transcript[mask]
# na_transcript = pd.read_table('var_transcript_data/na_transcript.csv',sep=',').astype(str)
# mask = na_transcript['Name'].isna()
# na_transcript = na_transcript[mask]

def get_trans_details(na_transcript:pd.DataFrame):
    for i in range(na_transcript.shape[0]):
        try:
            rs_id = na_transcript.iloc[i]['ID']
            gene = na_transcript.iloc[i,7].split('=')[1]
            ALT = na_transcript.iloc[i]['ALT']
            REF = na_transcript.iloc[i]['REF']
            pos = na_transcript.iloc[i]['CHROM'][3:] + ':' + str(na_transcript.iloc[i]['POS'])
            change = na_transcript.iloc[i]['REF'] + '>' + na_transcript.iloc[i]['ALT']
            query = f"{rs_id}" if rs_id !='.' else f"?term=/{pos}"
            print(f'没有找到转录本信息，尝试从NCBI爬取: {rs_id} {pos}')
        # 4. 这里可以调用爬虫脚本
            df = get_transcript_by_worm(query)
            if df is not None:
                transcript = df.loc[df[ALT].str.split('.').str[0].isin(transcripts_genome), ALT].iat[0]
                trans_nm = transcript.split(':')[0]
                trans_change = transcript.split(':')[1]
                print('找到转录本信息:', transcript)
                protein_change = df.loc[df[ALT].str.startswith('NP_'), ALT].iat[0]
                protein_change = protein_change.split(':')[1]
                print('找到蛋白质变化信息:', protein_change)
                na_transcript.iloc[i,11] = f'{trans_nm}({gene}):{trans_change}({protein_change})'
        except Exception as e:
            print('没有找到转录本信息', e)

    na_transcript.to_csv('var_transcript_data/na_transcript.csv', index=False)

get_trans_details(na_transcript)  