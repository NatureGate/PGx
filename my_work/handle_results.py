import pandas as pd
import re
import datetime
from bs4 import BeautifulSoup
import json
import os

def merge_drug_data():
    
    # 读取drugs_1013.csv文件
    df_drugs = pd.read_csv('reporter/drugs_1013.csv', sep='\t')
    
    # 读取drugs_union.json文件
    with open('drugs_union.json', 'r', encoding='utf-8') as f:
        drugs_json = json.load(f)
    jsondrugs_en = set(drugs_json.keys())
    df_drugs_en = set(df_drugs['英文名称'])
    # 找出drugs_union.json中没有的英文名称
    en_not_in_drugs = jsondrugs_en - df_drugs_en
    en_not_in_drugs = list(filter(lambda x: not (x.__contains__(' / ') or x.__contains__(' and ')), en_not_in_drugs))
    # df = pd.DataFrame('', index=range(n), columns=[f'c{i}' for i in range(m)])
    # pd.DataFrame(columns=['英文名称', '中文名称'])
    df = pd.DataFrame('', index=range(len(en_not_in_drugs)), columns=['英文名称', '中文名称'])
    for i in range(len(en_not_in_drugs)):
        df.iloc[i, 0] = en_not_in_drugs[i]
        df.iloc[i, 1] = drugs_json[en_not_in_drugs[i]]
    # 合并数据，保留drugs_union.json中drugs_1013.csv没有的数据
    merged_df = pd.concat([df_drugs, df], ignore_index=True)
    
    # 重新输出drugs_1013.csv
    merged_df.to_csv('reporter/drugs_1013.csv', sep='\t', index=False)
    

# 调用合并函数
# merge_drug_data()




def process_data():
    # 读取Test.pgx.tsv文件
    try:
        df_tsv = pd.read_csv('Test.pgx.tsv', sep='\t')
    # 读取reporter/drugs_1013.xlsx文件，假设该文件包含药物中英文对照
   
        df_drugs = pd.read_csv('reporter/drugs_1013.csv',sep='\t')
        print(df_drugs.columns)
    except FileNotFoundError:
        print("未找到reporter/drugs_1013.csv文件，请检查文件路径。")
        return
    
    merged_for_english = pd.merge(df_tsv, df_drugs, left_on='drug', right_on='中文名称', how='left')
    merged_for_english.to_csv('merged_for_english.csv',index=False)
    try:
        df_excel = pd.read_excel('reporter/results_guidelie_shorts_grouped.xlsx')
    except FileNotFoundError:
        print("未找到reporter/results_guidelie_shorts_grouped.xlsx文件，请检查文件路径。")
        return
    except Exception as e:
        print(f"读取reporter/results_guidelie_shorts_grouped.xlsx文件时出错: {e}")
        return
    
   

    # 定义处理函数
    def process_suggest(text):
        if pd.isna(text):
            return 'pp'
        text = str(text).strip()
        
        text = text.replace('&quot;', '"').replace('&lt;', '<').replace('&gt;', '>')
        if str(text).startswith('<p>'):
            soup = BeautifulSoup(text, 'html.parser')
            first_p = soup.find('p')
            if first_p:
                text = first_p.get_text()
        # elif text.__contains__('ul'):
        #     text = text[:text.index('ul')]
        ## 当text里面没有p标签时，应该直接进行下面的操作
        # 只保留a-zA-Z0-9的字符
        filtered_text = re.sub(r'[^a-zA-Z0-9]', '', str(text))
        if 'h4idother' in filtered_text:
            index = filtered_text.index('h4idother')
            filtered_text = filtered_text[:index]
        filtered_text = filtered_text.replace('Seelabelformoreinformation', '')
        # 在字符前后加上p字符
        return f"p{filtered_text}p"
        

    # 对suggest列进行处理，新增一列
    if 'suggest' in merged_for_english.columns:
        merged_for_english['processed_suggest'] = merged_for_english['suggest'].apply(process_suggest)
        merged_for_english.to_excel('Test.pgx.tsv.xlsx', index=False)
    else:
        print("Test.pgx.tsv文件中未找到suggest列")
        return

    # 假设excel文件中药物名称列名为'drug_name'
    # 文件药物名称相同，新列内容相同的数据
    merged_df = pd.merge(merged_for_english, df_excel, left_on=['英文名称','processed_suggest'], right_on=['Related Chemicals','short_description'], how='left')
    return merged_df
    

# 调用函数
result = process_data()
if result is not None:
    today = datetime.datetime.now().strftime("%Y%m%d")
    result.to_csv(f'result_{today}.csv', index=False)

