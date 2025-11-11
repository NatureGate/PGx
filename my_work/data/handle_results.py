from numpy import result_type
import pandas as pd
import re
import datetime
from bs4 import BeautifulSoup
import json
import os

short_guideline_path = 'reporter/new_short_guideline_1018.xlsx'
drugs_name_path = 'reporter/drugs_1013.csv'

annotation_map = {
    'CPIC Guideline Annotation':'CPIC',
    'DPWG Guideline Annotation':'DPWG',
    'FDA Label Annotation':'FDA Label',
    'FDA PGx Association':'FDA PGx'
}

phenotypes_dict = {
    "Decreased Function": "功能降低",
    "Normal Function": "功能正常",
    "Poor Function": "功能缺失",
    "Malignant Hyperthermia Susceptibility": "恶性高热易感性",
    "Uncertain Susceptibility": "易感性不确定",
    "ivacaftor non-responsive CFTR sequence": "对伊伐卡托无反应的CFTR序列",
    "ivacaftor non-responsive in CF patients": "CF患者中对伊伐卡托无反应",
    "ivacaftor responsive in CF patients": "CF患者中对伊伐卡托有反应",
    "Indeterminate": "不确定",
    "Intermediate Metabolizer": "中间代谢型",
    "Likely Intermediate Metabolizer": "很可能的中间代谢型",
    "Likely Poor Metabolizer": "很可能的不良代谢型",
    "Normal Metabolizer": "正常代谢型",
    "Poor Metabolizer": "不良代谢型",
    "Rapid Metabolizer": "快速代谢型",
    "Ultrarapid Metabolizer": "超快速代谢型",
    "Deficient": "缺乏",
    "Deficient with CNSHA": "伴有CNSHA的缺乏",
    "Variable": "变异",
    "increased risk of aminoglycoside-induced hearing loss": "氨基糖苷类药物诱导听力损失风险增加",
    "normal risk of aminoglycoside-induced hearing loss": "氨基糖苷类药物诱导听力损失风险正常",
    "uncertain risk of aminoglycoside-induced hearing loss": "氨基糖苷类药物诱导听力损失风险不确定",
    "Possible Intermediate Metabolizer": "可能的中间代谢型",
    "Decreased Function": "功能降低",
    "Increased Function": "功能增加",
    "Possible Decreased Function": "可能的功能降低"
}




def process_data(pd_results):
    # 读取Test.pgx.tsv文件
    try:
        df_tsv = pd_results
    # 读取reporter/drugs_1013.xlsx文件，假设该文件包含药物中英文对照
   
        df_drugs = pd.read_csv(drugs_name_path,sep='\t')
        print(df_drugs.columns)
    except FileNotFoundError:
        print(f"未找到{drugs_name_path}文件，请检查文件路径。")
        return
    
    merged_for_english = pd.merge(df_tsv, df_drugs, left_on='drug', right_on='name_en', how='left')
    merged_for_english.to_csv('merged_for_english.csv',index=False)
    try:
        df_excel = pd.read_excel(short_guideline_path)
        print(df_excel.columns)
    except FileNotFoundError:
        print(f"未找到{short_guideline_path}文件，请检查文件路径。")
        return
    except Exception as e:
        print(f"读取{short_guideline_path}文件时出错: {e}")
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
    merged_df = pd.merge(merged_for_english, df_excel, left_on=['name_en','processed_suggest'], right_on=['Related Chemicals','short_description'], how='left')
    return merged_df
    

# 调用函数
def get_short_guide(pd_results):
    result = process_data(pd_results)
    result.to_csv('org_result.csv',index=False)
    # sample	drug	gatkscore	depth	readsratio	gene	diplotype	cHGVS	pHGVS	zyg	guide_x	advice_x	effect	suggest	ref_guide	drugid	rsID	location	name_en	name_zh	processed_suggest	short_description	Related Chemicals	Guideline Name	Name	Text	Implications	Lookup Key	text	guide_y	advice_y	zh_guide	Guideline Type
    result = result[['sample','name_zh','gatkscore','depth','readsratio','gene','diplotype','cHGVS','pHGVS','zyg','guide','effect','advice','zh_guide','annotation_type_y','name_en','drugid','rsID','location']]
    result['effect'] = result['effect'].map(phenotypes_dict).fillna('.')
    result = result[~result['name_en'].isin(['amphetamine','hydrocodone','oliceridine','dronabinol','lesinurad','mavacamten'])]
    result.rename(columns={
        
        'name_zh': 'drug',
        'guide': 'guide',
        'advice': 'advice',
        'zh_guide':'suggest',
        'annotation_type_y':'ref_guide'
    }, inplace=True)
    result['ref_guide'] = result['ref_guide'].map(annotation_map).fillna('.')
    result = result[['sample','name_en','drug','gatkscore','depth','readsratio','gene','diplotype','cHGVS','pHGVS','zyg','guide','effect','advice','suggest','ref_guide','drugid','rsID','location']]
    cols_to_fill = ['gatkscore','depth','readsratio','cHGVS','pHGVS','zyg','rsID','location']
    result[cols_to_fill] = result[cols_to_fill].fillna('.')
    guideline_fill = ['advice','suggest']
    result[guideline_fill] = result[guideline_fill].fillna('常规用药')
    mask = result['drug'].str.contains('/') | result['drug'].str.contains(' and ')|result['drug'].isna()
    result = result[~mask]
    if result is not None:
        today = datetime.datetime.now().strftime("%Y%m%d")
        result.to_csv(f'result_{today}.csv', index=False)
        return result

    
    
