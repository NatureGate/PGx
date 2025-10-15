import pandas as pd
from bs4 import BeautifulSoup
import re

def process_suggest(text):
        if pd.isna(text):
            return 'pp'
        text = str(text).strip()
        text = text.replace('&quot;', '"').replace('&lt;', '').replace('&gt;', '')
        tmp_text = text
        if str(text).startswith('<p>'):
            index = text.index('</p>') 
            soup = BeautifulSoup(text, 'html.parser')
            first_p = soup.find('p')
            if first_p:
                text = first_p.get_text()
                 
                text = text+tmp_text[index+3:]      
            
        # elif text.__contains__('ul'):
        #     text = text[:text.index('ul')]
        
        # 只保留a-zA-Z0-9的字符
        filtered_text = re.sub(r'[^a-zA-Z0-9]', '', str(text))
        if 'h4idother' in filtered_text:
            index = filtered_text.index('h4idother')
            filtered_text = filtered_text[:index]
        filtered_text = filtered_text.replace('Seelabelformoreinformation', '')
        # 在字符前后加上p字符
        return f"p{filtered_text}p"
# 读取 reporter 下的 results_guidelie_shorts_grouped.xlsx 文件
file_path = 'reporter/results_guidelie_shorts_grouped.xlsx'

df = pd.read_excel(file_path)
print(df.shape)
# 找到 Text 列中包含 <ul> 的数据
result = df[df['Text'].str.contains('<ul>', na=False)]
# 对 results 数据的 Text 列进行 process_suggest 处理
result['Processed_Text'] = result['Text'].apply(process_suggest)
print(result.shape)
# 初始化一个空列表用于存储需要新增的行
new_rows = []

# 遍历 result 数据，比对 Processed_Text 和 short_description
for index, row in result.iterrows():
    
    if row['Processed_Text'] != str(row['short_description']).strip():
        # print(index, row['Processed_Text'], row['short_description'])
        new_row = row.copy()
        new_row['short_description'] = row['Processed_Text']
        new_rows.append(new_row)
    else:
        print(index, row['Processed_Text'], row['short_description'])
print(len(new_rows))
# 将新增的行添加到 result 中
if new_rows:
    new_df = pd.DataFrame(new_rows)
    result = pd.concat([df, new_df], ignore_index=True)
    # short_description	Related Chemicals	Guideline Name	Name	Text	Implications	Lookup Key	text	guide	advice	原始内容	中文总结	Guideline Type

    result = result[['short_description', 'Related Chemicals', 'Guideline Name', 'Name', 'Text', 'Implications', 'Lookup Key', 'text', 'guide', 'advice', '中文总结', 'Guideline Type']]
    result.to_excel('reporter/new_results_guidelie_shorts_grouped.xlsx', index=False)
print(result.shape)





