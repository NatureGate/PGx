import pandas as pd
import re

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


def process_excel(file_path):
    """
    原地修改 Excel：
    short_description == 'pnotfoundp' 且 Text != '.' 的行，
    将其 short_description 替换为 fix_short_desc(Text) 的返回值。
    """
    # 读取
    df = pd.read_excel(file_path)

    # 定位目标行
    mask = (df['short_description'] == 'pnotfoundp') & (df['Text'] != '.')

    # 应用自定义方法
    df.loc[mask, 'short_description'] = df.loc[mask, 'Text'].apply(process_suggest)

    # 写回同名文件（原地修改）
    df.to_excel('reporter/new_short_guideline_1018.xlsx', index=False)


# 示例调用
if __name__ == '__main__':
    process_excel('reporter/new_results_guidelie_shorts_grouped.xlsx')