import pandas as pd
import re

# 读取recommendations.csv文件
def merge_recommendations_drug():
    recommendations_path = 'reporter/recommendations.csv'
    recommendations_df = pd.read_csv(recommendations_path)

    # 读取药物修改.xlsx文件
    drug_modification_path = 'reporter/药物修改.xlsx'
    drug_modification_df = pd.read_excel(drug_modification_path)

    # 给药物修改数据的第一列开头和结尾加上p标签
    first_col_name = drug_modification_df.columns[0]
    drug_modification_df[first_col_name] = '<p>' + drug_modification_df[first_col_name].astype(str) + '</p>'

    # 对两个表格进行merge，这里默认使用inner join，可根据实际需求修改
    merged_df = pd.merge(recommendations_df, drug_modification_df, how='left',left_on='short_description',right_on=first_col_name)
    output_path = 'reporter/merged_data.csv'
    merged_df.to_csv(output_path, index=False)
    
def merge_guidelie_shorts():
    guidelie_path = 'reporter/药物指南翻译.xlsx'
    guidelie_df = pd.read_excel(guidelie_path)
    merged_df = pd.read_csv('reporter/merged_data.csv')
    # 给药物指南翻译数据的第一列开头和结尾加上p标签
    first_col_name = guidelie_df.columns[0]
    guidelie_df[first_col_name] = '<p>' + guidelie_df[first_col_name].astype(str) + '</p>'
    # 对两个表格进行merge，这里默认使用outer join
    merged_df = pd.merge(merged_df, guidelie_df, how='outer',left_on='text',right_on=first_col_name)
    output_path = 'reporter/results_guidelie_shorts.csv'
    merged_df.to_csv(output_path, index=False)
    
# merge_guidelie_shorts()
def process_empty_chinese_summary():
    # 读取 reporter/results_guidelie_shorts.csv 文件
    file_path = 'reporter/results_guidelie_shorts.csv'
    df = pd.read_csv(file_path)
def clean_string(text):
    """只保留A-Za-z0-9的字符，其他字符全都删掉"""
    if pd.isna(text):
        return text
    text = re.sub(r'[^A-Za-z0-9]', '', str(text))
    text = text.replace('Seelabelformoreinformation', '')
    return text

# 修改merge_recommendations_drug函数
def merge_recommendations_drug():
    recommendations_path = 'reporter/recommendations.csv'
    recommendations_df = pd.read_csv(recommendations_path)

    drug_modification_path = 'reporter/药物修改.xlsx'
    drug_modification_df = pd.read_excel(drug_modification_path)

    first_col_name = drug_modification_df.columns[0]
    drug_modification_df[first_col_name] = '<p>' + drug_modification_df[first_col_name].astype(str) + '</p>'

    # 处理merge的列，只保留A-Za-z0-9的字符
    recommendations_df['short_description'] = recommendations_df['short_description'].apply(clean_string)
    drug_modification_df[first_col_name] = drug_modification_df[first_col_name].apply(clean_string)

    merged_df = pd.merge(recommendations_df, drug_modification_df, how='left',left_on='short_description',right_on=first_col_name)
    output_path = 'reporter/merged_data.csv'
    merged_df.to_csv(output_path, index=False)

# 修改merge_guidelie_shorts函数
def merge_guidelie_shorts():
    guidelie_path = 'reporter/药物指南翻译.xlsx'
    guidelie_df = pd.read_excel(guidelie_path)
    merged_df = pd.read_csv('reporter/merged_data.csv')

    first_col_name = guidelie_df.columns[0]
    guidelie_df[first_col_name] = '<p>' + guidelie_df[first_col_name].astype(str) + '</p>'

    # 处理merge的列，只保留A-Za-z0-9的字符
    merged_df['short_description'] = merged_df['short_description'].apply(clean_string)
    guidelie_df[first_col_name] = guidelie_df[first_col_name].apply(clean_string)

    merged_df = pd.merge(merged_df, guidelie_df, how='outer',left_on='short_description',right_on=first_col_name)
    output_path = 'reporter/results_guidelie_shorts.csv'
    merged_df.to_csv(output_path, index=False)
    # 找到中文总结列为空的数据
    # empty_summary_df = df[df['中文总结'].isna()]
    
    # # 对空数据进行 groupby 操作，这里假设以short_description列本身进行分组
    # grouped_df = empty_summary_df.groupby('short_description')
    
    # # 取这一列数据
    # empty_summary_series = grouped_df['short_description'].first()
    
    # # 输出到文件中
    # output_path = 'reporter/empty_chinese_summary.csv'
    # empty_summary_series.to_csv(output_path, index=False)

# 调用函数
merge_guidelie_shorts()

