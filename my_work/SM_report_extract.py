import pandas as pd
from bs4 import BeautifulSoup

# 读取 HTML 文件
with open('SM.processed.report.html', 'r', encoding='utf-8') as file:
    html_content = file.read()

# 使用 BeautifulSoup 解析 HTML 内容
soup = BeautifulSoup(html_content, 'html.parser')

# 查找 Section II 的内容
section_ii = soup.find('section', id='section-ii')

# 初始化数据列表
data = []

# 提取每个药物部分的数据
for guideline in section_ii.find_all('section', class_='guideline'):
    drug_name = guideline.find('h3').text.strip() if guideline.find('h3') else ''

    # 遍历每个表格行
    for table_row in guideline.find_all('tr'):
        # 跳过表头行
        if table_row.find('th'):
            continue

        cells = table_row.find_all('td')
        if len(cells) < 5:
            continue  # 跳过不符合要求的行

        # 提取每个单元格的数据
        source_cell = cells[0]
        source = source_cell.find('b').text.strip() if source_cell.find('b') else ''
        population = source_cell.find('p').text.strip() if source_cell.find('p') else ''

        tags = [tag.text.strip() for tag in source_cell.find_all('div', class_='tag')]
        genes_info = cells[1].text.strip()

        implications = cells[2].text.strip() if len(cells) > 2 else ''
        recommendation = cells[3].text.strip() if len(cells) > 3 else ''
        classification = cells[4].text.strip() if len(cells) > 4 else ''

        # 合并多行推荐内容
        if len(cells) > 5:
            recommendation = ' '.join([cells[3].text.strip()] + [cell.text.strip() for cell in cells[4:-1]])
            classification = cells[-1].text.strip()

        # 将数据添加到列表
        data.append([
            drug_name, source, population, genes_info, implications,
            recommendation, classification, ', '.join(tags)
        ])

# 创建 DataFrame
df = pd.DataFrame(data, columns=[
    'Drug', 'Source', 'Population', 'Genes', 'Implications',
    'Recommendation', 'Classification', 'Tags'
])

# 保存为 CSV 文件
csv_file_path = "section_ii_data.csv"
df.to_csv(csv_file_path, index=False, encoding='utf-8')

print(f"Section II data has been successfully written to {csv_file_path}")