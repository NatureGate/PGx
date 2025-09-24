import pandas as pd
from bs4 import BeautifulSoup
UNKOWN_GENOTYPE_GENES = ['CYP2C19', 'VKORC1',  'SLCO1B1']
# 读取 HTML 文件
with open('SM.processed.report.html', 'r', encoding='utf-8') as file:
    html_content = file.read()

# 使用 BeautifulSoup 解析 HTML 内容
soup = BeautifulSoup(html_content, 'html.parser')

# 查找 Section II
section_ii = soup.find('section', id='section-ii')

# 初始化一个列表来存储数据
data = []

# 遍历每个药物部分
for guideline in section_ii.find_all('section', class_='guideline'):
    drug_name = guideline.h3.text.strip() if guideline.h3 else ''

    # 遍历表格中的每一行
    for row in guideline.find_all('tr'):
        # 跳过表头行
        if row.find('th'):
            continue

        # 提取每个单元格的数据
        cols = row.find_all('td')
        if len(cols) < 5:
            continue  # 跳过不符合要求的行

        # 提取源信息
        source = cols[0].find('b').text.strip() if cols[0].find('b') else ''
        
        tags = [tag.text.strip() for tag in cols[0].find_all('div', class_='tag')]

        # 提取基因信息中的 Genotype 和 Phenotype
        genes_info = cols[1].text.strip()
        genotype_parts = []
        phenotype_parts = []
        spans = cols[1].find_all('span', class_='rx-dip')
        for span in spans:
            gene_genotype = span.text
            gene = gene_genotype.split(':')[0].strip()
            if gene in UNKOWN_GENOTYPE_GENES:
                continue
            genotype = gene_genotype.split(':')[1].strip()
            
        break
        # 检查基因信息是否包含 Genotype 和 Phenotype
        # if '<div class="hint">Genotype</div>' in str(cols[1]):
        #     genotypes = cols[1].find_all('li')
        #     for genotype in genotypes:
        #         genotype_parts.append(genotype.text.strip())
        #     # 查找 Phenotype
        #     phenotype_div = cols[1].find('p', class_='rx-phenotype')
        #     if phenotype_div:
        #         phenotype_parts.append(phenotype_div.text.strip())

        # 提取影响和推荐信息
        implications = cols[2].text.strip() if len(cols) > 2 else ''
        recommendation = cols[3].text.strip() if len(cols) > 3 else ''
        classification = cols[4].text.strip() if len(cols) > 4 else ''

        # 合并多行推荐内容
        if len(cols) > 5:
            recommendation = ' '.join([cols[3].text.strip()] + [cell.text.strip() for cell in cols[4:-1]])
            classification = cols[-1].text.strip()

        # 将数据添加到列表，如果有多个 Genotype 或 Phenotype，则拆分为多行
        if genotype_parts:
            for genotype in genotype_parts:
                data.append([
                    drug_name, source, population, genotype, ', '.join(phenotype_parts),
                    implications, recommendation, classification, ', '.join(tags)
                ])
        else:
            data.append([
                drug_name, source, population, genes_info, ', '.join(phenotype_parts),
                implications, recommendation, classification, ', '.join(tags)
            ])

# 创建 DataFrame
df = pd.DataFrame(data, columns=[
    'Drug', 'Source', 'Population', 'Genotype', 'Phenotype',
    'Implications', 'Recommendation', 'Classification', 'Tags'
])

# 保存为 CSV 文件
csv_file_path = "section_ii_data.csv"
df.to_csv(csv_file_path, index=False, encoding='utf-8')

print(f"数据已成功保存到 {csv_file_path}")