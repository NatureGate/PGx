
### 示例代码
import re
import requests
from bs4 import BeautifulSoup
import pandas as pd
pd.set_option('display.max_columns', None)  # 显示所有列
pd.set_option('display.max_rows', None)     # 显示所有行
import time

BASE_URL = "https://www.ncbi.nlm.nih.gov/snp/"



def get_transcript_by_worm(query):
    # 目标网页的URL
    url = BASE_URL + query
    print(f"Fetching URL: {url}")
    # 发起请求
    response = requests.get(url)
    time.sleep(1)  # 避免请求过于频繁
    # 检查请求是否成功
    if response.status_code == 200:
        
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # 查找包含HGVS信息的部分
        # 假设HGVS信息位于一个特定的HTML标签中，例如<p>或<div>，并且有特定的类名或ID
        # 这里需要根据实际页面结构进行调整
        if query.startswith("?term="):
            a = soup.find('a', string=re.compile(r'^rs'))
            text = a.get_text(strip=True)        # 'rs123456'
            query = f"{text}"
            return get_transcript_by_worm(query)
        else:
            tables = soup.select('div#hgvs table')
            dfs = [pd.read_html(str(tb), encoding='utf-8')[0] for tb in tables]
            big_df = pd.concat(dfs, ignore_index=True)
            # print(big_df)
            return big_df
    else:
        print(f"Failed to retrieve the webpage. Status code: {response.status_code}")


if __name__ == "__main__":
    # 示例查询
    query = "?term=6:18132163"  # 替换为你想查询的内容
    get_transcript_by_worm(query)