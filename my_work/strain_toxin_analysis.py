# @Time : 2023/3/22 22:28
# @Author : 龙锐 942121483@qq.com
# @File : strain_toxin_analysis.py
# @desc :

import pandas as pd

data = pd.read_table("strain_toxin_filter.txt")
print(data.head())
# count_data = data.groupby("strain", as_index=False).count()
# sort_data = count_data.sort_values("toxin", ascending=False)
# sort_data.to_excel("sort_data.xlsx",index=False)

grouped_data = data.groupby("strain", as_index=False)
lasta_data = grouped_data['toxin'].apply(lambda x: ",".join(x))


def get_size(row):
    return len(row[1].split(","))


lasta_data['size'] = lasta_data.apply(get_size,axis=1)
sort_data = lasta_data.sort_values("size", ascending=False)
sort_data.to_excel("sort_data.xlsx",index=False)

