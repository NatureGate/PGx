# @Time : 2023/3/19 20:41
# @Author : 龙锐 942121483@qq.com
# @File : get_counts.py
# @desc :
import  pandas as pd

data = pd.read_table("province_counts.txt",header=None)
#print(data)
print((list(data[1])))