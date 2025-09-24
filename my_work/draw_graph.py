# @Time : 2023/3/24 11:16
# @Author : 龙锐 942121483@qq.com
# @File : draw_graph.py
# @desc :
import matplotlib.pyplot as plt
import networkx as nx
g = nx.read_gml("pre_filt_graph.gml",label='annotation')
#annotations = nx.get_node_attributes(g,'annotation')
nx.draw_networkx(g,node_size=12,font_size=1,alpha=0.75,width=0.25,pos=nx.spring_layout(g))
#plt.figure(figsize=(25, 25))
plt.savefig("ba.pdf",format='pdf')
plt.show()