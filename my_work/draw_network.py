# @Time : 2023/3/24 12:01
# @Author : 龙锐 942121483@qq.com
# @File : draw_network.py
# @desc :
import matplotlib.pyplot as plt
import networkx as nx
g = nx.read_gml("final_graph.gml")
nx.draw(g)
plt.savefig("work.png")
plt.show()