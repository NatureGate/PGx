import json

# 读取两个 JSON 文件
with open('en2zh/drugs.json', 'r', encoding='utf-8') as f1:
    d1 = json.load(f1)

with open('en2zh/drugs_other.json', 'r', encoding='utf-8') as f2:
    d2 = json.load(f2)

# 合并两个字典（后覆盖前）
merged = {**d1, **d2}

# 保存为新的 JSON 文件
with open('en2zh/drugs_union.json', 'w', encoding='utf-8') as out:
    json.dump(merged, out, ensure_ascii=False, indent=2)