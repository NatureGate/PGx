import json
with open('data/test_new.report.json','r',encoding='utf-8') as f:
    data = json.load(f)
cpic_drugs = data['drugs']['CPIC Guideline Annotation']
dpwg_drugs = data['drugs']['DPWG Guideline Annotation']
fda_drugs = data['drugs']['FDA Label Annotation']
pgx_drugs = data['drugs']['FDA PGx Association']
cpic_drugs = set(cpic_drugs.keys())
dpwg_drugs = set(dpwg_drugs.keys())
fda_drugs = set(fda_drugs.keys())
pgx_drugs = set(pgx_drugs.keys())

print(len(cpic_drugs))
print(len(dpwg_drugs))
print(len(fda_drugs))
print(len(pgx_drugs))
print((cpic_drugs.union(dpwg_drugs).union(fda_drugs).union(pgx_drugs)))
