import pandas as pd

       
data = pd.read_excel('pgx_annotate_param.xlsx')
data['merged_phased_vcf'] = data['dir'] + '/merged_phased.vcf'
data['merged_unphased_vcf'] = data['dir'] + '/merged_unphased.vcf'
data['merged_phased_vcf_gz'] = data['dir'] + '/merged_phased.vcf.gz'
data['merged_phased_vcf_gz_index'] = data['dir'] + '/merged_phased.vcf.gz.csi'
data['merged_unphased_vcf_gz'] = data['dir'] + '/merged_unphased.vcf.gz'
data['merged_unphased_vcf_gz_index'] = data['dir'] + '/merged_unphased.vcf.gz.csi'
data['outcall'] = '/Files/ResultData/Workflow/W202507240001617/outcall.tsv'
data.to_excel('pgx_annotate_param_updated.xlsx', index=False)