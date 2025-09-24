import os

base_dir = '/Files/ResultData/Notebook/NB2025062915234244639928'
with open('paramters.csv','w') as f:
    f.write('Entity ID'+','+'SampleID'+','+'FASTQ1'+','+'FASTQ2'+','+'ref_dir'+','+'OutputSortMarkdupBam'+'\n')
    for i in range(323, 333):
        entity_id = f'ERR1955{i}'
        sample_id = f'ERR1955{i}'
        fastq1 = f'"{base_dir}/{sample_id}_1.fastq.gz"'
        fastq2 = f'"{base_dir}/{sample_id}_2.fastq.gz"'
        ref_dir = '/Files/sz_history/zhenghaihui/lush'
        output_bam = 'false'
        f.write(f'{entity_id},{sample_id},{fastq1},{fastq2},{ref_dir},{output_bam}\n')