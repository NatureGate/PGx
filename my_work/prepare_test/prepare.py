import pandas as pd
# bam_msg = pd.read_excel('prepare_test/bam.xlsx')
vcf_msg = pd.read_excel('prepare_test/vcf_path.xlsx')
# bam_msg['bam_file_path'] = bam_msg['bam_dir']+'/'+bam_msg['bam']
# mask = bam_msg['bam_file_path'].str.endswith('.bam')
# bam_data = bam_msg[mask]['bam_file_path'].tolist()
# bam_index_data = bam_msg[~mask]['bam_file_path'].tolist()
# print(bam_data)
# print(bam_index_data)
# bam_data_df = pd.DataFrame({'bam':bam_data,'bam_index':bam_index_data}).to_excel('prepare_test/bam_data.xlsx',index=False)

# vcf_mask = vcf_msg['WGS_type']=='WGS-H'
# vcf_id = vcf_msg[vcf_mask]['id'].tolist()

# vcf_msg['vcf_index'] = vcf_msg['vcf']+'.tbi'
# vcf_path = vcf_msg[vcf_mask]['vcf'].tolist()
# vcf_index = vcf_msg[vcf_mask]['vcf_index'].tolist()
# vcf_data_df = pd.DataFrame({'EntityID':vcf_id,'vcf_allsite_file':vcf_path,'vcf_allsite_file_index':vcf_index}).to_excel('prepare_test/vcf_data.xlsx',index=False)

vcf_data = pd.read_excel('prepare_test/vcf_data.xlsx')
vcf_data = vcf_data.iloc[201:1001,:]
# bam_file	hla_file	bam_index	positions	pharmcat_positions	pharmcat_positions_index	SampleID	HLA_dir	lush_ref_dir	HLA_ref_dir	sourceD	dbSNPDir	dbSNP	dbSNP_index
vcf_data['bam_file'] = '/Files/ResultData/Workflow/W202508130007117/24S08310630.bqsr.bam'
vcf_data['hla_file'] = '/Files/ResultData/Workflow/W202507160008563/E-S24707256971/hla/E-S24707256971.hla_typing.txt'
vcf_data['bam_index'] = '/Files/ResultData/Workflow/W202508130007117/24S08310630.bqsr.bam.bai'
vcf_data['positions'] = '/Files/sz_history/longrui/pgx/positions.bed'
vcf_data['pharmcat_positions'] = '/Files/sz_history/longrui/pgx/pharmcat_positions.vcf.gz'
vcf_data['pharmcat_positions_index'] = '/Files/sz_history/longrui/pgx/pharmcat_positions.vcf.gz.csi'
vcf_data['SampleID'] = vcf_data['EntityID']
vcf_data['HLA_dir'] = '/Files/sz_history/huangfei/BGE/hcWGS/HLA'
vcf_data['lush_ref_dir'] = '/Files/sz_history/zhenghaihui/lush'
vcf_data['HLA_ref_dir'] = '/Files/sz_history/huangfei/BGE/database'
vcf_data['sourceD'] = 'megabolt'
vcf_data['dbSNPDir'] = '/Files/sz_history/huangfei/BGE/database/genome/hg38_noalt_withrandom'
vcf_data['dbSNP'] = '/Files/sz_history/huangfei/BGE/database/genome/hg38_noalt_withrandom/All_20180418.new.vcf.gz'
vcf_data['dbSNP_index'] = '/Files/sz_history/huangfei/BGE/database/genome/hg38_noalt_withrandom/All_20180418.new.vcf.gz.csi'
# vcf_data = pd.read_csv('F:\downloads\20251111_PGx-vcf-multisplit\20251111_PGx-vcf-multisplit.vcf',sep='\t')
print(vcf_data.shape)
vcf_data.to_excel('prepare_test/vcf_data_1112.xlsx',index=False)
