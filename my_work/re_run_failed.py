import subprocess
from pathlib import Path


with open('failed.txt', 'r') as file:
    samples = [line.strip() for line in file]


command_template = (
    "/shapeit4/shapeit4/bin/shapeit4.2 "
    "--input {input_file} "
    "--output {output_file} "
    "--reference /zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr6.filtered.shapeit2-duohmm-phased.vcf.gz "
    "--sequencing "
    "--region chr6 "
    "--threads 1"
    "--map /data/input/W202506100012705_1/PGx-allsite-test/call-split_vcf/shard-5/execution/chr6.b38.gmap.gz "
)


for sample in samples:
    folder_path = Path("/data/input/{sample}_1")
    if not folder_path.exists():
        continue
    chr6_vcf_file='/data/input/{sample}_1/PGx-allsite-test/call-split_vcf/shard-5/execution/chr6.vcf.gz'
    input_file = chr6_vcf_file           
    output_file = f"phased/{sample}.phased.vcf"  
    
   
    command = command_template.format(
        input_file=input_file,
        output_file=output_file
    )
    
    print(f"\nProcessing sample: {sample}")
    print(f"Executing command: {command}")
    
    result = subprocess.run(
        command,
        shell=True,
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        print(f"Error occurred for {sample}:")
        print("ERROR MESSAGE:")
        print(result.stderr)
    else:
        print(f"Successfully phased {sample}")
        print("STDOUT:")
        print(result.stdout)

print("\n success")