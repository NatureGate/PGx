import subprocess
import concurrent.futures
import argparse



def run_shell_script(script_path, vcf_file, vcf_file_index, chr_genetic_map, reference_panel, reference_panel_index):    
    """
    Run a shell script with the provided parameters.
    """
    command = [
        'bash', script_path,
        chr,
        vcf_file,
        
    ]

    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"Script executed successfully: {result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Script execution failed: {e.stderr}")

def main():
    argparser = argparse.ArgumentParser(description="Run a shell script with multiple parameters in parallel.")
    argparser.add_argument('--vcf_file', type=str, required=True, help='Path to the VCF file')
    argparser.add_argument('--vcf_file_index', type=str, required=True, help='Path to the VCF file index')
    argparser.add_argument('--chr_genetic_map', type=str, required=True, help=' Path to the chromosome genetic map file')
    argparser.add_argument('--reference_panel', type=str, required=True, help='Path to  the reference panel file')
    argparser.add_argument('--reference_panel_index', type=str, required=True, help='Path to the reference panel index file')
    
    args = argparser.parse_args()
    all_vcf_file = args.vcf_file
    all_vcf_file_index = args.vcf_file_index    
    all_genetic_map = args.chr_genetic_map
    all_reference_panel = args.reference_panel  
    all_reference_panel_index = args.reference_panel_index 
    all_vcf_file_list = all_vcf_file.split(',')
    all_vcf_file_index_list = all_vcf_file_index.split(',')     
    all_genetic_map_list = all_genetic_map.split(',')
    all_reference_panel_list = all_reference_panel.split(',')
    all_reference_panel_index_list = all_reference_panel_index.split(',')  
    print(f"all_vcf_file_list: {all_vcf_file_list}")
    print(f"all_vcf_file_index_list: {all_vcf_file_index_list}")  
    print(f"all_genetic_map_list: {all_genetic_map_list}")
    print(f"all_reference_panel_list: {all_reference_panel_list}")
    print(f"all_reference_panel_index_list: {all_reference_panel_index_list}") 
    # 要执行的脚本路径
    script_path = "phase.sh"
    tasks = []
    # 将所有参数组合成任务
    chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
    for idx, chr in enumerate(chrs):
        for vcf_file in all_vcf_file_list:
            if vcf_file.__contains__(f'{chr}.vcf.gz'):
                tasks.append((chr, vcf_file))
    # 使用线程池并行执行任务
    with concurrent.futures.ThreadPoolExecutor(max_workers=24) as executor:
        futures = []
        for task in tasks:
            # 提交任务到线程池
            futures.append(executor.submit(run_shell_script, script_path, *task))
        
        # 等待所有任务完成并处理结果
        for future in concurrent.futures.as_completed(futures):
            try:
                data = future.result()
            except Exception as exc:
                print(f"任务引发异常: {exc}")
    print("所有任务已完成。")
    

if __name__ == "__main__":
    main()