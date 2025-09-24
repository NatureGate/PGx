import subprocess
import concurrent.futures
import argparse



def run_shell_script(script_path, chr, allsite_vcf_file, allsite_vcf_file_index):    
    """
    Run a shell script with the provided parameters.
    """
    command = [
        'bash', script_path,
        chr,
        allsite_vcf_file,
        allsite_vcf_file_index
    ]
    print(command)

    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"Script executed successfully: {result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Script execution failed:{chr} \n {e.stderr}")

def main():
    argparser = argparse.ArgumentParser(description="Run a shell script with multiple parameters in parallel.")
    argparser.add_argument('--vcf_file', type=str, required=True, help='Path to the VCF file')
    argparser.add_argument('--vcf_file_index', type=str, required=True, help='Path to the VCF file index', default="")
    args = argparser.parse_args()
    allsite_vcf_file = args.vcf_file
    allsite_vcf_file_index = args.vcf_file_index
    print(f"allsite_vcf_file: {allsite_vcf_file}, allsite_vcf_file_index: {allsite_vcf_file_index}")

    # 要执行的脚本路径
    script_path = "/home/stereonote/script/split.sh"
    tasks = []
    # 将所有参数组合成任务
    chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
    #chrs = ['chr22','chrX']
    for idx, chr in enumerate(chrs):
        tasks.append((chr, allsite_vcf_file, allsite_vcf_file_index))
    print(f'tasks:{tasks}')
    # 使用线程池并行执行任务
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
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