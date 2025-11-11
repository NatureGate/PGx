import csv
import sys
import argparse



def process_cyp2d6_file(input_file, output_file):
    # 读取TSV文件并处理数据
    try:
        with open(input_file, 'r') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            
            with open(output_file, 'w') as txtfile:
                for row in reader:
                    genotype = row['Genotype']
                    if genotype!='None':
                    # 如果genotype包含+但不包含空格+，则将+替换为空格+
                        if '+' in genotype and ' + ' not in genotype:
                            genotype = genotype.replace('+', ' + ')
                        txtfile.write(f"CYP2D6\t{genotype}")
                    else:
                        txtfile.write(f"CYP2D6\t*1/*1")
                    break  # 如果只提取第一个符合条件的行，可保留此行
    except FileNotFoundError:
        print(f"错误：找不到文件 '{input_file}'")
        sys.exit(1)
    except Exception as e:
        print(f"处理文件时出错：{e}")
        sys.exit(1)
    
if __name__ == "__main__":
    output_file = 'cyp2d6_result.txt'
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description='从 CYP2D6 类型文件中提取特定内容。')
    parser.add_argument('--input_file', help='输入文件名')
    args = parser.parse_args()
    # 调用处理函数
    process_cyp2d6_file(args.input_file, output_file)