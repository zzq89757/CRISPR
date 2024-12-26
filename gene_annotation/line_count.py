from multiprocessing import Pool, RLock
from os import system

import pandas as pd

import subprocess

def count_lines(file_path: str) -> str:
    try:
        command = f"cut -f1 {file_path} | sort -u | wc -l"
        
        # 使用 subprocess 运行管道命令
        result = subprocess.run(
            command,
            shell=True,  # 允许通过 shell 执行管道命令
            capture_output=True,
            text=True,
            check=True
        )
        # 提取输出结果并转换为整数
        line_count = result.stdout.strip()
        return line_count
    except Exception as e:
        print(f"统计行数时发生错误: {e}")
        return -1

def async_in_iterable_structure(fun, iterable_structure, cpus):
    def init(l):
        global lock
        lock = l

    lock = RLock()
    p = Pool(int(cpus), initializer=init, initargs=(lock,))
    # apply async in iterable structure
    for i in iterable_structure:
        p.apply_async(func=fun, args=(i,))
    p.close()
    p.join()

def run_cmd(i):
    # system(f"awk '$4==\"{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_{i}.tsv")
    # print(f"cut -f1 spCas9_Homo_{i}.tsv | sort -u|wc -l")
    # system(f"cut -f1 spCas9_Homo_{i}.tsv | sort -u|wc -l && echo {i}")
    file= open(f"line_count/spCas9_Homo_{i}.tsv",'w')
    line_num = count_lines(f"spCas9_Homo_{i}.tsv")
    print(f"{line_num}",file=file,end="\t")
    print(f"spCas9_Homo_{i}.tsv",file=file,end="\n")
    # system(f"echo spCas9_Homo_{i}.tsv >> line_count/spCas9_Homo_{i}.tsv")
    # system(f"cut -f1 spCas9_Homo_{i}.tsv | sort -u|wc -l >> line_count/spCas9_Homo_{i}.tsv")


# for i in range(1,23):
#     print(f"awk '$4==\"chr{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_chr{i}.tsv")

# for i in ['X', 'Y']:
#     print(f"awk '$4==\"chr{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_chr{i}.tsv")


def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    gtf_file = "GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.gtf"
    # gtf_file = "ZNF568.gtf"
    # num_count = Counter()
    # info_dict = defaultdict()
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    # nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    # nc_li = nc_df[0].tolist()
    chr_li = ["chr" + str(x) for x in (list(range(1, 23)) + ["X", "Y"])]
    # print(f"cat {" ".join([f"spCas9_Homo_{x}.tsv" for x in chr_li])} > spCas9_Homo\(WGS\)_gRNA-Gene\(NM\&NR\)")
    # system(f"cat {" ".join([f"spCas9_Homo_{x}.tsv" for x in chr_li])} > spCas9_Homo\(WGS\)_gRNA-Gene\(NM\&NR\)")
    # print(f"cat {" ".join([f"spCas9_Homo_{x}.tsv" for x in chr_li])} > spCas9_Homo\(WGS\)_gRNA-Gene\(NM\&NR\).tsv")
    # system(f"cat {" ".join([f"spCas9_Homo_{x}.tsv" for x in chr_li])} > spCas9_Homo\(WGS\)_gRNA-Gene\(NM\&NR\).tsv")
    # for i in chr_li:
        # print(f"wc -l spCas9_Homo_{i}.tsv")
        # print(f"spCas9_Homo_{i}.tsv uniq : ",end="")
        # system(f"cut -f1 spCas9_Homo_{i}.tsv | sort -u | wc -l")
    # print(f"wc -l spCas9_Homo\(WGS\)_gRNA-Gene\(NM\&NR\).tsv")
    # system(f"wc -l spCas9_Homo\(WGS\)_gRNA-Gene\(NM\&NR\).tsv")
    # print(f"cut -f1 spCas9_Homo\(WGS\)_gRNA-Gene\(NM\&NR\).tsv | sort -u | wc -l")
    # print(f"spCas9_Homo(WGS)_gRNA-Gene(NM&NR).tsv uniq : ",end="")
    # system(f"cut -f1 spCas9_Homo\(WGS\)_gRNA-Gene\(NM\&NR\).tsv | sort -u | wc -l")
    # async_in_iterable_structure(run_cmd,chr_li,24)
    system(f"cat {" ".join([f"line_count/spCas9_Homo_{x}.tsv" for x in chr_li])} > line_count/all")   

if __name__ == "__main__":
    main()