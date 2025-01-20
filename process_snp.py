import os
import time
import pandas as pd
import psutil
from utils.read_tsv import tsv2df
import numpy as np
# 用snp文件遍历gdb起止 似乎不可行 gdb交集区太多


def multi_indel_split(snp_df: pd.DataFrame) -> np.ndarray:
    # 将多位点indel拆分成多个item 并返回pos array
    pos_array = snp_df[0].to_numpy()
    ref_array = snp_df[1].to_numpy()
    alter_array = snp_df[2].to_numpy()
    snp_pos_array = np.array(object=[],dtype=np.int64)
    # 遍历寻找multi indel 位点
    for pos,ref,alt in zip(pos_array,ref_array,alter_array):
        # 可能出现len(ref)==len(alt)!=1的情况
        if len(ref) == len(alt) == 1:
            snp_pos_array = np.append(snp_pos_array,pos)
        else:
            ...
    
    print(snp_pos_array)


# 双指针遍历 snp和gdb
def snp_detective(gdb: pd.DataFrame, snp_df: pd.DataFrame) -> None:
    # print(gdb)
    # 提取guide+pam
    
    # print(snp_df)
    multi_indel_split(snp_df)
    return


def run_snp(nc_no: str) -> None:
    t1 = time.time()
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc2chr_dict = dict(zip(nc_df[0],nc_df[1]))
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/az_score/{nc_no}.tsv"
    # gdb_path = "/mnt/ntc_data/wayne/Repositories/CRISPR/ag_mark/NC_000024.10.tsv"
    # read gdb
    header_type_li = ["string", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category"]
    gdb_df = tsv2df(gdb_path, header_type_li)
    # print(gdb_df)
    print(f"load gdb<{nc_no}> time cost:{time.time() - t1}")
    # read snp file
    snp_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/vcf_split/filter/{nc_no}.tsv"
    snp_type_li = ['int32', 'string', 'string', 'string'] 
    snp_df = tsv2df(snp_path, snp_type_li)
    # 注释gdb
    t1 = time.time()
    snp_detective(gdb_df, snp_df)
    peak_memory_gb = memory_info.peak_wset / (1024**3) if hasattr(memory_info, 'peak_wset') else memory_info.rss / (1024**3)

    print(f"process {nc_no} peak memory cost: {peak_memory_gb:.2f} GB")
    print(f"process {nc_no} time cost:{time.time() - t1}")
    



def main() -> None:
    run_snp("NC_000024.10")
    return


if __name__ == "__main__":
    main()