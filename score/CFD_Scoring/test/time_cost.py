import os
import time
import pandas as pd
from os import system
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure

import psutil

def gdb2fa(gdb_df: pd.DataFrame, nc_no: str) -> None:
    seq = gdb_df[1].to_numpy() + gdb_df[2].to_numpy()
    # 切割序列文件 88w一切
    seq_li = [">"] * len(gdb_df[0]) + gdb_df[8].to_numpy() + ["_"] * len(gdb_df[0]) + gdb_df[0].astype(str).to_numpy() + ["\n"] * len(gdb_df[0]) + seq
    chunk_size = 880000  # 每块的大小（88w）

    # 切割列表
    chunks = [seq_li[i:i + chunk_size] for i in range(0, len(seq_li), chunk_size)]
    
    for idx, chunk in enumerate(chunks):
        with open(f"/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/fa/{nc_no}_{idx + 1}.fa", "w") as f:
            f.write("\n".join(chunk))


def flashfry_cmd(gdb_fa_path: str, ot_path: str, output_path: str) -> None:
    # print(f"java -Xmx8g -jar /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/FlashFry-assembly-1.15.jar  discover  --database {ot_path}  --fasta {gdb_fa_path}  --output {output_path} --maxMismatch 3 --maximumOffTargets 100")
    system(f"java -Xmx8g -jar /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/FlashFry-assembly-1.15.jar  discover  --database {ot_path}  --fasta {gdb_fa_path}  --output {output_path} --maxMismatch 3 --maximumOffTargets 100 --forceLinear")

# 将gdb的y染色体部分 分别扫各个染色体的ot库
def gdb2ot(output_fa: str, ot_path: str, output_path: str) -> pd.DataFrame:
    flashfry_cmd(output_fa, ot_path, output_path)
    

def map(nc_no:str):
    # output_fa = f"/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/fa/{nc_no}.fa"
    # output_fa = f"/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/fa/NC_000024.10.fa"
    output_fa = f"/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/fa/{nc_no}.fa"
    ot_path = "/mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/NCA_cas9_db"
    res_path = "/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/all_res/"
    
    res_out_path = f"{res_path}/{nc_no}.tsv"

    flashfry_cmd(output_fa,ot_path,res_out_path)
    
def run_map(nc_no: str) -> None:
    t1 = time.time()
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/cds_mark/{nc_no}.tsv"
    header_type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category"]
    gdb_df = tsv2df(gdb_path, header_type_li)
    print(f"load gdb<{nc_no}> time cost:{time.time() - t1}")
    t1 = time.time()
    
    gdb2fa(gdb_df,nc_no)
    print(f"gdb2fa <{nc_no}> time cost:{time.time() - t1}")
    t1 = time.time()  
    # map(nc_no)
    

    print(f"run flashfry <{nc_no}> time cost:{time.time() - t1}")
    peak_memory_gb = memory_info.peak_wset / (1024**3) if hasattr(memory_info, 'peak_wset') else memory_info.rss / (1024**3)

    print(f"process {nc_no} peak memory cost: {peak_memory_gb:.2f} GB")
    
    
    
    
def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # nc_li = nc_li[1:-1]
    # run_map(nc_li[0])
    async_in_iterable_structure(run_map,nc_li,24)
    
    
if __name__ == "__main__":
    main()