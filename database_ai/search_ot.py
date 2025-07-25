import os
from pathlib import Path
import time
import pandas as pd
from os import system
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


def gdb2fa(gdb_df: pd.DataFrame, nc_no: str) -> None:
    fa_dir = "/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/fa/"
    Path(fa_dir).mkdir(exist_ok=True,parents=True)
    seq = gdb_df[1].to_numpy() + gdb_df[2].to_numpy()
    # 切割序列文件 88w一切
    seq_li = [">"] * len(gdb_df[0]) + gdb_df[7].to_numpy() + ["_"] * len(gdb_df[0]) + gdb_df[0].astype(str).to_numpy() + ["\n"] * len(gdb_df[0]) + seq
    chunk_size = 880000  # 每块的大小（88w）

    # 切割列表
    chunks = [seq_li[i:i + chunk_size] for i in range(0, len(seq_li), chunk_size)]
    
    for idx, chunk in enumerate(chunks):
        with open(f"{fa_dir}/{nc_no}_{idx + 1}.fa", "w") as f:
            f.write("\n".join(chunk))
            f.write("\n")


def flashfry_cmd(gdb_fa_path: str, ot_path: str, output_path: str) -> None:
    # print(f"java -Xmx8g -jar /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/FlashFry-assembly-1.15.jar  discover  --database {ot_path}  --fasta {gdb_fa_path}  --output {output_path} --maxMismatch 3 --maximumOffTargets 100")
    # system(f"java -Xmx8g -jar /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/FlashFry-assembly-1.15.jar  discover  --database {ot_path}  --fasta {gdb_fa_path}  --output {output_path} --maxMismatch 3 --maximumOffTargets 100 --forceLinear")
    system(f"java -Xmx8g -jar /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/FlashFry-assembly-1.15.jar  discover  --database {ot_path}  --fasta {gdb_fa_path}  --output {output_path} --maxMismatch 3 --maximumOffTargets 100 --forceLinear")
    

def run_flashfry_cmd(nc_no:str):
    ot_path = "/mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/db_NGG/NCA_cas9_db"
    
    output_fa_li = list(Path(f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/fa/").glob(f"{nc_no}*fa"))
    
    res_path = "/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/flashfry_out/"
    Path(res_path).mkdir(exist_ok=True,parents=True)
    res_out_path_li = [f"{res_path}/{fa.name.replace("fa","tsv")}" for fa in output_fa_li]
    
    for output_fa,res_out_path in zip(output_fa_li, res_out_path_li): 
        if Path(res_out_path).exists():
            print(f"file <{res_out_path}> already exists, skip")
            continue   
        flashfry_cmd(output_fa,ot_path,res_out_path)
    
def run_map(nc_no: str) -> None:
    t1 = time.time()
    
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/gene_annotated/{nc_no}.tsv"
    header_type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category"]
    header_type_li = []
    gdb_df = tsv2df(gdb_path, header_type_li)
    print(f"load gdb<{nc_no}> time cost:{time.time() - t1}")
    t1 = time.time()
    
    gdb2fa(gdb_df,nc_no)
    print(f"gdb2fa <{nc_no}> time cost:{time.time() - t1}")
    t1 = time.time()  
    run_flashfry_cmd(nc_no)
    

    print(f"run flashfry <{nc_no}> time cost:{time.time() - t1}")
    
    
    
    
def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # run_map(nc_li[-1])
    async_in_iterable_structure(run_map,nc_li,24)
    
    
if __name__ == "__main__":
    main()