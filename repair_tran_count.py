import time
import pandas as pd
from collections import defaultdict
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from cds_mark import gene_ori_dict
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure
from os import system


def repair(nc_no: str, dir_path: str) -> None:
    gdb_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/{dir_path}/{nc_no}.tsv"
    tmp_gdb_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/{dir_path}/{nc_no}.tmp.tsv"
    df = tsv2df(gdb_file,[]) if dir_path != "utr_mark" else pd.read_csv(gdb_file,sep="\t",header=None,dtype={24: str})
    left_numbers = df[14].apply(lambda x: len(x.split(';')))
    df[16] = left_numbers.astype(str) + '/' + df[16].str.split('/').str[1]
    df.to_csv(tmp_gdb_file,sep="\t",header=None,index=False)
    system(f"mv {tmp_gdb_file} {gdb_file}")
    

def run_repair(nc_no: str) -> None:
    dir_li = ["tran_count", "ag_mark", "pam_filter", "cfd_filter", "flank_fill", "az_score", "snp_mark", "utr_mark"]
    # for dir in dir_li:
    #     repair(nc_no,dir)
    repair(nc_no,dir_li[-1])


def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_repair,nc_li,24)
    # run_repair(nc_li[-1])
    

if __name__ == "__main__":
    main()