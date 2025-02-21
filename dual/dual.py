from pathlib import Path
import pandas as pd
from sys import path
from os import system
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


def dual(raw_db: str) -> pd.DataFrame:
    # 读取filter20 数据库
    df = tsv2df(raw_db,[])
    # 按照基因分组 
    for gene, sub_df in df.groupby(9):
        grna_num = len(sub_df)
        # 若只有一条 无法成对 跳过
        if grna_num == 1:
            continue
        # 若只有两条 强行配对并输出信息
        if grna_num == 2:
            print(sub_df)
        # 两条以上的情况
        
            
    exit()


def run_dual(nc_no: str) -> None:
    raw_db = f"/mnt/ntc_data/wayne/Repositories/CRISPR/filter_20/{nc_no}.tsv"
    output = f"/mnt/ntc_data/wayne/Repositories/CRISPR/dual/{nc_no}.tsv"
    if Path(output).exists():
        print(f"{nc_no} exists !!!")
        return
    new_df = dual(raw_db)
    # 保存为tsv
    new_df.to_csv(output, sep="\t", header=None, index=None)

def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # async_in_iterable_structure(run_dual, nc_li, 24)
    run_dual(nc_li[-1])


if __name__ == "__main__":
    main()