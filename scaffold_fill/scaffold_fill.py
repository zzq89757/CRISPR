from pathlib import Path
import pandas as pd
from sys import path
from os import system
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


def scaffold_fill(raw_db: str) -> pd.DataFrame:
    # 读取标记Low后的原始文件
    df = tsv2df(raw_db, type_li=[])
    # 选择gRNA name,gene id,guide seq 列并添加scaffold和full sequence列
    df = df[[0, 10, 2]]
    df[11] = ['GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'] * len(df)
    df[12] = df[2] + df[11]
    return df


def run_fill(nc_no: str) -> None:
    raw_db = f"/mnt/ntc_data/wayne/Repositories/CRISPR/low_mark/{nc_no}.tsv"
    output = f"/mnt/ntc_data/wayne/Repositories/CRISPR/scaffold_fill/{nc_no}.tsv"
    if Path(output).exists():
        print(f"{nc_no} exists !!!")
        return
    new_df = scaffold_fill(raw_db)
    # 保存为tsv
    new_df.to_csv(output, sep="\t", header=None, index=None)


def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_fill, nc_li, 24)
    # run_fill(nc_li[-1])
    system(f"cat NC*tsv > all.tsv")


if __name__ == "__main__":
    main()