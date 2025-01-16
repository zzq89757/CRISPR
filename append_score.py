from pathlib import Path
import pandas as pd
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure
# 联合 pam_filter 和 score/CFD_Scoring/test/all_score_ngg/ 中的结果 filter by cfd score and append cfd score


def score_res2df(nc_no: str) -> pd.DataFrame:
    file_li = sorted(list(Path(f"/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/all_score_ngg/").glob(f"{nc_no}*.tsv")),key=lambda x: int(x.stem.split('_')[-1]))
    if len(file_li) == 1:
        return pd.read_csv(file_li[0],sep="\t",usecols=["contig","DoenchCFD_specificityscore","otCount"])
    all_res = pd.DataFrame([])
    for file in file_li:
        tmp_res = pd.read_csv(file,sep="\t",usecols=["contig","DoenchCFD_specificityscore","otCount"])
        all_res = pd.concat([all_res, tmp_res])
    
    return all_res


def gdb_reset_sample_id(gdb_df: pd.DataFrame) -> None:
    new_id = gdb_df[8].to_numpy() + ["_"] * len(gdb_df[0]) + gdb_df[0].astype(str).to_numpy()
    gdb_df[0] = new_id


def run_append(nc_no: str) -> None:
    # gdb 2 df 
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/pam_filter/{nc_no}.tsv"
    type_li =  ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category"]
    gdb_df = tsv2df(gdb_path,type_li)
    
    # score result 2 df
    score_df = score_res2df(nc_no)
    
    # gdb sample id alter
    gdb_reset_sample_id(gdb_df)
    
    # filter gdb by new sample id
    gdb_df = gdb_df[gdb_df[0].isin(score_df["contig"])]
    
    # append score by new id
    merged_table = pd.merge(
    gdb_df,
    score_df,
    left_on=0,  # 表一的主键列
    right_on='contig',  # 表二的主键列
    how='left'  # 使用左连接，确保保留表一中的所有数据
)   
    del(merged_table["contig"])
    
    append_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/cfd_filter/{nc_no}.tsv"
    merged_table.to_csv(append_path,sep="\t",header=None,index=None)
    
    
    
    

def main() -> None:
    # run_append("NC_000024.10")
    # run_append("NC_000001.11")
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_append,nc_li,24)


if __name__ == "__main__":
    main()