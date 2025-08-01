from os import system
import pandas as pd
from pathlib import Path
from collections import Counter
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


def ot_res2df(ot_res: str) -> pd.DataFrame:
    type_li = ["string", "category", "category", "string", "category", "category", "category", "int32","string"]
    type_dict = dict(enumerate(type_li))

    df = pd.read_csv(
        ot_res,
        sep="\t",
        header=0,
        dtype=type_dict,
    )
    return df


def run_score_cmd(nc_no: str) -> None:
    file_li = Path("/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/flashfry_filter_out/").glob(f"{nc_no}*tsv")
    for file_name in file_li:
        mem = 8
        # print(f"java -Xmx{mem}g -jar FlashFry-assembly-1.15.jar  score  --database NCA_cas9_db  --input {str(file_name)}  --output {str(file_name).replace("all_res_filter","all_score")}  --scoringMetrics doench2016cfd")
        system(f"java -Xmx{mem}g -jar /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/FlashFry-assembly-1.15.jar  score  --database /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/db_NGG/NCA_cas9_db  --input {str(file_name)}  --output {str(file_name).replace("out","score_out")}  --scoringMetrics doench2016cfd")


# 读取扫ot库的结果 获取保留的gRNA list 用于算分等操作

def filtered_gRNA_list(ot_res_df: pd.DataFrame) -> list:
    # 过滤PAM为NAG的item
    ot_res_df = ot_res_df[ot_res_df["target"].str.endswith("GG")]
    # print(f"NGG filter,current num {len(ot_res_df)}")
    # 过滤RVS
    ot_res_df = ot_res_df[ot_res_df["orientation"]=="FWD"]
    # print(f"RVS filter,current num {len(ot_res_df)}")

    # 过滤OVERFLOW
    ot_res_df = ot_res_df[ot_res_df["overflow"]=="OK"]
    # print(f"OV filter,current num {len(ot_res_df)}")

    # 过滤 N_0 (N>1)
    ot_res_df["N_0"] = ot_res_df["offTargets"].str.findall(r'(\d+)_0').apply(lambda x: sum(map(int, x)))  # 转换为整数并求和
    # 将 N_0 转换为整数类型
    ot_res_df["N_0"] = pd.to_numeric(ot_res_df["N_0"], errors="coerce")

    # N >= 2 pass
    ot_res_df = ot_res_df[ot_res_df["N_0"]==1]
    # print(f"N0 filter,current num {len(ot_res_df)}")
    del(ot_res_df["N_0"])
    # print(ot_res_df)
    # ot_res_df.to_csv("ttt.tsv",sep="\t",index=False)
    return ot_res_df


def run_filter(nc_no) -> None:
    # glob file
    ot_res_li = Path("/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/flashfry_out/").glob(f"{nc_no}*tsv")
    Path("/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/flashfry_filter_out/").mkdir(exist_ok=True,parents=True)
    for ot_res in ot_res_li:
        ot_res_df = ot_res2df(ot_res)
        filtered_gRNA_list(ot_res_df).to_csv(str(ot_res).replace("out","filter_out"),sep="\t",index=False)
    
    
    
def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # run_filter(nc_li[-1])
    async_in_iterable_structure(run_filter,nc_li,24)
    Path("/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/flashfry_filter_score_out/").mkdir(exist_ok=True,parents=True)
    # print("ot search result filtered,scoring...")
    # run_score_cmd("/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/all_res_filter_ngg/NC_000024.10_1.tsv")
    # file_li = Path("/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/all_res_filter/").glob("*tsv")
    async_in_iterable_structure(run_score_cmd,nc_li,24)
    # run_score_cmd(nc_li[-1])
    


if __name__ == "__main__":
    main()

