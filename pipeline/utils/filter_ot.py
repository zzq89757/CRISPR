from os import system
import pandas as pd
from pathlib import Path



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


def run_score_cmd(file_name) -> None:
    mem = 8
    print(f"java -Xmx{mem}g -jar FlashFry-assembly-1.15.jar  score  --database NCA_cas9_db  --input {str(file_name)}  --output {str(file_name).replace("all_res_filter","all_score")}  --scoringMetrics doench2016cfd")
    system(f"java -Xmx{mem}g -jar /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/FlashFry-assembly-1.15.jar  score  --database /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/db_NGG/NCA_cas9_db  --input {str(file_name)}  --output {str(file_name).replace("all_res_filter","all_score")}  --scoringMetrics doench2016cfd")


# 读取扫ot库的结果 获取保留的gRNA list 用于算分等操作

def filtered_gRNA_list(ot_res_df: pd.DataFrame) -> list:
    # 过滤RVS
    ot_res_df = ot_res_df[ot_res_df["orientation"]=="FWD"]
    # print(f"RVS filter,current num {len(ot_res_df)}")

    # 过滤OVERFLOW
    ot_res_df = ot_res_df[ot_res_df["overflow"]=="OK"]
    # print(f"OVERFLOW filter,current num {len(ot_res_df)}")

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


def filter_ot(project_dir: str, nc_no:str) -> None:
    # file path
    ot_res_file = f"{project_dir}/flashfry_out/raw_out/{nc_no}.tsv"
    filtered_file = f"{project_dir}/flashfry_out/filter_out/{nc_no}.tsv"
    ot_res_df = ot_res2df(ot_res_file)
    filtered_gRNA_list(ot_res_df).to_csv(filtered_file,sep="\t",index=False)
    
    
    
def main() -> None:
    nc2chr_file = "/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/nc_li"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    filter_ot("/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus", nc_li[-1])
    # async_in_iterable_structure(run_filter,nc_li,24)
    
    


if __name__ == "__main__":
    main()

