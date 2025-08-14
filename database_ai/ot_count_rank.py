import pandas as pd
from pathlib import Path
from collections import Counter


def count_rank(nc_no: str) -> None:
    ot_res_li = Path("/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/flashfry_out/").glob(f"{nc_no}*tsv")
    # Counter
    rank_counter = Counter()
    # read tsv
    type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string"]
    # type_li = ["string", "int32", "i", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string"]
    for res in ot_res_li:
        df = pd.read_csv(res,sep="\t",usecols=[0, 5, 6, 7])
        for contig, sub_df in df.groupby("contig"):
            sub_df = sub_df[sub_df["orientation"]=="FWD"]
            ot_count = sum(sub_df["otCount"]) - 1
            of_li = sub_df["overflow"].to_list()
            if ot_count > 99 or "OVERFLOW" in of_li:
                rank_counter['100+'] += 1
                continue
            rank_start = ot_count // 10 * 10
            rank_end = ot_count // 10 * 10 + 9
            rank_counter[f"{rank_start}-{rank_end}"] += 1
            # print(ot_count)
            # print(rank_start)
            # print(rank_end)
    rank_df = pd.DataFrame([rank_counter])
    header = ['0-9',  '10-19',  '20-29',    '30-39',  '40-49',   '50-59',  '60-69',  '70-79',  '80-89', '90-99',  '100+']
    rank_df = rank_df[header]
    rank_df.reset_index()
    rank_df.index=[nc_no]
    # print(rank_df)
    return rank_df
    # rank_df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/score/CFD_Scoring/test/rank/{nc_no}.xls",sep="\t",index=False)

def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    
    res = pd.DataFrame([])
    for nc in nc_li:
        print(f"{nc} satrt !!!")
        sub_res = count_rank(nc)
        res = pd.concat([res,sub_res])
        print(f"{nc} finished !!!")
    # sort by nc_li
    res = res.sort_index()
    
    res.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/all_ot_rank.xls",sep="\t",index=True)

if __name__ == "__main__":
    main()