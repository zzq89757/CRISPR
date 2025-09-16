import pandas as pd
from utils.read_tsv import tsv2df

   
    
def run_count(nc_no: str) -> None:
    # 读取gdb
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/ag_mark/{nc_no}.tsv"
    type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category"]
    
    gdb_df = tsv2df(gdb_path, type_li)
    
    no_df= gdb_df[(gdb_df[15]=="yes") & (gdb_df[17]=="no")]
    no_len = len(set(no_df[0].to_list()))
    ag_df= gdb_df[(gdb_df[15]=="yes") & (gdb_df[17]=="AG")]
    ag_len = len(set(ag_df[0].to_list()))
    gg_df= gdb_df[(gdb_df[15]=="yes") & (gdb_df[17]=="GG")]
    gg_len = len(set(gg_df[0].to_list()))

    print(f"{no_len}\t{ag_len}\t{gg_len}",end="\t")

def main() -> None:
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()

    for idx, nc_no in enumerate(nc_li):
        run_count(nc_no)
        if idx%8 == 7:
            print("\n",end="")
    
    
if __name__ == "__main__":
    main()