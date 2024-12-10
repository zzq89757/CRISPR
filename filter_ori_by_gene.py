import pandas as pd


def annotation_gdb(gdb: pd.DataFrame, gene_pos_df: pd.DataFrame) -> pd.DataFrame:
    # 需要分别考虑切点和起止 三个位置
    
    return pd.DataFrame([])


def main() -> None:
    nc_no = "NC_000024.10"
    type_li = ["string", "int32", "int32", "category", "int32", "category"]
    type_dict = dict(enumerate(type_li))
    # 读取gRNA db 文件 
    gdb_df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/{chr}")
    # 读取基因位置信息文件
    gene_pos_df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/{nc_no}/Gene_list.tsv", sep="\t", header=None, low_memory=False, dtype=type_dict)
    # gdb 注释
    annotated_df = annotation_gdb(gdb_df, gene_pos_df)
    # 保存为tsv文件
    annotated_df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/{chr}")
    

if __name__ == "__main__":
    main()
    