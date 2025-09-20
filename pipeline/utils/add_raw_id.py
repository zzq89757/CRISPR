import pandas as pd

def add_raw_id(nc_idx: int, lc_file_path: str, raw_db: str, output_file: str) -> None:
    # 读取lc文件并存为列表
    lc_li = pd.read_csv(lc_file_path,sep="\t",header=None)[0].to_list()
    forward_lc = sum(lc_li[:nc_idx]) + 1 # index 从0开始 +1进行补偿
    # 读取raw_db 并将index 加上 forward_lc
    raw_db_df = pd.read_csv(raw_db,sep="\t",header=None)
    raw_db_df.index += forward_lc
    # 导出为新的数据库
    raw_db_df.to_csv(output_file,sep="\t",header=None)
    


if __name__ == "__main__":
    add_raw_id(1,"/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/ref_scan/lc_all","/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/ref_scan/NC_000068.8.tsv","/mnt_data/Wayne/Repositories/CRISPR/pipeline/NC_000068.8.tsv")