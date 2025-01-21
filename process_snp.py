import time
import pandas as pd
from utils.read_tsv import tsv2df
import numpy as np
# 用snp文件遍历gdb起止 似乎不可行 gdb交集区太多


def guide_append_pam(gdb_df: pd.DataFrame) -> np.ndarray:
    return gdb_df[1].to_numpy() + gdb_df[2].to_numpy()


def multi_indel_split(snp_df: pd.DataFrame) -> np.ndarray:
    # 将多位点indel拆分成多个item 并返回pos array
    pos_array = snp_df[0].to_numpy()
    ref_array = snp_df[1].to_numpy()
    alter_array = snp_df[2].to_numpy()
    snp_pos_list = []
    # 遍历寻找multi indel 位点
    for pos,ref,alt in zip(pos_array,ref_array,alter_array):
        ref_len = len(ref)
        # 常规1->1 直接append
        if len(alt) == 1 and ref_len == 1:
            snp_pos_list.append(pos)
            continue     
        # 逗号分隔alt
        for alt_seg in alt.split(","):
            alt_len = len(alt_seg)
            duplex_len = min(ref_len,alt_len)
            full_len = max(ref_len,alt_len)
            # duplex 判断碱基是否一致 
            for i in range(duplex_len):
                if ref[i] == alt_seg[i]:continue
                snp_pos_list.append(pos + i)
            # gap 直接判定为indel
            snp_pos_list.extend(range(pos + duplex_len, pos + full_len))
    return np.unique(snp_pos_list)


# 双指针遍历 snp和gdb
def snp_detective(gdb_df: pd.DataFrame, snp_df: pd.DataFrame) -> None:
    # 提取guide+pam
    guide_pam_seq_array = guide_append_pam(gdb_df)
    
    # 多位点indel拆分
    snp_pos_array = multi_indel_split(snp_df)
    
    # 双指针遍历
    
    return


def run_snp(nc_no: str) -> None:
    t1 = time.time()
    # read gdb
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/az_score/{nc_no}.tsv"
    header_type_li = ["string", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category"]
    gdb_df = tsv2df(gdb_path, header_type_li)
    print(f"load gdb<{nc_no}> time cost:{time.time() - t1}")
    t1 = time.time()
    # read snp file
    snp_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/vcf_split/filter/{nc_no}.tsv"
    snp_type_li = ['int32', 'string', 'string', 'string'] 
    snp_df = tsv2df(snp_path, snp_type_li)
    print(f"load snp <{nc_no}> time cost:{time.time() - t1}")
    # process
    t1 = time.time()
    snp_detective(gdb_df, snp_df)
    print(f"process {nc_no} time cost:{time.time() - t1}")
    



def main() -> None:
    run_snp("NC_000001.11")
    return


if __name__ == "__main__":
    main()