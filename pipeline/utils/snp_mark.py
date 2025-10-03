import pandas as pd
import numpy as np
from .filter_intron import gdb2df
# 用snp文件遍历gdb起止 似乎不可行 gdb交集区太多


def append_middle_pos(gdb_df: pd.DataFrame) -> pd.DataFrame:
    # 将guide + Pam 23bp的中点 用于snp覆盖判断
    offset_array = gdb_df[4].astype(str) + "3"
    start_array = gdb_df[5].to_numpy()
    end_array = gdb_df[6].to_numpy() + offset_array.astype(int).to_numpy()
    gdb_df['middle'] = (start_array + end_array)/2
    gdb_df = gdb_df.sort_values("middle")
    gdb_df['middle'] = gdb_df['middle'].where(~gdb_df['middle'].duplicated(keep='first'), 0)
    gdb_df['middle'] = gdb_df['middle'].astype(int)
    return gdb_df

def multi_indel_split(snp_df: pd.DataFrame) -> np.ndarray:
    # 将多位点indel拆分成多个item 并返回pos array
    pos_array = snp_df[1].to_numpy()
    ref_array = snp_df[3].to_numpy()
    alter_array = snp_df[4].to_numpy()
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


def double_pointer(gdb_df: pd.DataFrame, snp_pos_array: np.ndarray) -> pd.DataFrame:
    middle_pos_array = gdb_df['middle'].to_numpy()
    range_val = 11
    result = np.full(len(middle_pos_array), -1)  # 初始化结果数组为 -1
    pointer2 = 0
    dup_idx_li = []
    for i, val1 in enumerate(middle_pos_array):
        # 如果 middle_pos_array 中的值为 0，跳过判断
        if int(val1) == 0:
            dup_idx_li.append(i)
            continue

        # 确保 pointer2 指向 snp_pos_array 中第一个大于等于 val1 - 11 的位置
        while pointer2 < len(snp_pos_array) and snp_pos_array[pointer2] < val1 - range_val:
            pointer2 += 1

        # 检查 snp_pos_array 中的数字是否在 [val1-11, val1+11] 范围内
        covered = False
        while pointer2 < len(snp_pos_array) and snp_pos_array[pointer2] <= val1 + range_val:
            covered = True
            break

        # 根据覆盖情况设置结果
        result[i] = 1 if covered else 0
    # 填充重复值
    for i in dup_idx_li:
        result[i] = result[i - 1]
    gdb_df['snp'] = result
    del(gdb_df["middle"])
    gdb_df = gdb_df.sort_index()
    return gdb_df

# 双指针遍历 snp和gdb
def snp_detective(gdb_df: pd.DataFrame, snp_df: pd.DataFrame) -> pd.DataFrame:
    # 提取guide+pam
    gdb_df = append_middle_pos(gdb_df)
    
    # 多位点indel拆分
    snp_pos_array = multi_indel_split(snp_df)
    # 双指针遍历
    gdb_df = double_pointer(gdb_df, snp_pos_array)  
    
    return gdb_df


def snp_mark(project_dir: str, nc_no: str) -> None:
    # read gdb
    gdb_path = f"{project_dir}/rs2_score/{nc_no}.tsv"
    gdb_df = gdb2df(gdb_path)
    # read snp file
    snp_path = f"{project_dir}/GCF/vcf/{nc_no}.vcf"
    snp_df = pd.read_csv(snp_path,sep="\t",header=None,usecols=[1,3,4])
    # snp detective
    gdb_df = snp_detective(gdb_df, snp_df)
    # 调整表头顺序
    # headers = list(range(17)) + ["snp",18,19,] 
    # 保存为tsv文件
    output_path = f"{project_dir}/snp_mark/{nc_no}.tsv"
    gdb_df.to_csv(output_path,sep="\t",header=None,index=None)
    



def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # async_in_iterable_structure(run_snp,nc_li,24)
    # run_snp("NC_000024.10")
    return


if __name__ == "__main__":
    main()