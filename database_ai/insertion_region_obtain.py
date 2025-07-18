import pandas as pd
from collections import defaultdict

from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from generate_split_ori import async_in_iterable_structure

# 根据基因方向寻找TSS上游1kb 并获取交集区域 记录落在交集区能影响的转录本

# 自定义：合并区间
def merge_intervals(intervals):
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return ",".join(["-".join([str(x) for x in interval]) for interval in merged])



def found_insertion(nc_no: str):
    # 读取基因的tss tran 结果
    df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/tss_tran/{nc_no}.tsv",sep="\t",header=None)
    # 根据偏移量获取3'位置
    df["offset"] = (df[1].astype(str) + '1').astype(int)
    df["start"] = df[2] - df["offset"] * 1000
    df["end"] = df[2] - df["offset"]
    # 生成区间列
    df["region_start"] = df[["start", "end"]].min(axis=1)
    df["region_end"] = df[["start", "end"]].max(axis=1)
    # 以基因为单位合并区间并映射至新列
    result_dict = (
        df.groupby(0,sort=False)
        .apply(lambda g: merge_intervals(list(zip(g["region_start"], g["region_end"]))))
        .to_dict()
    )
    df["regions"] = df[0].map(result_dict)
    # 去除重复行和多余列
    df = df[[0,1,"regions"]].drop_duplicates()
    # 写入文件
    df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/tss_regions/{nc_no}.tsv",sep="\t",header=None,index=False)
        


def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(found_insertion,nc_li,24)
    # found_insertion(nc_li[-1])
    

if __name__ == "__main__":
    main()