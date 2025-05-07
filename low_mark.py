from collections import defaultdict
from pathlib import Path
import time
import pandas as pd
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


def gene_pos_group(start: int, end: int, cut_site: int) -> int:
    # 将基因划分五等分 并根据gRNA切点位置将其归至区间
    gene_len = end - start + 1
    region_len = gene_len / 5
    left_distance = cut_site - start
    return int(left_distance // region_len)


def low_mark(raw_db: str, gene_pos_dict: defaultdict) -> pd.DataFrame:
    # 合并所有子 DataFrame
    sub_dfs = []
    sub20_dfs = []
    mark_dfs = []
    # 读取标记SNP后的原始文件
    df = pd.read_csv(
        raw_db,
        sep="\t",
        header=None,
        dtype={
            24: str
            },
    )
    # 还原首列编号
    df[0] = df[0].apply(lambda x: x.split("_")[1])
    # 按照基因分组排序命名同时标记low score
    for gene, sub_df in df.groupby(by=8, sort=False):
        # 已挑选的id 列表
        picked_id_li = []
        # 添加gRNA名 格式为 hSLC16A1[gRNA3975]
        sub_df = sub_df.reset_index(drop=True)
        sub_df.index += 1
        sub_df[25] = sub_df[8] + "[gRNA" + sub_df.index.astype(str) + "]"
        # 标记 low score 若候选中包含UTR 或 low score(CFD score ≤ 0.1 && RS2 score ≤ 0.3)
        sub_df["_has_0"] = sub_df[24].astype(str).str.contains("0")
        sub_df["_low_quality"] = ((sub_df[18] <= 0.1) | (sub_df[22] <= 0.3)).astype(int)
        # 根据基因起止添加分组信息
        gene_start, gene_end = gene_pos_dict[gene]
        sub_df['group'] = sub_df[7].apply(lambda x:gene_pos_group(gene_start, gene_end, x))
        # 将得分改为百分制并保留两位小数
        sub_df[18] = (sub_df[18] * 100).round(2)
        sub_df[22] = (sub_df[22] * 100).round(2)
        # 将gRNA按照CFD score + RS2 score总分按从高到低排序
        sub_df['sum_score'] = sub_df[18] + sub_df[22]
        # 总分相同时，按照打靶转录本比例从高到低排序
        sub_df['tran_num'] = sub_df[16].apply(
            lambda x: x.split("/")[0]).astype(int)
        # 打靶转录本比例相同时，将CFD score和RS2 score差值从小到大排序
        sub_df['diff_score'] = abs(sub_df[18] - sub_df[22])
        # 按照给定的规则进行排序
        sub_df = sub_df.sort_values(
            by=["_has_0", "_low_quality",'sum_score', 'tran_num', 'diff_score'],
            ascending=[False, True, False, False, True]  # True 表示 False（不包含/不满足）在前
        ).drop(columns=['sum_score', 'tran_num', 'diff_score',"_has_0"]).reset_index(drop=True)
        # 合并所有子 DataFrame
        mark_dfs.append(sub_df)
        # 按照区域分组pick 将各个组别前十进行标记 
        picked_id_li += sub_df[sub_df['group']==0][0][:10].to_list()
        picked_id_li += sub_df[sub_df['group']==1][0][:10].to_list()
        picked_id_li += sub_df[sub_df['group']==2][0][:10].to_list()
        picked_id_li += sub_df[sub_df['group']==3][0][:10].to_list()
        picked_id_li += sub_df[sub_df['group']==4][0][:10].to_list()
        sub_df['picked'] = sub_df[0].isin(picked_id_li)
        # 若挑选的数目为50 直接筛出 否则回补
        sub_df = pd.concat([sub_df[sub_df['picked']==True],sub_df[sub_df['picked']==False].head(50-len(picked_id_li))])
        # 删除临时列
        sub_df = sub_df.drop(columns=['group', 'picked'])
        # 选取前二十个
        sub20_df = sub_df.head(20)
        # 合并所有子 DataFrame
        sub_dfs.append(sub_df)
        sub20_dfs.append(sub20_df)
    new_df = pd.concat(sub_dfs, ignore_index=True)
    new20_df = pd.concat(sub20_dfs, ignore_index=True)
    mark_df = pd.concat(mark_dfs, ignore_index=True)
    # 调整列顺序
    new_header = [25] + list(range(21)) + [22, 21, "_low_quality", 23, 24]
    return mark_df[new_header], new_df[new_header], new20_df[new_header]


def run_mark(nc_no: str) -> None:
    t1 = time.time()
    # print(f"{nc_no} start !!!")
    raw_db = f"/mnt/ntc_data/wayne/Repositories/CRISPR/utr_mark/{nc_no}.tsv"
    mark_output = f"/mnt/ntc_data/wayne/Repositories/CRISPR/low_mark/{nc_no}.tsv"
    filter_output = f"/mnt/ntc_data/wayne/Repositories/CRISPR/filter_50/{nc_no}.tsv"
    filter20_output = f"/mnt/ntc_data/wayne/Repositories/CRISPR/filter_20/{nc_no}.tsv"
    # if Path(filter_output).exists():
    #     print(f"{nc_no} exists !!!")
    #     return
    # 读取基因位置信息并存为字典
    gene_info_dict = defaultdict(list)
    gene_info_df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/Gene_list.tsv",sep="\t",header=None)
    for gene, start, end in zip(gene_info_df[0], gene_info_df[1], gene_info_df[2]):
        gene_info_dict[gene] = [start, end]
    mark_df, new_df, new20_df = low_mark(raw_db,gene_info_dict)
    # 保存为tsv
    mark_df.to_csv(mark_output, sep="\t", header=None, index=None)
    new_df.to_csv(filter_output, sep="\t", header=None, index=None)
    new20_df.to_csv(filter20_output, sep="\t", header=None, index=None)
    print(f"{nc_no} finished,time cost:{time.time() - t1}!!!")
    

def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_mark, nc_li, 24)
    # run_mark(nc_li[-1])


if __name__ == "__main__":
    main()
