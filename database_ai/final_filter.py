from pathlib import Path
from collections import defaultdict
import time
import pandas as pd
from sys import path
<<<<<<< HEAD
path.append("/mnt_data/Wayne/Repositories/CRISPR/")
=======
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
>>>>>>> fce1903c581ad706ef1336b7872026b2ae6b2929
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure

def gene_pos_group(gene_pos_str: str) -> int:
    # 根据gRNA端点距离较小的位置 划分为三个等级：<=100bp;101-300bp;>300 bp，
    gene_pos_li = gene_pos_str.replace("(","").replace(")","").split("-")
    gene_pos_li = [int(x) for x in gene_pos_li if x]
    distance = min(gene_pos_li)
    if distance <= 100:return 1
    if distance <= 300:return 2
    return 3



def low_mark_group(raw_db: str) -> pd.DataFrame:
    # 合并所有子 DataFrame
    sub_dfs = []
<<<<<<< HEAD
    # 读取re_num后的原始文件
=======
    # 读取标记SNP后的原始文件
>>>>>>> fce1903c581ad706ef1336b7872026b2ae6b2929
    df = pd.read_csv(
        raw_db,
        sep="\t",
        header=None,
<<<<<<< HEAD
        
    )
    # 还原首列编号
    # df[0] = df[0].apply(lambda x: x.split("_")[1])
    # 按照基因分组排序命名同时标记low score
    for gene, sub_df in df.groupby(by=8, sort=False):
        # 已挑选的id 列表
        picked_id_li = []
        # 添加gRNA名 格式为 hSLC16A1[gRNA3975]
        # sub_df = sub_df.reset_index(drop=True)
        # sub_df.index += 1
        # sub_df["new_name"] = sub_df[7] + "[gRNA" + sub_df.index.astype(str) + "]"
        # 根据基因起止添加分组信息
        sub_df['group'] = sub_df[12].apply(gene_pos_group)
        # 将得分改为百分制并保留两位小数
        sub_df[14] = (sub_df[14] * 100).round(2)
        sub_df[18] = (sub_df[18] * 100).round(2)
        # 将gRNA按照组别和总分从高到低排序
        sub_df['sum_score'] = sub_df[14] + sub_df[18]
=======
        dtype={
            24: str
            },
    )
    # 还原首列编号
    df[0] = df[0].apply(lambda x: x.split("_")[1])
    # 按照基因分组排序命名同时标记low score
    for gene, sub_df in df.groupby(by=7, sort=False):
        # 已挑选的id 列表
        picked_id_li = []
        # 添加gRNA名 格式为 hSLC16A1[gRNA3975]
        sub_df = sub_df.reset_index(drop=True)
        sub_df.index += 1
        sub_df["new_name"] = sub_df[7] + "[gRNA" + sub_df.index.astype(str) + "]"
        # 根据基因起止添加分组信息
        sub_df['group'] = sub_df[11].apply(gene_pos_group)
        # 将得分改为百分制并保留两位小数
        sub_df[13] = (sub_df[13] * 100).round(2)
        sub_df[17] = (sub_df[17] * 100).round(2)
        # 将gRNA按照组别和总分从高到低排序
        sub_df['sum_score'] = sub_df[13] + sub_df[17]
>>>>>>> fce1903c581ad706ef1336b7872026b2ae6b2929
        sub_df = sub_df.sort_values(
            by=['group', 'sum_score'],  # 排序优先级：和 -> tran_num -> 差值
            ascending=[True, False]       # 和降序，tran_num降序，差值升序
        )
        # 按照区域分组pick 将各个组别前十进行标记 
        picked_id_li += sub_df[sub_df['group']==1][0][:20].to_list()
        picked_id_li += sub_df[sub_df['group']==2][0][:20].to_list()
        picked_id_li += sub_df[sub_df['group']==3][0][:10].to_list()
        sub_df['picked'] = sub_df[0].isin(picked_id_li)
        # 若挑选的数目为50 直接筛出 否则回补
        sub_df = pd.concat([sub_df[sub_df['picked']==True],sub_df[sub_df['picked']==False].head(50-len(picked_id_li))])
        # 删除临时列
        sub_df = sub_df.drop(columns=['sum_score', 'group', 'picked'])
        # 合并所有子 DataFrame
        sub_dfs.append(sub_df)
    new_df = pd.concat(sub_dfs, ignore_index=True)
    # 调整列顺序
<<<<<<< HEAD
    # new_header = ["new_name"] + list(range(19))
    # return new_df[new_header]
    return new_df
=======
    new_header = ["new_name"] + list(range(19))
    return new_df[new_header]
>>>>>>> fce1903c581ad706ef1336b7872026b2ae6b2929


def run_mark(nc_no: str) -> None:
    t1 = time.time()
    # print(f"{nc_no} start !!!")
<<<<<<< HEAD
    raw_db = f"/mnt_data/Wayne/Repositories/CRISPR/database_ai/re_num/{nc_no}.tsv"
    filter_output = f"/mnt_data/Wayne/Repositories/CRISPR/database_ai/final_filter/{nc_no}.tsv"
    # if Path(filter_output).exists():
    #     print(f"{nc_no} exists !!!")
    #     return
=======
    raw_db = f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/snp_mark/{nc_no}.tsv"
    filter_output = f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/final_filter/{nc_no}.tsv"
    if Path(filter_output).exists():
        print(f"{nc_no} exists !!!")
        return
>>>>>>> fce1903c581ad706ef1336b7872026b2ae6b2929
    new_df = low_mark_group(raw_db)
    # 保存为tsv
    new_df.to_csv(filter_output, sep="\t", header=None, index=None)
    print(f"{nc_no} finished,time cost:{time.time() - t1},filter_50 num {len(new_df)}!!!")
    

def main() -> None:
<<<<<<< HEAD
    nc2chr_file = "/mnt_data/Wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_mark, nc_li, 48)
    # run_mark(nc_li[11])
=======
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_mark, nc_li, 24)
    # run_mark(nc_li[-1])
>>>>>>> fce1903c581ad706ef1336b7872026b2ae6b2929


if __name__ == "__main__":
    main()
