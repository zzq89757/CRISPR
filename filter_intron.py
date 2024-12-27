import time
import pandas as pd
import psutil
import os
from sys import path

path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")

from process_border_withid import scaffold_detective_numpy
from generate_split_ori import async_in_iterable_structure


def gdb2df(gdb_path: str) -> pd.DataFrame:
    type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32"]

    type_dict = dict(enumerate(type_li))

    gdb_df = pd.read_csv(
        gdb_path,
        sep="\t",
        header=None,
        dtype=type_dict,
        # low_memory=False,
    )
    # print(gdb_df)
    return gdb_df


def region2df(region_file: str) -> pd.DataFrame:
    type_li = ["string", "string", "int32", "int32"]

    type_dict = dict(enumerate(type_li))

    region_df = pd.read_csv(
        region_file,
        sep="\t",
        header=None,
        dtype=type_dict,
        # low_memory=False,
    )
    # print(region_df)
    return region_df


def merge_tran_exon(df: pd.DataFrame, nc_no: str) -> pd.DataFrame:
    df.columns = ["Gene", "Transcript", "Start", "End", "Exon"]
    # 为每个 `[Start, End]` 分配唯一组号
    # print(df[df["Gene"]=="GTPBP6"])
    df['Group'] = df.groupby(['Gene', 'Start', 'End']).ngroup() + 1
    # 转录本和 Exon 合并
    df['Transcript_Exon'] = df['Transcript'] + '#' + df['Exon'].astype(str)
    # df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/exon_filter/{nc_no}_region.tsv",header=None,index=False,sep="\t")

    # 按组聚合数据
    grouped_df = df.groupby('Group').agg({
        'Gene': 'first',  # 每组基因取第一个
        'Start': 'first',  # 每组起始位置取第一个
        'End': 'first',    # 每组结束位置取第一个
        'Transcript_Exon': lambda x: ';'.join(x)  # 转录本#外显子编号以分号分隔
    }).reset_index(drop=True)
    # grouped_df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/exon_filter/{nc_no}_region.tsv",header=None,index=False,sep="\t")
    grouped_df = grouped_df.sort_values("Start")
    grouped_df.columns = [0,1,2,3]
    # 查看结果
    # print(grouped_df)
    return grouped_df


def filter_intron(nc_no :str, gdb_df: pd.DataFrame, exon_df: pd.DataFrame, output_file: str) -> None:
    # output_file = "exon_filter/exu.tsv"
    output_handle = open(output_file,'w')
    
    # 根据gdb 的基因和切点 找到cut的转录本、外显子以及cds并集前2/3位置
    # group_order = pd.Categorical(exon_df[1], categories=exon_df[1].unique(), ordered=True)
    # exon_df = exon_df.sort_values(by=[2, 3], key=lambda col: group_order if col.name == 2 else col)
    exon_df[4] = exon_df.groupby(1).cumcount() + 1
    # exon_df = exon_df.sort_values(2)
    exon_df = merge_tran_exon(exon_df, nc_no)
    scaffold_pos_li = scaffold_detective_numpy(exon_df)
    # print(scaffold_pos_li)
    # exon_df.to_csv("exp.tsv",sep="\t",header=None,index=False)
    exon_array = exon_df.to_numpy()
    # 提取所需列并矢量化计算
    gdb_ori_array = gdb_df[4].astype(str) + "0.5"
    gdb_cut_pos_array = gdb_df[7].to_numpy()
    gdb_gene_name_array = gdb_df[8].to_numpy()
    gdb_array = gdb_df.to_numpy()
    
    gene_start_array = exon_df[1].to_numpy()
    gene_end_array = exon_df[2].to_numpy()

    # 计算 g_pos 的矢量化结果
    g_pos_array = gdb_cut_pos_array + gdb_ori_array.astype(float)

    # 筛选 exon 数据
    result = []
    res_count = 0
    idx = 0
    current_gene_pos_idx = 0
    scaffold_pos_len = len(scaffold_pos_li)
    gene_pos_len = len(gene_start_array)
    current_scaffold_pos_idx = 0
    for g_idx, g_pos in enumerate(g_pos_array):
        # idx += 1
        # # 一次性筛选出满足条件的 exon 数据
        # sub_exon_df = exon_df[
        #     (exon_df[0] == gene_name) & (g_pos > exon_df[2]) & (g_pos < exon_df[3])
        # ]
        # print(f"{idx}/{len(gdb_gene_name_array)} finished !!!",end="\r",flush=True)
        # if sub_exon_df.empty:
        #     res_count += 1
        gene_start = gene_start_array[current_gene_pos_idx]
        gene_end = gene_end_array[current_gene_pos_idx]
        if g_pos < gene_start:
            continue
        # 跨多个基因时的情况
        while (current_gene_pos_idx < gene_pos_len -1  and g_pos > gene_end):
            current_gene_pos_idx += 1
            # 同步更新scaffold_pos_idx
            border_right = scaffold_pos_li[current_scaffold_pos_idx][1] if current_scaffold_pos_idx < scaffold_pos_len else scaffold_pos_li[-1][1]
            if current_gene_pos_idx > border_right and current_scaffold_pos_idx <= scaffold_pos_len:
                current_scaffold_pos_idx += 1
            gene_start = gene_start_array[current_gene_pos_idx]
            gene_end = gene_end_array[current_gene_pos_idx]
        
        # 大于max end时 结束
        if current_gene_pos_idx == gene_pos_len -1 and g_pos > gene_end_array[-1]:
            break
        
        # 处于scaffold中的gRNA处理
        if current_scaffold_pos_idx < scaffold_pos_len and scaffold_pos_li[current_scaffold_pos_idx][0] <= current_gene_pos_idx <= scaffold_pos_li[current_scaffold_pos_idx][1]:
            tran_li = []
            for j in range(current_gene_pos_idx, scaffold_pos_li[current_scaffold_pos_idx][1] + 1):
                # if gene_start_array[j] < g_pos < gene_end_array[j]:
                if gene_start_array[j] < g_pos < gene_end_array[j]:
                    if gdb_array[g_idx][8] == exon_array[j][0]:
                        # print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                        # print(exon_array[j][3],end="\n",file=output_handle)
                        tran_li.append(exon_array[j][3])
                    
                    # print(g_pos,end="\t")
                        ...
                    # print(gene_start_array[j],end="\t")
                    # print(gene_end_array[j],end="\t")
                    # print(gene_pos_array[j])
            if tran_li:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                print(";".join(tran_li),end="\n",file=output_handle)
            continue
        
        
        # 如果 gRNA 位于当前基因范围内，记录信息 chr1:1335323 跨两个或以上的位点会被跳过
        if gene_start < g_pos < gene_end:
            # print(f"{g_pos} in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}", flush=True, end="\r")
            # print(f"{g_pos} in {gene_start}-{gene_end}")
            # g_item = gdb.iloc[g_idx]
            # g_item = gdb_array[g_idx]
            if gdb_array[g_idx][8] == exon_array[current_gene_pos_idx][0]:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                print(exon_array[current_gene_pos_idx][3],end="\n",file=output_handle)
            #     print(gdb_array[g_idx][8],end="\t")
            #     print(exon_array[current_gene_pos_idx][0],end="\t")
            #     print(g_pos)
            # ...
            # print(g_item)
        else:
            # print(f"{g_pos} not in {gene_start}-{gene_end} and {gene_start_array[current_gene_pos_idx - 1]}-{gene_end_array[current_gene_pos_idx - 1]}")
            # g_item = gdb_array[g_idx]
            # print(g_item)
            ...
        
        
    return 


def run_filter(nc_no: str) -> None:
    t1 = time.time()
    process = psutil.Process(os.getpid())
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc2chr_dict = dict(zip(nc_df[0],nc_df[1]))
    
    # 读取gdb
    gdb_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/gene_annotation/spCas9_Homo_chrY.tsv"
    gdb_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/gene_annotation/spCas9_Homo_chr1.tsv"
    gdb_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/gene_annotation/spCas9_Homo_{nc2chr_dict[nc_no]}.tsv"
    gdb_df = gdb2df(gdb_file)
    # gdb_df = pd.DataFrame([])
    # 读取对应的exon表
    exon_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/NC_000024.10/EXON.tsv"
    exon_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/NC_000001.11/EXON.tsv"
    exon_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/EXON.tsv"
    output_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/exon_filter/{nc_no}.tsv"
    exon_df = region2df(exon_file)
    print(f"load gdb time cost:{time.time() - t1}")
    # 注释gdb
    t1 = time.time()
    filter_intron(nc_no, gdb_df, exon_df, output_file)
    memory_info = process.memory_info()
    peak_memory_gb = memory_info.peak_wset / (1024**3) if hasattr(memory_info, 'peak_wset') else memory_info.rss / (1024**3)

    print(f"\n峰值内存使用: {peak_memory_gb:.2f} GB")
    print(f"annotation time cost:{time.time() - t1}")
    


def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_filter,nc_li,8)
    # run_filter(nc_li[23])

if __name__ == "__main__":
    main()