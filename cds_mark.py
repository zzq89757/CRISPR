from multiprocessing import process
import os
import time
import pandas as pd
from sys import path
import psutil
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")




from process_border_withid import scaffold_detective_numpy
from generate_split_ori import async_in_iterable_structure
from filter_intron import  region2df



def gdb2df(gdb_path: str) -> pd.DataFrame:
    type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string"]

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




def merge_tran_exon(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = ["Gene", "Transcript", "Start", "End", "Exon"]
    # 为每个 `[Start, End]` 分配唯一组号
    df['Group'] = df.groupby(['Gene', 'Start', 'End']).ngroup() + 1

    # 转录本和 Exon 合并
    df['Transcript_Exon'] = df['Transcript'] + '#' + df['Exon'].astype(str)

    # 按组聚合数据
    grouped_df = df.groupby('Group').agg({
        'Gene': 'first',  # 每组基因取第一个
        'Start': 'first',  # 每组起始位置取第一个
        'End': 'first',    # 每组结束位置取第一个
        'Transcript_Exon': lambda x: ';'.join(x)  # 转录本#外显子编号以分号分隔
    }).reset_index(drop=True)
    grouped_df = grouped_df.sort_values("Start")
    grouped_df.columns = [0,1,2,3]
    # 查看结果
    # print(grouped_df)
    return grouped_df


def common_cds_front_region(cds_df: pd.DataFrame) -> pd.DataFrame:
    new_res = pd.DataFrame([])
    cds_df.columns = ["Gene", "Transcript", "CDS_start", "CDS_end"]
    # 按照转录本分组
    for tran, tran_df in cds_df.groupby("Transcript"):
        result = []
        # 计算总长
        full_cds_len = sum(tran_df["CDS_end"] - tran_df["CDS_start"])
        # 找出前2/3 full_cds_len - 每行的区间长度 当剩余数字为负时 将终点改为该行终点加上剩余数字 丢弃后续cds区间
        target_length = full_cds_len * 2 / 3  # 目标长度（前2/3）

        accumulated_length = 0  # 累积长度
        for idx, row in tran_df.iterrows():
            start, end = row["CDS_start"], row["CDS_end"]
            interval_length = end - start

            # 判断是否超过目标长度
            if accumulated_length + interval_length > target_length:
                # 修正终点，取部分区间
                remaining_length = target_length - accumulated_length
                new_end = start + remaining_length
                result.append([row["Gene"], tran, start, int(new_end)])
                break
            else:
                # 累加区间
                result.append([row["Gene"], tran, start, end])
                accumulated_length += interval_length

        # 转为 DataFrame
        result_df = pd.DataFrame(result, columns=["Gene", "Transcript", "CDS_start", "CDS_end"])
        # 查看结果
        new_res = pd.concat([new_res, result_df])
    new_res = new_res.sort_values("CDS_start")
    # new_res.to_csv("cdt.tsv",header=None,index=False,sep="\t")
    # print(new_res)
    return new_res

    
def mark_cds(gdb_df: pd.DataFrame, cds_df: pd.DataFrame, output_file: str) -> None:
    # output_file = "exon_filter/exu.tsv"
    output_handle = open(output_file,'w')
    
    # 根据gdb 的基因和切点 找到cut的转录本、外显子以及cds并集前2/3位置
    # group_order = pd.Categorical(cds_df[1], categories=cds_df[1].unique(), ordered=True)
    # cds_df = cds_df.sort_values(by=[2, 3], key=lambda col: group_order if col.name == 2 else col)
    # cds_df[4] = cds_df.groupby(1).cumcount() + 1
    # cds_df = cds_df.sort_values(2)
    cds_df = common_cds_front_region(cds_df)
    cds_df.columns = [0,4,1,2]
    scaffold_pos_li = scaffold_detective_numpy(cds_df)
    # print(scaffold_pos_li)
    # cds_df.to_csv("cdsp.tsv",sep="\t",header=None,index=False)
    exon_array = cds_df.to_numpy()
    # 提取所需列并矢量化计算
    gdb_ori_array = gdb_df[4].astype(str) + "0.5"
    gdb_cut_pos_array = gdb_df[7].to_numpy()
    gdb_gene_type_array = gdb_df[10].to_numpy()
    gdb_array = gdb_df.to_numpy()
    
    gene_start_array = cds_df[1].to_numpy()
    gene_end_array = cds_df[2].to_numpy()

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
        # sub_cds_df = cds_df[
        #     (cds_df[0] == gene_name) & (g_pos > cds_df[2]) & (g_pos < cds_df[3])
        # ]
        # print(f"{idx}/{len(gdb_gene_name_array)} finished !!!",end="\r",flush=True)
        # if sub_cds_df.empty:
        #     res_count += 1
        gene_start = gene_start_array[current_gene_pos_idx]
        gene_end = gene_end_array[current_gene_pos_idx]
        if  gdb_gene_type_array[g_idx] == "non_coding" or g_pos < gene_start or g_pos > gene_end_array[-1]:
            print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
            print("no",end="\n",file=output_handle)
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
        # if current_gene_pos_idx == gene_pos_len -1 and g_pos > gene_end_array[-1]:
        #     break
        
        # 处于scaffold中的gRNA处理
        if current_scaffold_pos_idx < scaffold_pos_len and scaffold_pos_li[current_scaffold_pos_idx][0] <= current_gene_pos_idx <= scaffold_pos_li[current_scaffold_pos_idx][1]:
            fall_flag = 0
            for j in range(current_gene_pos_idx, scaffold_pos_li[current_scaffold_pos_idx][1] + 1):
                # if gene_start_array[j] < g_pos < gene_end_array[j]:
                if gene_start_array[j] < g_pos < gene_end_array[j]:
                    if gdb_array[g_idx][8] == exon_array[j][0]:
                        fall_flag = 1
                        print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                        # print(exon_array[current_gene_pos_idx][2],end="\t",file=output_handle)
                        # print(exon_array[current_gene_pos_idx][3],end="\n",file=output_handle)
                        print("yes",end="\n",file=output_handle)
                    break
            
            if fall_flag == 0:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                print("no",end="\n",file=output_handle)
            continue
        
        
        # 如果 gRNA 位于当前基因范围内，记录信息 chr1:1335323 跨两个或以上的位点会被跳过
        if gene_start < g_pos < gene_end:
            
            # print(f"{g_pos} in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}", flush=True, end="\r")
            # print(f"{g_pos} in {gene_start}-{gene_end}")
            # g_item = gdb.iloc[g_idx]
            # g_item = gdb_array[g_idx]
            if gdb_array[g_idx][8] == exon_array[current_gene_pos_idx][0]:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                # print(exon_array[current_gene_pos_idx][2],end="\t",file=output_handle)
                # print(exon_array[current_gene_pos_idx][3],end="\n",file=output_handle)
                print("yes",end="\n",file=output_handle)
            else:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                print("no",end="\n",file=output_handle)
            #     print(gdb_array[g_idx][8],end="\t")
            #     print(exon_array[current_gene_pos_idx][0],end="\t")
            #     print(g_pos)
            # ...
            # print(g_item)
        else:
            print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
            print("no",end="\n",file=output_handle)
            # print(f"{g_pos} not in {gene_start}-{gene_end} and {gene_start_array[current_gene_pos_idx - 1]}-{gene_end_array[current_gene_pos_idx - 1]}")
            # g_item = gdb_array[g_idx]
            # print(g_item)
            ...
        
        
    return 

    
    
def run_mark(nc_no) -> None:
    t1 = time.time()
    print(f"{nc_no} mark cds start !!!")
    process = psutil.Process(os.getpid())
    gdb_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/exon_filter/NC_000024.10.tsv"
    gdb_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/exon_filter/{nc_no}.tsv"
    gdb_df = gdb2df(gdb_file)
    # 读取对应的CDS exon 表
    cds_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/NC_000024.10/CDS.tsv"
    cds_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/CDS.tsv"
    cds_df = region2df(cds_file)
    # cds_df.to_csv("./nc22.cds",header=None,sep="\t")
    # mark_cds(gdb_df,cds_df,"/mnt/ntc_data/wayne/Repositories/CRISPR/cds_mark/NC_000024.10.tsv")
    mark_cds(gdb_df,cds_df,f"/mnt/ntc_data/wayne/Repositories/CRISPR/cds_mark/{nc_no}.tsv")
    memory_info = process.memory_info()
    peak_memory_gb = memory_info.peak_wset / (1024**3) if hasattr(memory_info, 'peak_wset') else memory_info.rss / (1024**3)

    print(f"峰值内存使用: {peak_memory_gb:.2f} GB")
    print(f"{nc_no} mark cds time cost:{time.time() - t1}")



def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_mark,nc_li,24)
    # run_mark(nc_li[21])
    
    

if __name__ == "__main__":
    main()