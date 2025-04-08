from pathlib import Path
import time
import pandas as pd
from collections import defaultdict
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from cds_mark import gene_ori_dict
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


def merge_intervals(intervals: list) -> list:
    if not intervals:
        return []

    # 按起始位置排序
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]

    for current in intervals[1:]:
        prev_start, prev_end = merged[-1]
        curr_start, curr_end = current

        if curr_start <= prev_end:  # 存在重叠
            merged[-1] = (prev_start, max(prev_end, curr_end))  # 合并
        else:
            merged.append(current)

    return merged


# 根据exon（含cds和UTR）和cds 计算5'UTR和3'UTR位置
def utr_region_obtain(exon_file_path: str, cds_file_path: str, ori_dict: dict) -> defaultdict:
    # 读取exon和cds region 文件
    exon_df = tsv2df(exon_file_path,[])
    cds_df = tsv2df(cds_file_path,[])  
    exon_region_dict = defaultdict(lambda:defaultdict(list))
    interval_exon_region_dict = defaultdict(lambda:defaultdict(list))
    # 按照基因分组
    for gene, sub_exon_df in exon_df.groupby(0,sort=False):
        # 若exon无cds（NR） 跳过
        sub_cds_df = cds_df[cds_df[0]==gene]
        if sub_cds_df.empty:continue
        gene_ori = ori_dict[gene]
        # 根据转录本分组 
        # for (exon_tran, sub_exon_tran_df),(cds_tran, sub_cds_tran_df) in zip(sub_exon_df.groupby(1,sort=False), sub_cds_df.groupby(1,sort=False)):
        for exon_tran, sub_exon_tran_df in sub_exon_df.groupby(1,sort=False):
            sub_cds_tran_df = sub_cds_df[sub_cds_df[1]==exon_tran]
            # 若exon对应转录本无cds（NR） 跳过
            if sub_cds_tran_df.empty:continue
            # 记录所有cds区域
            cds_start_li = sub_cds_tran_df[2].to_numpy()
            cds_end_li = sub_cds_tran_df[3].to_numpy()
            exon_start_li = sub_exon_tran_df[2].to_numpy()
            exon_end_li = sub_exon_tran_df[3].to_numpy()
            [exon_region_dict[gene]['CDS'].append((int(cds_start_li[i]), int(cds_end_li[i]))) for i in range(len(cds_end_li))]
            # 根据方向分别寻找UTR起始终止             
            utr5_start = int(exon_start_li[0])
            utr5_end = int(cds_start_li[0]) - 1
            utr3_start = int(cds_end_li[-1]) + 1
            utr3_end = int(exon_end_li[-1])
            
            if gene_ori == "-":
                utr5_start = int(cds_end_li[0]) + 1
                utr5_end = int(exon_end_li[0])
                utr3_start = int(exon_start_li[-1])
                utr3_end = int(cds_start_li[-1]) - 1
                
            # 将utr region 存入字典
            exon_region_dict[gene]["UTR5"].append((utr5_start, utr5_end))
            exon_region_dict[gene]["UTR3"].append((utr3_start, utr3_end))          
    
    # 分别对UTR5 UTR3 CDS 求区域并集 
    for gene in exon_region_dict.keys():
        interval_exon_region_dict[gene]["CDS"] = merge_intervals(exon_region_dict[gene]["CDS"])
        interval_exon_region_dict[gene]["UTR5"] = merge_intervals(exon_region_dict[gene]["UTR5"])
        interval_exon_region_dict[gene]["UTR3"] = merge_intervals(exon_region_dict[gene]["UTR3"])  
    # print(interval_exon_region_dict)
    return interval_exon_region_dict


def region_code2str(region_code: int) -> str:
    '''
    1->utr5,2->cds,4->utr3
    -1->utr5,0->cds,1->utr3
    '''
    bin_region_code = '{:03b}'.format(region_code)[::-1]
    mark_li = ["-1", "0", "1"]
    res_li = [mark_li[i] for i in range(3) if int(bin_region_code[i])]
    return ",".join(res_li)


def region_classify(gene_name: str, cut_pos: int, grna_ori: str, utr_pos_dict: defaultdict) -> str:
    '''
    1->utr5,2->cds,4->utr3
    -1->utr5,0->cds,1->utr3
    '''
    region_code = 0
    cut_pos_offset = cut_pos + float(grna_ori + "0.5")
    cds_region_li = utr_pos_dict[gene_name]["CDS"]
    utr5_region_li = utr_pos_dict[gene_name]["UTR5"]
    utr3_region_li = utr_pos_dict[gene_name]["UTR3"]
    
    for start, end in cds_region_li:
        if start < cut_pos_offset < end:
            region_code += 2
            break
    
    for start, end in utr5_region_li:
        if start < cut_pos_offset < end:
            region_code += 1
            break
    
    for start, end in utr3_region_li:
        if start < cut_pos_offset < end:
            region_code += 4
            break
    return region_code2str(region_code)


def search_regions(grna_table: str, region_dict: defaultdict) -> pd.DataFrame:
    gdb_df = tsv2df(grna_table,[])
    gdb_df["Target Region"] = gdb_df.apply(lambda x:region_classify(x[9],x[8],x[5],region_dict),axis=1)
    return gdb_df
    



# 根据UTR位置 添加新列 region 标记UTR和CDS
def utr_mark(nc_no: str) -> None:
    t1 = time.time()
    print(f"{nc_no} start:{t1}")
    exon_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/EXON.tsv"
    cds_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/CDS.tsv"
    gdb_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/low_mark/{nc_no}.tsv"
    res_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/utr_mark/{nc_no}.tsv"
    if Path(res_file).exists():
        print(f"{nc_no} is exist,exit !!!")
        return
    ori_dict = gene_ori_dict(nc_no)
    utr_pos_dict = utr_region_obtain(exon_file, cds_file, ori_dict)
    res_df = search_regions(gdb_file, utr_pos_dict)
    res_df.to_csv(res_file,sep="\t",header=False,index=False)
    print(f"{nc_no} annotation time cost:{time.time() - t1}")

def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(utr_mark,nc_li,24)
    # utr_mark(nc_no)


if __name__ == "__main__":
    main()
