from collections import defaultdict
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from cds_mark import gene_ori_dict
from utils.read_tsv import tsv2df


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
            [exon_region_dict[gene]['CDS'].append((int(cds_start_li[i]), int(cds_end_li[i]))) for i in range(len(cds_end_li))]
            # 根据方向分别寻找UTR起始终止
            cds_start_idx = 0 if gene_ori == "+" else -1
            cds_end_idx = -1 - cds_start_idx
            cds_start = cds_start_li[cds_start_idx]
            cds_end = cds_end_li[cds_end_idx]
            exon_start = sub_exon_tran_df[2].to_numpy()[cds_start_idx]
            exon_end = sub_exon_tran_df[3].to_numpy()[cds_end_idx]
            utr5_start = int(exon_start)
            utr5_end = int(cds_start) - 1
            utr3_start = int(cds_end) + 1
            utr3_end = int(exon_end)
            # 将utr region 存入字典
            exon_region_dict[gene]["UTR5"].append((utr5_start, utr5_end))
            exon_region_dict[gene]["UTR3"].append((utr3_start, utr3_end))          
    
    # 分别对UTR5 UTR3 CDS 求区域并集 
    for gene in exon_region_dict.keys():
        exon_region_dict[gene]["CDS"] = merge_intervals(exon_region_dict[gene]["CDS"])
        exon_region_dict[gene]["UTR5"] = merge_intervals(exon_region_dict[gene]["UTR5"])
        exon_region_dict[gene]["UTR3"] = merge_intervals(exon_region_dict[gene]["UTR3"])

    # 拆分为 -1->5UTR,0->CDS,1->3UTR (真utr : utr 并集 - cds并集) 落在多个区域时用逗号分隔
    



# 根据UTR位置 添加新列 region 标记UTR和CDS
def utr_mark(nc_no: str) -> None:
    exon_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/EXON.tsv"
    cds_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/CDS.tsv"
    ori_dict = gene_ori_dict(nc_no)
    utr_pos_dict = utr_region_obtain(exon_file, cds_file, ori_dict)
    ...
    

def main() -> None:
    nc_no = "NC_000024.10"
    
    utr_mark(nc_no)


if __name__ == "__main__":
    main()
