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
            
            # 开始分情况计算UTR区域
            cds_region_count = len(cds_start_li)
            exon_region_count = len(exon_start_li)
            
            # 记录cds 5' 和 3' 位置
            cds_terminal5 = cds_start_li[0] if gene_ori == "+" else cds_end_li[0]
            cds_terminal3 = cds_end_li[-1] if gene_ori == "+" else cds_start_li[-1]
            
            # 方向为正的情况(默认)
            utr5_start = int(exon_start_li[0])
            utr5_end = int(cds_start_li[0]) - 1
            utr3_start = int(cds_end_li[cds_region_count - 1]) + 1 if cds_region_count <= exon_region_count else 0
            print(exon_tran)
            # print(exon_end_li)
            utr3_end = int(exon_end_li[cds_region_count - 1]) if cds_region_count <= exon_region_count else 0
            
            # 用 5' 和 3' 位置分别从Exon起始终止开始扫描
            for i in range(0,exon_region_count):
                exon_start = exon_start_li[i]
                exon_end = exon_end_li[i]
                # 若无UTR5（cds_terminal5 == exon_start OR cds_terminal5 == exon_end） 则不存入区间
                if cds_terminal5 == exon_start or cds_terminal5 == exon_end:
                    # print(exon_region_dict[gene]["UTR5"])
                    # print(exon_tran)
                    break
                # 直至寻找到cds 5' 所在的Exon区间 
                if exon_start < cds_terminal5 < exon_end:
                    # 计算UTR5起始终止
                    utr5_start = exon_start if gene_ori == "+" else cds_terminal5 + 1
                    utr5_end = cds_terminal5 - 1 if gene_ori == "+" else exon_end
                    exon_region_dict[gene]["UTR5"].append((utr5_start, utr5_end))
                    break
                # 将已经遍历过的区间加入extra utr5 region
                exon_region_dict[gene]["UTR5"].append((exon_start, exon_end))
            
            stop_codon_index = 0
            stop_codon_pos = cds_terminal3 + 3 if gene_ori == "+" else cds_terminal3 - 3
            # 3'需要考虑终止密码子是否与cds3’落在同一外显子区间
            for j in range(exon_region_count - 1, -1, -1):
                exon_start = exon_start_li[j]
                exon_end = exon_end_li[j]
                
                # 记录cds_terminal3 所在区间 
                if exon_start <= cds_terminal3 <= exon_end:
                    # 判断终止密码子落点是否在当前Exon区间
                    if cds_region_count == exon_region_count:
                        utr3_start = stop_codon_pos + 1 if gene_ori == "+" else exon_start
                        utr3_end = exon_end if gene_ori == "+" else stop_codon_pos - 1
                        # 计算UTR3起始终止
                        exon_region_dict[gene]["UTR3"].append((utr3_start, utr3_end))
                        break
                    else:
                        if exon_start < stop_codon_pos < exon_end:
                            utr3_start = stop_codon_pos + 1 if gene_ori == "+" else exon_start
                            utr3_end = exon_end if gene_ori == "+" else stop_codon_pos - 1
                            # 计算UTR3起始终止
                            exon_region_dict[gene]["UTR3"].append((utr3_start, utr3_end))
                            break
                        else:
                            
                            # print(exon_tran)
                            # 终止密码子落点不在当前Exon区间
                            # 计算区间剩余
                            over_len = stop_codon_pos - exon_end if gene_ori == "+" else exon_start - stop_codon_pos
                            
                            # 若无UTR3 （over_len==0）
                            if over_len == 0:
                                break
                            
                            # 添加跨区碱基
                            
                            stop_codon_index = j + 1
                            # 定位终止密码子区间位置
                            utr3_start = cds_terminal3 + 1 if gene_ori == "+" else exon_start
                            if gene_ori == "+":
                                utr3_start = exon_start_li[stop_codon_index] + over_len - 1
                                utr3_end = exon_end_li[stop_codon_index]
                            else:
                                utr3_start = exon_start_li[stop_codon_index]
                                utr3_end = exon_end_li[stop_codon_index] - over_len # +1 -1 抵消
                            exon_region_dict[gene]["UTR3"].remove((exon_start_li[stop_codon_index], exon_end_li[stop_codon_index]))                           
                            exon_region_dict[gene]["UTR3"].append((utr3_start, utr3_end))
                            # print(stop_codon_pos)
                            # print(f"current:({utr3_start},{utr3_end})")
                            # print(f"all:({exon_region_dict[gene]["UTR3"]}")
                            break
                # 将已经遍历过的区间加入extra utr3 region
                exon_region_dict[gene]["UTR3"].append((exon_start, exon_end))  
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
    # 若为0 直接返回空字符串
    if not region_code:return ""
    # 将代表切点落在的区域之和的十进制数字(1->utr5,2->cds,4->utr3)转为二进制
    bin_region_code = '{:03b}'.format(region_code)[::-1]
    # -1->utr5,0->cds,1->utr3
    mark_li = ["-1", "0", "1"]
    res_li = [mark_li[i] for i in range(3) if int(bin_region_code[i])]
    return ",".join(res_li)


def region_classify(gene_name: str, cut_pos: int, grna_ori: str, utr_pos_dict: defaultdict) -> str:
    '''
    1->utr5,2->cds,4->utr3
    -1->utr5,0->cds,1->utr3
    '''
    region_code = 0
    # 根据方向添加切点偏移量
    cut_pos_offset = cut_pos + float(grna_ori + "0.5")
    cds_region_li = utr_pos_dict[gene_name]["CDS"]
    if not cds_region_li:return ""
    utr5_region_li = utr_pos_dict[gene_name]["UTR5"]
    utr3_region_li = utr_pos_dict[gene_name]["UTR3"]
    # 遍历每个区域 进行十进制数字累加
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
    gdb_df["Target Region"] = gdb_df.apply(lambda x:region_classify(x[8],x[7],x[4],region_dict),axis=1)
    return gdb_df   


# 根据UTR位置 添加新列 region 标记UTR和CDS
def utr_mark(nc_no: str) -> None:
    t1 = time.time()
    print(f"{nc_no} start !!!")
    exon_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/EXON.tsv"
    cds_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/CDS.tsv"
    gdb_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/snp_mark/{nc_no}.tsv"
    res_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/utr_mark/{nc_no}.tsv"
    # if Path(res_file).exists():
    #     print(f"{nc_no} is exist,exit !!!")
    #     return
    ori_dict = gene_ori_dict(nc_no)
    utr_pos_dict = utr_region_obtain(exon_file, cds_file, ori_dict) 
    res_df = search_regions(gdb_file, utr_pos_dict)
    res_df.to_csv(res_file,sep="\t",header=False,index=False)
    print(f"{nc_no} time cost:{time.time() - t1}")


def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(utr_mark,nc_li,24)
    # utr_mark(nc_li[0])


if __name__ == "__main__":
    main()
