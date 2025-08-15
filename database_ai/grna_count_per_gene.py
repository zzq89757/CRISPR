

# 读取Gene list ，将gdb group by gene 若去除一柱和蓝色后仍有结果 则将其从gene list中去除 同时统计grna数目

# 原始结果统计

from collections import defaultdict,Counter
import pandas as pd
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df


gene_li = pd.read_csv("/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/Gene_list.tsv",sep="\t",header=None)[0].to_list()

left_gene_li = defaultdict(list)

left_gene_li['raw'] = [x for x in gene_li]
left_gene_li['low_filter'] = [x for x in gene_li]
left_gene_li['high_filter'] = [x for x in gene_li]
counter_raw = Counter()
counter_low_filter = Counter()
counter_high_filter = Counter()
# each nc_no
nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
nc_li = nc_df[0].tolist()

max_num = 0 
for nc_no in nc_li:
    print(f"{nc_no} start !!!")
    # nc_no = "NC_000024.10"
    type_li =  ["string", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category", "float", "int32"]
    gdb_df = tsv2df(f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/rs2_score/{nc_no}.tsv",[])

    for gene, sub_df in gdb_df.groupby(7):
        # print(gene)
        # raw sub df
        left_gene_li["raw"].remove(gene)
        # counter raw
        raw_count = len(sub_df)
        
        if raw_count > max_num:
            max_num = raw_count
        
        raw_rank_bottom = (raw_count // 10) * 10 + 1
        raw_rank_top = ((raw_count // 10) + 1) * 10
        rank_str = f"{raw_rank_bottom}-{raw_rank_top}" if raw_rank_bottom < 50 else "50+"
        counter_raw[rank_str] += 1
        # low filter sub df
        filtered_sub_df = sub_df[(sub_df[13]>0.3)&(sub_df[17]>0.1)]
        
        # counter low filter
        filter_count = len(filtered_sub_df)
        if filter_count != 0:
            # counter_filter['0'] += 1
            # continue
        
            left_gene_li["low_filter"].remove(gene)
            
            raw_rank_bottom = (filter_count // 10) * 10 + 1
            raw_rank_top = ((filter_count // 10) + 1) * 10
            rank_str = f"{raw_rank_bottom}-{raw_rank_top}" if raw_rank_bottom < 50 else "50+"
            counter_low_filter[rank_str] += 1
        # filter sub df
        filtered_sub_df = sub_df[(sub_df[13]>=0.5)&(sub_df[17]>=0.5)]
        
        
        # counter filter
        filter_count = len(filtered_sub_df)
        if filter_count != 0:
            # counter_filter['0'] += 1
            # continue
        
            left_gene_li["high_filter"].remove(gene)
            
            raw_rank_bottom = (filter_count // 10) * 10 + 1
            raw_rank_top = ((filter_count // 10) + 1) * 10
            rank_str = f"{raw_rank_bottom}-{raw_rank_top}" if raw_rank_bottom < 50 else "50+"
            counter_high_filter[rank_str] += 1


# output 0 gene li
# raw_0_li = open("raw_5.sli",'w')
# filter_0_li = open("filter_5.sli",'w')

# raw_0_li.write("\n".join(left_gene_li['raw']))
# filter_0_li.write("\n".join(left_gene_li['filter']))

# add 0 gene count 
counter_raw['0'] += len(left_gene_li['raw'])
counter_low_filter['0'] += len(left_gene_li['low_filter'])
counter_high_filter['0'] += len(left_gene_li['high_filter'])
print(counter_raw)
print(counter_low_filter)
print(counter_high_filter)
# print(f"Max num is {max_num}")