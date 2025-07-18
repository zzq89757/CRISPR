

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
# left_gene_li['filter'] = [x for x in gene_li]
counter_raw = Counter()
# counter_filter = Counter()
# each nc_no
nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
nc_li = nc_df[0].tolist()
for nc_no in nc_li:
    print(f"{nc_no} start !!!")
    # nc_no = "NC_000024.10"
    type_li =  []
    gdb_df = tsv2df(f"/mnt/ntc_data/wayne/Repositories/CRISPR/dual_cd/{nc_no}.tsv",type_li)

    for gene, sub_df in gdb_df.groupby(4):
        # raw sub df
        left_gene_li["raw"].remove(gene)
        # counter raw
        raw_count = len(sub_df)
        raw_rank_bottom = (raw_count // 10) * 10 + 1
        raw_rank_top = ((raw_count // 10) + 1) * 10
        rank_str = f"{raw_rank_bottom}-{raw_rank_top}" if raw_rank_bottom < 40 else "40+"
        counter_raw[rank_str] += 1
        # filter sub df
        # filtered_sub_df = sub_df[(sub_df[18]>0.1)&(sub_df[22]>0.3)]
        # counter filter
        # filter_count = len(filtered_sub_df)
        # if filter_count == 0:
            # counter_filter['0'] += 1
            # continue
        
        # left_gene_li["filter"].remove(gene)
        
        # raw_rank_bottom = (filter_count // 10) * 10 + 1
        # raw_rank_top = ((filter_count // 10) + 1) * 10
        rank_str = f"{raw_rank_bottom}-{raw_rank_top}" if raw_rank_bottom < 40 else "40+"
        # counter_filter[rank_str] += 1


# output 0 gene li
raw_0_li = open("raw_0.sli",'w')
# filter_0_li = open("filter_0.sli",'w')

raw_0_li.write("\n".join(left_gene_li['raw']))
# filter_0_li.write("\n".join(left_gene_li['filter']))

# add 0 gene count 
counter_raw['0'] += len(left_gene_li['raw'])
# counter_filter['0'] += len(left_gene_li['filter'])
print(counter_raw)
# print(counter_filter)