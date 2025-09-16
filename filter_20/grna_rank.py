

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
counter_raw = Counter()
# each nc_no
nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
nc_li = nc_df[0].tolist()

for nc_no in nc_li:
    print(f"{nc_no} start !!!")
    # nc_no = "NC_000024.10"
    type_li =  ["string", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category", "float", "int32"]
    gdb_df = tsv2df(f"/mnt/ntc_data/wayne/Repositories/CRISPR/filter_20/{nc_no}.tsv",[])

    for gene, sub_df in gdb_df.groupby(9):
        # raw sub df
        # print(gene)
        left_gene_li["raw"].remove(gene)
        # counter raw
        raw_count = len(sub_df)        
        counter_raw[raw_count] += 1
        


# output 0 gene li
raw_0_li = open("raw_5.sli",'w')

raw_0_li.write("\n".join(left_gene_li['raw']))

# add 0 gene count 
counter_raw[0] = len(left_gene_li['raw'])

for i in range(21):
    print(f"{i}:{counter_raw[i]}",end=", ")