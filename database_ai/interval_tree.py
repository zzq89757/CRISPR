import pandas as pd
from intervaltree import Interval, IntervalTree
from time import time

t1 = time()
# === 1. 读取基因区间文件 ===
gene_df = pd.read_csv("tss_regions/NC_000024.10.tsv", sep="\t", header=None, names=["gene", "strand", "regions"])

# 构建区间树
tree = IntervalTree()
for idx, row in gene_df.iterrows():
    gene_name = row["gene"]
    strand = row["strand"]
    region_list = row["regions"].split(",")
    for region in region_list:
        start, end = map(int, region.split("-"))
        tree.add(Interval(start, end + 1, gene_name))  # 注意end+1


# === 2. 读取gRNA文件 ===
grna_df = pd.read_csv("../split_out/sorted/no_head/spCas9_Homo_chrY.tsv", sep="\t", header=None)
# 假设位点在最后一列
pos_col_idx = grna_df.shape[1] - 1


# === 3. 批量判断命中的基因 ===
def check_gene(pos):
    hits = tree[pos]
    if hits:
        return ",".join(sorted(set(iv.data for iv in hits)))
    else:
        return "None"


grna_df["Hit_Gene"] = grna_df[pos_col_idx].apply(check_gene)


# === 4. 保存或显示结果 ===
print(grna_df)

# 可选：保存结果
grna_df.to_csv("grna_annotated.tsv", sep="\t", index=False)
print(f"time cost: {time() - t1}")