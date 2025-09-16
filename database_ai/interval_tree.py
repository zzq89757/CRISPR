import pandas as pd
from intervaltree import Interval, IntervalTree
from time import time
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from process_border_withid import async_in_iterable_structure, relative_pos_calc

def run(nc_no: str) -> None:
    t1 = time()
    # === 1. 读取基因区间文件 ===
    gene_df = pd.read_csv(f"tss_tran/{nc_no}.tsv", sep="\t", header=None, names=["trans", "tss", "gene", "gene_ori", "gene_start", "gene_end", "gene_id", "gene_type"])

    # 构建区间树
    tree = IntervalTree()
    for idx, row in gene_df.iterrows():
        tran_str = row["trans"]
        gene = row["gene"]
        gene_id = row["gene_id"]
        gene_type = row["gene_type"]
        strand = row["gene_ori"]
        gene_start = row["gene_start"]
        gene_end = row["gene_end"]
        left = row["tss"] - 1000 - 11 if strand == "+" else row["tss"] + 1 - 11
        right = row["tss"] - 1 + 11 if strand == "+" else row["tss"] + 1000 + 11
        hit_str = f"{tran_str}\t{gene_start}\t{gene_end}\t{gene}\t{gene_id}\t{strand}\t{gene_type}"
        tree.add(Interval(left, right, hit_str))  # 注意end+1


    # === 2. 读取gRNA文件 ===
    grna_df = pd.read_csv(f"nag_remove/{nc_no}.tsv", sep="\t", header=None)



    # === 3. 批量判断命中的基因 ===
    def check_gene(pos):
        hits = tree[pos]
        if hits:
            return ";".join(sorted(set(iv.data for iv in hits)))
        else:
            return None
    # 找到grna中点
    grna_offset = (grna_df[4].astype(str) + "3").astype(int)
    grna_df["grna_middle"] = ((grna_df[6] + grna_offset + grna_df[5]) / 2).astype(int)

    grna_df["Hit_Gene"] = grna_df["grna_middle"].apply(check_gene)
    # 过滤未hit的
    grna_df = grna_df[grna_df["Hit_Gene"].notna()]
    # 如果 Hit_Gene 是以逗号分隔的字符串，先拆分成列表
    grna_df["Hit_Gene"] = grna_df["Hit_Gene"].str.split(";")

    # 拆成多行
    grna_df = grna_df.explode("Hit_Gene").reset_index(drop=True)

    # 拆分信息列
    grna_df[["trans", "gene_start", "gene_end", "gene", "gene_id", "gene_strand", "gene_type"]] = grna_df["Hit_Gene"].str.split("\t", expand=True)
    grna_df["gene_start"] = grna_df["gene_start"].astype(int)
    grna_df["gene_end"] = grna_df["gene_end"].astype(int)

    # 计算相对基因的位置、方向、切点相对位置
    grna_df["re_dir"] = grna_df.apply(
        lambda row: "fwd" if row[4] == row["gene_strand"] else "rev", axis=1
    )


    grna_df["grna_left"] = grna_df["grna_middle"] - 11
    grna_df["grna_right"] = grna_df["grna_middle"] + 11

    grna_df["re_loc"] = pd.DataFrame(
        grna_df.apply(
            lambda row: relative_pos_calc(row[5], row[6], row[7], row[4], row["gene_start"], row["gene_end"], row["gene_strand"])[0],
            axis=1
        ).tolist(),  # 重要：先转成列表，再传给DataFrame
        index=grna_df.index  # 保证索引一致
    )
    grna_df["re_loc"] = grna_df["re_loc"].apply(
        lambda x: x.replace("--", "-(-") + ")" if "--" in x else x
    )
    # 重设表头
    new_header = list(range(7)) + ["gene", "gene_id", "gene_type", "re_dir", "re_loc", "trans"]
    grna_df = grna_df[new_header]
    # 合并相同id 相同基因的行
    grna_df = grna_df.groupby(list(range(7)) + ["gene", "gene_id", "gene_type", "re_dir", "re_loc"], as_index=False).agg({
    "trans": lambda x: ",".join(sorted(set(x))),
    # 其他列也按需添加，假设这些列在每组中是一致的
})
    # === 4. 保存或显示结果 ===
    # print(grna_df)

    # 可选：保存结果
    grna_df.to_csv(f"gene_annotated/{nc_no}.tsv", sep="\t", index=False, header=None)
    print(f"<{nc_no}> annotated,time cost: {time() - t1}")



def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run,nc_li,24)
    # run(nc_li[-1])

if __name__ == "__main__":
    main()