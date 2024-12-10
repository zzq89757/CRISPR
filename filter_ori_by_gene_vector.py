import time
import pandas as pd


def annotation_gdb(gdb: pd.DataFrame, gene_pos_df: pd.DataFrame) -> pd.DataFrame:
    # 需要分别考虑切点和起止 三个位置 超出基因末端时 loc 采用(+/- len(over))
    # 提取必要列为 NumPy 数组
    gene_start = gene_pos_df[1].to_numpy()
    gene_end = gene_pos_df[2].to_numpy()
    gRNA_pos = gdb[6].to_numpy()
    gene_pos_len = len(gene_start)
    current_gene_pos_idx = 0

    ## 遍历 gRNA
    for g_idx, g_pos in enumerate(gRNA_pos):
        # 跳过位于当前基因起始位置之前以及基因间的的 gRNA
        if g_pos < gene_start[current_gene_pos_idx]:
            continue
      
        # 跨多个基因时的情况
        while (current_gene_pos_idx < gene_pos_len -1  and g_pos > gene_end[current_gene_pos_idx]):
            current_gene_pos_idx += 1
        
        # 大于max end时 结束
        if current_gene_pos_idx == gene_pos_len -1 and g_pos > gene_end[-1]:
            break
        
        # 如果 gRNA 位于当前基因范围内，记录信息 chr1:1335323 跨两个或以上的位点会被跳过
        if g_pos <= gene_end[current_gene_pos_idx]:
            # print(f"{g_pos} in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}", flush=True, end="\r")
            print(f"{g_pos} in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}")
            ...
        else:
            print(f"{g_pos} not in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}")

    return pd.DataFrame([])


def main() -> None:
    t = time.time()
    nc_no = "NC_000024.10"
    chr = "chr1"

    nc_table = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    df = pd.read_csv(nc_table, sep="\t", header=None)
    chr2nc_dict = df.set_index(1)[0].to_dict()
    type_li = ["string", "int32", "int32", "category", "string", "category"]
    type_dict = dict(enumerate(type_li))
    # 读取基因位置信息文件
    gene_pos_df = pd.read_csv(
        f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/{chr2nc_dict[chr]}/Gene_list_cut_insertion.tsv",
        sep="\t",
        header=None,
        low_memory=False,
        dtype=type_dict,
    )
    # 读取gRNA db 文件
    type_li = ["string", "string", "category", "category", "int32", "int32", "int32"]
    type_dict = dict(enumerate(type_li))
    gdb_df = pd.read_csv(
        f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/spCas9_Homo_{chr}.tsv",
        header=None,
        sep="\t",
        dtype=type_dict
    )
    # gdb_df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/y_10w.tsv", header=None, sep="\t")
    # gdb 注释
    annotated_df = annotation_gdb(gdb_df, gene_pos_df)
    # 保存为tsv文件
    # annotated_df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/{chr}_gene.tsv", index=False, header=None, sep="\t")
    print(f"time cost:{time.time() - t}")


if __name__ == "__main__":
    main()
