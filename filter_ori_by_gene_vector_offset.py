import time
import pandas as pd


def annotation_gdb(gdb: pd.DataFrame, gene_pos_df: pd.DataFrame) -> pd.DataFrame:
    # 需要分别考虑切点和起止 三个位置 超出基因末端时 loc 采用(+/- len(over))
    # 提取必要列为 NumPy 数组
    gene_start_array = gene_pos_df[1].to_numpy()
    gene_end_array = gene_pos_df[2].to_numpy()
    gRNA_ori_array = gdb[3].to_numpy()
    gRNA_pos_array = gdb[6].to_numpy()
    gdb_array = gdb.to_numpy()
    gene_pos_array = gene_pos_df.to_numpy()
    gene_pos_len = len(gene_start_array)
    current_gene_pos_idx = 0
    max_end = max(gene_end_array)
    min_start = min(gene_end_array)
    ## 遍历 gRNA <考虑gRNA不同方向时切点位置不同>
    for g_idx, g_raw_pos in enumerate(gRNA_pos_array):
        g_ori = gRNA_ori_array[g_idx]
        # 分类考虑正负链的情况
        g_pos_offset = float(g_ori + "0.5")
        g_pos = g_raw_pos + g_pos_offset
        # 跳过位于当前基因起始位置之前以及基因间的的 gRNA
        if g_pos < min_start or g_pos > max_end:
            continue

        # 跨多个基因时的情况
        # while current_gene_pos_idx < gene_pos_len - 1 and g_pos > gene_end:
        #     current_gene_pos_idx += 1
        #     gene_start = gene_start_array[current_gene_pos_idx]
        #     gene_end = gene_end_array[current_gene_pos_idx]

        # 大于max end时 结束
        # if current_gene_pos_idx == gene_pos_len - 1 and g_pos > max_end:
        #     break

        # 如果 gRNA 位于当前基因范围内，记录信息 chr1:1335323 跨两个或以上的位点会被跳过
        for current_gene_pos_idx in range(gene_pos_len):
            gene_start = gene_start_array[current_gene_pos_idx]
            gene_end = gene_end_array[current_gene_pos_idx]
            if gene_start < g_pos < gene_end:
                # print(f"{g_pos} in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}", flush=True, end="\r")
                # print(f"{g_pos} in {gene_start}-{gene_end}")
                # g_item = gdb.iloc[g_idx]
                g_item = gdb_array[g_idx]
                # 浮动的切点可能错过子区间间区而不是基因间区
                if gene_pos_array[current_gene_pos_idx][0].find(",") != -1:
                    print(g_item)
                ...
            else:
                # print(f"{g_pos} {gene_start_array[current_gene_pos_idx - 1]}-{gene_end_array[current_gene_pos_idx - 1]} and {gene_start}-{gene_end} and {gene_start_array[current_gene_pos_idx + 1]}-{gene_end_array[current_gene_pos_idx + 1]}")
                ...

    return pd.DataFrame([])


def main() -> None:
    t = time.time()
    nc_no = "NC_000024.10"
    chr = "chrY"

    nc_table = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    df = pd.read_csv(nc_table, sep="\t", header=None)
    chr2nc_dict = df.set_index(1)[0].to_dict()
    type_li = ["string", "int32", "int32", "category", "string", "category"]
    type_dict = dict(enumerate(type_li))
    # 读取基因位置信息文件
    gene_pos_df = pd.read_csv(
        # f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/{chr2nc_dict[chr]}/Gene_list_cut_insertion.tsv",
        f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/{chr2nc_dict[chr]}/Gene_list.tsv",
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
    # gdb_df = pd.read_csv(
    #     f"/mnt/ntc_data/wayne/Repositories/CRISPR/y_10w.tsv",
    #     header=None,
    #     sep="\t",
    #     dtype=type_dict,
    # )
    # gdb 注释
    annotated_df = annotation_gdb(gdb_df, gene_pos_df)
    # 保存为tsv文件
    # annotated_df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/{chr}_gene.tsv", index=False, header=None, sep="\t")
    print(f"time cost:{time.time() - t}")


if __name__ == "__main__":
    main()
