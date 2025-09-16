import time
import pandas as pd


def gdb2df(gdb_path: str) -> pd.DataFrame:
    type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32"]

    type_dict = dict(enumerate(type_li))

    gdb_df = pd.read_csv(
        gdb_path,
        sep="\t",
        header=None,
        dtype=type_dict,
        # low_memory=False,
    )
    # print(gdb_df)
    return gdb_df


def region2df(region_file: str) -> pd.DataFrame:
    type_li = ["string", "string", "int32", "int32"]

    type_dict = dict(enumerate(type_li))

    region_df = pd.read_csv(
        region_file,
        sep="\t",
        header=None,
        dtype=type_dict,
        # low_memory=False,
    )
    # print(region_df)
    return region_df


def annotation_gdb(gdb_df: pd.DataFrame, cds_df: pd.DataFrame, exon_df: pd.DataFrame) -> None:
    # 根据gdb 的基因和切点 找到cut的转录本、外显子以及cds并集前2/3位置
    # 提取所需列并矢量化计算
    gdb_ori_array = gdb_df[4].astype(str) + "0.5"
    gdb_cut_pos_array = gdb_df[7].to_numpy()
    gdb_gene_name_array = gdb_df[8].to_numpy()

    # 计算 g_pos 的矢量化结果
    g_pos_array = gdb_cut_pos_array + gdb_ori_array.astype(float)

    # 筛选 exon 数据
    result = []
    res_count = 0
    idx = 0
    for gene_name, g_pos in zip(gdb_gene_name_array, g_pos_array):
        idx += 1
        # 一次性筛选出满足条件的 exon 数据
        sub_exon_df = exon_df[
            (exon_df[0] == gene_name) & (g_pos > exon_df[2]) & (g_pos < exon_df[3])
        ]
        print(f"{idx}/{len(gdb_gene_name_array)} finished !!!",end="\r",flush=True)
        if sub_exon_df.empty:
            res_count += 1
    print(f"not fall in exon num:{res_count}/{len(gdb_gene_name_array)} !!!")


def main() -> None:
    t1 = time.time()
    nc2chr_file = "nc2chr.tsv"
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc2chr_dict = dict(zip(nc_df[0],nc_df[1]))
    # 读取gdb
    gdb_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/gene_annotation/spCas9_Homo_chrY.tsv"
    gdb_df = gdb2df(gdb_file)
    # 读取对应的CDS exon 表
    cds_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/NC_000024.10/CDS.tsv"
    cds_df = region2df(cds_file)
    exon_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/NC_000024.10/EXON.tsv"
    exon_df = region2df(exon_file)
    print(f"load gdb time cost:{time.time() - t1}")
    # 注释gdb
    t1 = time.time()
    annotation_gdb(gdb_df, cds_df, exon_df)
    print(f"annotation time cost:{time.time() - t1}")
    
    


if __name__ == "__main__":
    main()