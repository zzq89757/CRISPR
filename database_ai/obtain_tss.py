import pandas as pd
import numpy as np
from collections import defaultdict

from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from generate_split_ori import async_in_iterable_structure



def obtain_gene_info_dict(nc_no: str) -> dict:
    gene_df: pd.DataFrame = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/Gene_list.tsv",sep="\t",header=None)
    # gene_df["5_end"] = np.where(gene_df.iloc[:, 3] == '+', gene_df.iloc[:, 1], gene_df.iloc[:, 2])
    # gene_df["3_end"] = np.where(gene_df.iloc[:, 3] == '+', gene_df.iloc[:, 2], gene_df.iloc[:, 1])
    gene_dict = defaultdict(list)

    for _, row in gene_df.iterrows():
        symbol = row[0]
        gene_dict[symbol] = [row[3], row[1], row[2], row[4], row[5]]
    return gene_dict


# 根据split gtf下对应染色体的Gene_list.tsv(获取基因方向)和TRAN.tsv(获取基因对应转录本及其起止) 提取基因对应的tss
def get_gene_transcript_starts_from_df(nc_no: str, gene_ori_dict: dict) -> defaultdict:
    """
    根据基因方向，从DataFrame中返回该基因所有转录本及其起始位点，合并重复并排序。
    """
    gene_tss_dict = defaultdict(dict)
    df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/TRAN.tsv",sep="\t",header=None)
    for gene, sub_df in df.groupby(0,sort=False):
        ori = gene_ori_dict[gene][0]
        sub_df["5_end"] = sub_df[2] if ori == "+" else sub_df[3]
        gene_tss_dict[gene] = list(
                sub_df.groupby("5_end")[1]
                .apply(list)
                .items()
            )
    return gene_tss_dict


def export_tss_dict(nc_no: str, tss_dict: defaultdict, gene_info_dict: defaultdict) -> None:
    # 构建DataFrame用的列表
    records = []

    for gene, tuples in tss_dict.items():
        for position, transcript_list in tuples:
            transcript_str = ",".join(transcript_list)
            gene_ori, gene_5end_pos, gene_3end_pos, gene_id, gene_type = gene_info_dict[gene]
            records.append([transcript_str, position, gene, gene_ori, gene_5end_pos, gene_3end_pos, gene_id, gene_type])

    # 转为DataFrame
    df = pd.DataFrame(records)
    # 按照tss排序
    df.sort_values(1,inplace=True)
    # 导出TSV
    df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/tss_tran/{nc_no}.tsv", index=False, sep="\t", header=None)


def extract_tss(nc_no: str) -> defaultdict:
    # 获取基因方向和起点的字典
    gene_info_dict = obtain_gene_info_dict(nc_no)
    # print(gene_info_dict)
    # 获取基因TSS位点
    tss_dict = get_gene_transcript_starts_from_df(nc_no, gene_info_dict)
    # 构造为df并导出
    export_tss_dict(nc_no, tss_dict, gene_info_dict)
    

def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(extract_tss,nc_li,24)
    # extract_tss(nc_li[-1])
    

if __name__ == "__main__":
    main()
    