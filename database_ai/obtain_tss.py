import pandas as pd
from collections import defaultdict

from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from generate_split_ori import async_in_iterable_structure
from cds_mark import gene_ori_dict


# 根据split gtf下对应染色体的Gene_list.tsv(获取基因方向)和TRAN.tsv(获取基因对应转录本及其起止) 提取基因对应的tss
def get_gene_transcript_starts_from_df(nc_no: str, gene_ori_dict: dict) -> defaultdict:
    """
    根据基因方向，从DataFrame中返回该基因所有转录本的起始位点，去重并排序。
    """
    gene_tss_dict = defaultdict(dict)
    df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/TRAN.tsv",sep="\t",header=None)
    for gene, sub_df in df.groupby(0,sort=False):
        ori = gene_ori_dict[gene]
        if ori == '+':
            # 构造字典
            gene_tss_dict[gene] = list(
                sub_df.groupby(2)[1]
                .apply(list)
                .items()
            )
            # positions = sub_df[2].unique()
        else:
            gene_tss_dict[gene] = list(
                sub_df.groupby(3)[1]
                .apply(list)
                .items()
            )
            # positions = sub_df[3].unique()
        # gene_tss_dict[gene] = sorted(positions.tolist())
    return gene_tss_dict


def export_tss_dict(nc_no: str, tss_dict: defaultdict) -> None:
    # 构建DataFrame用的列表
    records = []

    for gene, tuples in tss_dict.items():
        for position, transcript_list in tuples:
            transcript_str = ",".join(transcript_list)
            records.append([gene, position, transcript_str])

    # 转为DataFrame
    df = pd.DataFrame(records)

    # 查看
    # print(df)

    # 导出CSV（可选）
    df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/tss_tran/{nc_no}.tsv", index=False, sep="\t", header=None)


def extract_tss(nc_no: str) -> defaultdict:
    # 获取基因方向字典
    ori_dict = gene_ori_dict(nc_no)
    # 获取基因TSS位点
    tss_dict = get_gene_transcript_starts_from_df(nc_no, ori_dict)
    # 构造为df并导出
    export_tss_dict(nc_no, tss_dict)
    # 取每个TSS的上游1kb并取并集
    # print(tss_dict)
    

def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(extract_tss,nc_li,24)
    # nc_no = "NC_000024.10"
    # extract_tss(nc_no)
    

if __name__ == "__main__":
    main()
    