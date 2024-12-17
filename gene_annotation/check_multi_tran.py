import pandas as pd
from collections import defaultdict




def gtf2df(gtf_file: str, nc_li: list) -> pd.DataFrame:
    """提前预设每列的数据类型并将gtf文件存入DataFrame"""
    use_col_li = [0, 2, 8]
    type_li = ["string", "category", "string"]

    type_dict = dict(zip(use_col_li, type_li))

    gtf_df = pd.read_csv(
        gtf_file,
        sep="\t",
        header=None,
        usecols=use_col_li,
        dtype=type_dict,
        # low_memory=False,
        comment="#",
    )
    gtf_df[9] = gtf_df[8].str.extract(r'gene "([^"]+)"')
    # print(gtf_df[9])
    gtf_df[10] = gtf_df[8].str.split('"').str[3]
    gtf_df = gtf_df[(gtf_df[0].isin(nc_li) & (gtf_df[2] == "transcript") & gtf_df[10].str.contains("N"))][[0, 9, 10]]
    # print(gtf_df)
    return gtf_df


def check_multi_transcript(gtf_df: pd.DataFrame) -> None:
    gene_tran_dict = defaultdict(list)
    for gene, tran in zip(gtf_df[9],gtf_df[10]):
        gene_tran_dict[gene].append(tran)
    two_more_tran_dict = {i:v for i,v in gene_tran_dict.items() if len(v) > 1}
    print(two_more_tran_dict)
    
    


def main() -> None:
    gtf_file = "GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.gtf"
    nc2chr_file = "nc2chr.tsv"
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].to_list()
    gtf_df = gtf2df(gtf_file, nc_li)
    check_multi_transcript(gtf_df)

if __name__ == "__main__":
    main()