import time
import pandas as pd
from collections import defaultdict, Counter
from pathlib import Path
from multiprocessing import Pool, RLock
from sys import argv


def async_in_iterable_structure(fun, iterable_structure, cpus):
    def init(l):
        global lock
        lock = l

    lock = RLock()
    p = Pool(int(cpus), initializer=init, initargs=(lock,))
    # apply async in iterable structure
    # for i in iterable_structure:
    #     p.apply_async(func=fun, args=(i,))
    results = [p.apply_async(func=fun, args=(i,)) for i in iterable_structure]
    p.close()
    p.join()
    # return [result.get() for result in results]


def gtf2df(gtf_file: str) -> None:
    """提前预设每列的数据类型并将gtf文件存入DataFrame"""
    use_col_li = [0, 2, 3, 4, 6, 8]
    type_li = ["string", "category", "int32", "int32", "category", "string"]

    type_dict = dict(zip(use_col_li, type_li))

    gtf_df = pd.read_csv(
        gtf_file,
        sep="\t",
        header=None,
        usecols=use_col_li,
        dtype=type_dict,
        low_memory=False,
        comment="#",
    )

    return gtf_df


def filter_notchr(nc2chr_dict: dict, gtf_df: pd.DataFrame) -> None:
    """过滤掉不在24个染色体上的gene item"""
    gtf_df = gtf_df[gtf_df[0].isin(nc2chr_dict.keys())]


def append_gene_id_col(gtf_df: pd.DataFrame) -> None:
    """为gtf添加gene id 和 tran id 列"""
    gtf_df[9] = gtf_df[8].str.split('"').str[1]
    gtf_df[10] = gtf_df[8].str.split('"').str[3]


def obtain_gene_id_li(gtf_df: pd.DataFrame) -> defaultdict:
    """挑出gene行,只保留gene biotype 为protein_coding和ncRNA的gene id 列表并统计数目"""
    gene_df = gtf_df[gtf_df[2] == "gene"]
    contains_bool_li1 = gene_df[8].str.contains(r"gene_biotype \"protein_coding\";")
    contains_bool_li2 = gene_df[8].str.contains(r"gene_biotype \"ncRNA\";")
    contains_bool_li3 = gene_df[8].str.contains(r"gene_biotype \"lncRNA\";")
    contains_bool_li = [x or y or z for x, y, z in zip(contains_bool_li1, contains_bool_li2, contains_bool_li3)]
    gene_df = gene_df[contains_bool_li]

    gene_id_li = gene_df[9].to_list()
    # 创建拆分的gtf的路径
    # for chr_name, sub_df in gene_df.groupby(0):
        # 先不创建 因为有的基因无NM和NR
        # for gene in sub_df[9]:
        #     Path(f"split_gtf/{chr_name}/{gene}").mkdir(exist_ok=1,parents=1)
        # gene_df[[9,3,4,6]].to_csv(f"split_gtf/{chr_name}/{chr_name}.tsv",header=None,index=False,sep="\t")
    return [gene_id_li, sum(contains_bool_li1), sum(contains_bool_li2)]


def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    gtf_file = "GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.gtf"
    # gtf_file = "1w.gtf"
    num_count = Counter()
    info_dict = defaultdict()
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc2chr = dict(zip(nc_df[0], nc_df[1]))
    # 记录程序开始时间
    t = time.time()
    # 读取gtf文件 生成df
    gtf_df = gtf2df(gtf_file)
    # 过滤非染色体的条目
    filter_notchr(nc2chr, gtf_df)
    # 添加gene id 和tran id列
    append_gene_id_col(gtf_df)
    # 找到coding和uncoding的gene id 列表
    gene_id_li, num_count["coding_num"], num_count["noncoding_num"] = (
        obtain_gene_id_li(gtf_df)
    )
    # 按照gene id group by 只保留gene id 列表中且包含NM和NR的转录本（统计数目）
    # gb = gtf_df.groupby(0,sort=False)
    # print(sub_df)
    # return 
    def cut_gtf(i):
        sub_df = gtf_df[gtf_df[0]==i]
        chr_name = i
        # sub_df = gb.get_group(chr_name)
    # for chr_name, sub_df in gtf_df.groupby(0,sort=False):
        a=[]
        for gene_name, sub_gene_df in sub_df.groupby(9,sort=False):
            # print(sub_gene_df)
            append_flag = 1
            # 跳过coding和nc之外的基因
            if gene_name not in gene_id_li:
                continue
            # 跳过没有转录本的基因
            if len(sub_gene_df) == 1:
                print(f"{gene_name} has no transcript found !!!")
                num_count["no_tran_num"] += 1
                continue
            # 记录基因信息
            for tran_id, sub_tran_df in sub_gene_df.groupby(10,sort=False):
                # 跳过基因行
                if len(sub_tran_df) == 1:
                    continue
                # 统计转录本类型数目
                tran_prefix = tran_id.split("_")[0]
                num_count[tran_prefix] += 1
                # 选取含NM和NR转录本的gene 创建路径
                if tran_prefix.startswith("N") and append_flag:
                    # a = pd.concat([a,sub_gene_df[[9,3,4,6]].iloc[0]],axis=0)
                    # 处理基因信息
                    gene_item = sub_gene_df[[9, 3, 4, 6]].iloc[0]
                    gene_item[10] = sub_gene_df.iloc[0][8].split('GeneID:')[1].split('"')[0] # ID
                    gene_type_raw = sub_gene_df.iloc[0][8].split('gene_biotype "')[1].split('"')[0]
                    gene_item[11] = "protein coding" if gene_type_raw.find("RNA") == -1 else "non-coding" # type
                    # gene_item = gene_item[[9, 10, 3, 4, 6, 11]]
                    a.append(gene_item)
                    append_flag = 0
                    
                    # 保存exon 和 cds 信息
                    exon_df = sub_tran_df[sub_tran_df[2]=="exon"]
                    cds_df = sub_tran_df[sub_tran_df[2]=="CDS"]
                    Path(f"split_gtf/{chr_name}/{gene_name}").mkdir(exist_ok=1,parents=1)
                    exon_df[[3,4,6]].to_csv(f"split_gtf/{chr_name}/{gene_name}/{tran_id}.exon",header=None,index=False,sep="\t")
                    cds_df[[3,4,6]].to_csv(f"split_gtf/{chr_name}/{gene_name}/{tran_id}.cds",header=None,index=False,sep="\t")
        # 生成只含NM NR 的gene 表
        pd.DataFrame(a).to_csv(f"split_gtf/{chr_name}/Gene_list.tsv",header=None,index=False,sep="\t")
    
    # async_in_iterable_structure(cut_gtf,["NC_000001.11","NC_000002.12"],2)
    # for i in ["NC_000001.11","NC_000002.12"]:
    #     cut_gtf(i)arg
    cut_gtf(argv[1])

    # print(num_count)

    # 字典记录gene->tran->exon->cds 信息

    # 计算程序运行时间
    print(f"time cost:{time.time() - t}")


if __name__ == "__main__":
    main()
