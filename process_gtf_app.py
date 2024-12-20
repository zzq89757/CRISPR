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
    for i in iterable_structure:
        p.apply_async(func=fun, args=(i,))
    p.close()
    p.join()


def gtf2df(gtf_file: str, nc_no: str) -> pd.DataFrame:
    """提前预设每列的数据类型并将gtf文件存入DataFrame"""
    use_col_li = [0, 2, 3, 4, 6, 8]
    type_li = ["category", "category", "int32", "int32", "category", "string"]

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

    return gtf_df[gtf_df[0]==nc_no]


def append_gene_id_col(gtf_df: pd.DataFrame) -> None:
    """为gtf添加gene id 和 tran id 列"""
    gtf_df[[10 ,9]] = gtf_df[8].str.extract(
        r'transcript_id "([^"]*)".*gene "([^"]*)";'
    )
    gtf_df[9] = gtf_df[9].fillna('')
    gtf_df[10] = gtf_df[10].fillna('')


def obtain_gene_id_li(gtf_df: pd.DataFrame) -> defaultdict:
    """挑出gene行,只保留gene biotype 为protein_coding和ncRNA的gene id 列表并统计数目"""
    gene_df = gtf_df[gtf_df[2] == "gene"]
    contains_bool_li1 = gene_df[8].str.contains(r"gene_biotype \"protein_coding\";")
    contains_bool_li2 = gene_df[8].str.contains(r"gene_biotype \"ncRNA\";")
    contains_bool_li3 = gene_df[8].str.contains(r"gene_biotype \"lncRNA\";")
    contains_bool_li = [x or y or z for x, y, z in zip(contains_bool_li1, contains_bool_li2, contains_bool_li3)]
    gene_df = gene_df[contains_bool_li]
    gene_id_li = gene_df[9].to_list()
    return gene_id_li


def cut_gtf(nc_no, gtf_df, gene_id_li):
    sub_df = gtf_df[gtf_df[0]==nc_no]
    chr_name = nc_no
    Path(f"split_gtf/extract/{chr_name}").mkdir(exist_ok=1,parents=1)
    all_gene_li=[]
    tran_exon_li = pd.DataFrame([])
    tran_cds_li = pd.DataFrame([])
    gene_tran_li = pd.DataFrame([])
    for gene_name, sub_gene_df in sub_df.groupby(9,sort=False):
        append_flag = 1
        # 跳过coding和nc之外的基因
        if gene_name not in gene_id_li:
            continue
        # 跳过没有转录本的基因
        if len(sub_gene_df) == 1:
            print(f"{gene_name} has no transcript found !!!")
            continue
        # 记录基因信息
        for tran_id, sub_tran_df in sub_gene_df.groupby(10,sort=False):
            # print(sub_tran_df)
            # 跳过基因行(基因行的转录本 id 为空 分组后仅一行)
            if len(sub_tran_df) == 1:
                continue
            # 统计转录本类型数目
            tran_prefix = tran_id.split("_")[0]
            # 选取含NM和NR转录本的gene 创建路径
            if tran_prefix.startswith("N"): #  此处可以在添加基因和转录本信息后直接过滤
                # 处理基因信息
                if append_flag:
                    gene_item = sub_gene_df[[9, 3, 4, 6]].iloc[0]
                    gene_item[10] = sub_gene_df.iloc[0][8].split('GeneID:')[1].split('"')[0] # ID
                    gene_type_raw = sub_gene_df.iloc[0][8].split('gene_biotype "')[1].split('"')[0]
                    if gene_type_raw not in ["protein_coding", "lncRNA", "ncRNA"]:
                        print(sub_gene_df)
                    gene_item[11] = "protein_coding" if gene_type_raw.find("RNA") == -1 else "non_coding" # type
                    # gene_item = gene_item[[9, 10, 3, 4, 6, 11]]
                    all_gene_li.append(gene_item)
                    append_flag = 0
                
                # 保存exon 和 cds 信息
                exon_df = sub_tran_df[sub_tran_df[2]=="exon"][[9, 10, 3, 4]]
                cds_df = sub_tran_df[sub_tran_df[2]=="CDS"][[9, 10, 3, 4]]
                tran_cds_li = pd.concat([tran_cds_li, cds_df])
                tran_exon_li = pd.concat([tran_exon_li, exon_df])
    # 生成只含NM NR 的 gene表 cds表 和 exon表
    pd.DataFrame(all_gene_li).to_csv(f"split_gtf/extract/{chr_name}/Gene_list_re.tsv",header=None,index=False,sep="\t")
    pd.DataFrame(tran_cds_li).to_csv(f"split_gtf/extract/{chr_name}/CDS_re.tsv",header=None,index=False,sep="\t")
    pd.DataFrame(tran_exon_li).to_csv(f"split_gtf/extract/{chr_name}/EXON_re.tsv",header=None,index=False,sep="\t")


def run_cut(nc_no) -> None:
    gtf_file = f"split_gtf/raw_gtf_split/{nc_no}.gtf"
    # 记录程序开始时间
    t = time.time()
    # 读取gtf文件 生成df
    gtf_df = gtf2df(gtf_file, nc_no)
    # 添加gene id 和tran id列
    append_gene_id_col(gtf_df)
    # 找到coding和uncoding的gene id 列表
    gene_id_li = obtain_gene_id_li(gtf_df)
    # 按照gene id group by 只保留gene id 列表中且包含NM和NR的转录本（统计数目）  
    cut_gtf(nc_no,gtf_df,gene_id_li)
    # 计算程序运行时间
    print(f"time cost:{time.time() - t}")


def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    chr_li = ["chr" + str(x) for x in (list(range(1, 23)) + ["X", "Y"])]
    # async_in_iterable_structure(run_cut,nc_li,24)
    run_cut(chr_li[-1])

if __name__ == "__main__":
    main()
