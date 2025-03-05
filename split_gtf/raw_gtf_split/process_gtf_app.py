import time
import pandas as pd
from collections import defaultdict
from pathlib import Path
from multiprocessing import Pool, RLock

def async_in_iterable_structure(fun, iterable_structure, cpus) -> None:
    def init(l) -> None:
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


def cut_gtf(nc_no: str, gtf_df: pd.DataFrame) -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    chr_li = ["chr" + str(x) for x in (list(range(1, 23)) + ["X", "Y"])]
    nc_chr_dict = dict(zip(nc_li, chr_li))
    gtf_df = gtf_df[(gtf_df[2]=="gene")|(gtf_df[2]=="transcript")]
    chr_name = nc_no
    all_gene_li=[]
    for gene_name, sub_gene_df in gtf_df.groupby(9,sort=False):
        # 记录基因信息
        gene_item = sub_gene_df[[0, 9, 3, 4, 6]].iloc[0]
        # 提取GeneID和基因类型信息
        gene_item[10] = sub_gene_df.iloc[0][8].split('GeneID:')[1].split('"')[0] # ID
        gene_item[11] = sub_gene_df.iloc[0][8].split('gene_biotype "')[1].split('"')[0]
        # 添加chr列
        gene_item[12] = nc_chr_dict[sub_gene_df.iloc[0][0]]
        # 获取基因对应转录本列表
        tran_li = sub_gene_df[10].to_list()[1:]
        # 拼接转录本字符串，记录转录本数目
        tran_str = ",".join(tran_li)
        tran_num = len(tran_li)
        gene_item[13] = tran_str
        gene_item[14] = tran_num
        # 重新对列进行排序
        all_gene_li.append(gene_item[[9, 10, 11, 0, 12, 6, 3, 4, 13, 14]])
    # 生成gene info表
    pd.DataFrame(all_gene_li).to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/raw_gtf_split/{chr_name}.tsv",header=None,index=False,sep="\t")



def run_cut(nc_no: str) -> None:
    gtf_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/raw_gtf_split/{nc_no}.gtf"
    # 记录程序开始时间
    t = time.time()
    # 读取gtf文件 生成df
    gtf_df = gtf2df(gtf_file, nc_no)
    # 添加gene id 和tran id列
    append_gene_id_col(gtf_df)
    # 按照gene id group by 只保留gene id 列表中且包含NM和NR的转录本（统计数目）  
    cut_gtf(nc_no,gtf_df)
    # 计算程序运行时间
    print(f"<{nc_no}> time cost:{time.time() - t}")


def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_cut,nc_li,24)
    # run_cut(nc_li[-1])

if __name__ == "__main__":
    main()
