from multiprocessing import Pool, RLock
from os import system

import pandas as pd

def async_in_iterable_structure(fun, iterable_structure, cpus) -> None:
    def init(l):
        global lock
        lock = l

    lock = RLock()
    p = Pool(int(cpus), initializer=init, initargs=(lock,))
    p = Pool(int(cpus), initializer=init, initargs=(lock,))
    # apply async in iterable structure
    for i in iterable_structure:
        p.apply_async(func=fun, args=(i,))
    p.close()
    p.join()

def run_cmd(i) -> None:
    nc2chr_file = "../nc2chr.tsv"
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    chr_li = list(range(1, 23)) + ["X", "Y"]
    print(f"awk '$1==\"{chr_li[i]}\"' ../GCF_000001405.40/00-common_all.vcf > ../vcf_split/raw/{nc_li[i]}.tsv")
    system(f"awk '$1==\"{chr_li[i]}\"' ../GCF_000001405.40/00-common_all.vcf > ../vcf_split/raw/{nc_li[i]}.tsv")
    # system(f"awk '$4==\"{chr_li[i]}\"' ../GCF_000001405.40/*vcf > ../vcf_split/{nc_li[i]}.tsv")
    # print(f"mkdir ./tmp_{i} && java -Xmx8g -jar FlashFry-assembly-1.15.jar  index  --tmpLocation ./tmp_{i}  --database {i}_cas9_db --reference /mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/split_fa/{i}.fa --enzyme spcas9")
    # system(f"mkdir ./tmp_{i} && java -Xmx8g -jar FlashFry-assembly-1.15.jar  index  --tmpLocation ./tmp_{i}  --database {i}_cas9_db --reference /mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/split_fa/{i}.fa --enzyme spcas9")


# for i in range(1,23):
#     print(f"awk '$4==\"chr{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_chr{i}.tsv")

# for i in ['X', 'Y']:
#     print(f"awk '$4==\"chr{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_chr{i}.tsv")


def main() -> None:
    nc2chr_file = "../nc2chr.tsv"
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    chr_li = ["chr" + str(x) for x in (list(range(1, 23)) + ["X", "Y"])]
    async_in_iterable_structure(run_cmd,list(range(24)),22)

if __name__ == "__main__":
    main()