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

def run_cmd(i):
    # system(f"awk '$4==\"{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_{i}.tsv")
    print(f"mkdir ./tmp_{i} && java -Xmx120g -jar FlashFry-assembly-1.15.jar  index  --tmpLocation ./tmp_{i}  --database {i}_cas9_db --reference /mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/split_fa/{i}.fa --enzyme spcas9")
    system(f"mkdir ./tmp_{i} && java -Xmx120g -jar FlashFry-assembly-1.15.jar  index  --tmpLocation ./tmp_{i}  --database {i}_cas9_db --reference /mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/split_fa/{i}.fa --enzyme spcas9")
    system(f"java -Xmx8g -jar FlashFry-assembly-1.15.jar  discover  --database NCA_cas9_db  --fasta NC_000001.11.fa  --output 1a_sg.tsv --maxMismatch 3 --maximumOffTargets 100 --forceLinear")


# for i in range(1,23):
#     print(f"awk '$4==\"chr{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_chr{i}.tsv")

# for i in ['X', 'Y']:
#     print(f"awk '$4==\"chr{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_chr{i}.tsv")


def main() -> None:
    nc2chr_file = "../../nc2chr.tsv"
    gtf_file = "GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.gtf"
    # gtf_file = "ZNF568.gtf"
    # num_count = Counter()
    # info_dict = defaultdict()
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    nc_li.remove("NC_000001.11")
    nc_li.remove("NC_000024.10")
    chr_li = ["chr" + str(x) for x in (list(range(1, 23)) + ["X", "Y"])]
    # async_in_iterable_structure(run_cmd,nc_li,22)
    run_cmd("NCA")

if __name__ == "__main__":
    main()