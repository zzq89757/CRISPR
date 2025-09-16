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

def run_cmd(nc_no: str)->None:
    print(f"awk \'$3~/GG$/\' ag_mark/{nc_no}.tsv > pam_filter/{nc_no}.tsv")
    system(f"awk \'$3~/GG$/\' ag_mark/{nc_no}.tsv > pam_filter/{nc_no}.tsv")


def main() -> None:
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    
    async_in_iterable_structure(run_cmd,nc_li,24)
    

if __name__ == "__main__":
    main()