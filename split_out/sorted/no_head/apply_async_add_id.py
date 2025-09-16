from multiprocessing import Pool, RLock
from os import system

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

def run_cmd(i):
    print(f"awk '$4==\"{i}\"' 'spCas9_Homo(WGS)_gRNA-Original.tsv' > with_id/spCas9_Homo_{i}.tsv")
    system(f"awk '$4==\"{i}\"' 'spCas9_Homo(WGS)_gRNA-Original.tsv' > with_id/spCas9_Homo_{i}.tsv")


# for i in range(1,23):
#     print(f"awk '$4==\"chr{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_chr{i}.tsv")

# for i in ['X', 'Y']:
#     print(f"awk '$4==\"chr{i}\"' split_out/sorted/'spCas9_Homo(WGS)_gRNA-Original.tsv' > ori_split/ori_chr{i}.tsv")


def main() -> None:
    chr_li = ["chr" + str(x) for x in (list(range(1, 23)) + ["X", "Y"])]
    async_in_iterable_structure(run_cmd,chr_li,12)

if __name__ == "__main__":
    main()