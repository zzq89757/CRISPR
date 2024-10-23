from pysam import FastaFile
from multiprocessing import Pool, RLock


def async_in_iterable_structure(fun, iterable_structure, cpus):
    def init(l):
        global lock
        lock = l

    lock = RLock()
    p = Pool(int(cpus), initializer=init, initargs=(lock,))
    # apply async in iterable structure
    for i in iterable_structure:
        p.apply_async(func=fun, args=[i])
    p.close()
    p.join()


def reverse_complement(seq: str) -> str:
    trantab = str.maketrans("ACGTNacgtnRYMKrymkVBHDvbhd", "TGCANtgcanYRKMyrkmBVDHbvdh")
    return seq.translate(trantab)[::-1]


def generate_sgRNA_table():

    ref = FastaFile("/mnt/ntc_data/wayne/DataBase/Homo/hg38.fa")
    output = open("./sgRNA_sorted.tsv", "w")
    print("no\tchr\tstart\tend\tstrand\tsgRNA", file=output)
    # print(dir(ref))
    # print(ref.get_reference_length('chr1'))
    # print(ref.fetch('chr1',1,10))
    # ref name
    # print(ref.references)

    # chr1_str = ref.fetch('chr1')

    # print(len(chr1_str))
    # exit()
    no = 1

    sort_li = [str(x) for x in range(1, 23)] + ["X", "Y", "M", "Un"]

    chr_li = sorted(
        ref.references, key=lambda x: (sort_li.index(x.split("_")[0][3:]), len(x))
    )

    for chr in chr_li:
        chr_seq = ref.fetch(chr)
        chr_len = len(chr_seq)
        for i in range(chr_len):
            if chr_seq[i] == "N":
                continue

            # NGG or NAG forward
            if (
                i >= 22
                and chr_seq[i].upper() == "G"
                and (chr_seq[i - 1].upper() == "A" or chr_seq[i - 1].upper() == "G")
            ):
                sgRNA = chr_seq[i - 22 : i + 1].upper()
                # print(no, end="\t", file=output)
                print(chr, end="\t", file=output)
                print(i - 21, end="\t", file=output)
                print(i - 1, end="\t", file=output)
                print("+", end="\t", file=output)
                print(sgRNA, end="\n", file=output)
                # print(chr1_str[i - 2 :i + 1],end="\t")
                # print(chr1_str[i - 22 :i + 1],end="\n")
                no += 1

            # NGG or NAG reverse
            if (
                i <= chr_len - 23
                and chr_seq[i].upper() == "C"
                and (chr_seq[i + 1].upper() == "T" or chr_seq[i + 1].upper() == "C")
            ):
                sgRNA = reverse_complement(chr_seq[i : i + 23].upper())
                print(no, end="\t", file=output)
                print(chr, end="\t", file=output)
                print(i + 4, end="\t", file=output)
                print(i + 23, end="\t", file=output)
                print("-", end="\t", file=output)
                print(sgRNA, end="\n", file=output)
                # print(chr_seq[i + 3:i + 23])
                # print(chr1_str[i - 2 :i + 1],end="\t")
                # print(chr1_str[i - 22 :i + 1],end="\n")
                no += 1
                # exit()
