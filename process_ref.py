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


def obtain_chr_li(ref: str) -> list:
    ref_handle = FastaFile(ref)
    chr_li = ref_handle.references
    ref_handle.close()
    return chr_li


def generate_sgRNA_table(chr_name: str):
    # REFPATH = "/mnt/ntc_data/wayne/DataBase/Homo/hg38.fa"
    REFPATH = "/mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
    ref = FastaFile(REFPATH)
    output = open(
        f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/spCas9_Homo_{chr_name}.tsv", "w"
    )
    print(
        "gRNA Sequence\tPAM\tchr\tchr_strand\tgRNA_start(in chr)\tgRNA_end(in chr)\tgRNA_cut(in chr)",
        file=output,
    )

    chr_seq = ref.fetch(chr_name)
    chr_len = len(chr_seq)
    for i in range(chr_len):
        if chr_seq[i].upper() == "N":
            continue

        # NGG or NAG forward
        if (
            i >= 22
            and chr_seq[i].upper() == "G"
            and (chr_seq[i - 1].upper() == "A" or chr_seq[i - 1].upper() == "G")
        ):
            sgRNA = chr_seq[i - 22 : i - 2].upper()
            if sgRNA.find("TTTT") == -1:
                pam_seq = chr_seq[i - 2 : i + 1].upper()
                grna_start = i - 21
                grna_end = i - 2
                grna_cut = i - 5
                detail = f"{sgRNA}\t{pam_seq}\t{chr_name}\t+\t{grna_start}\t{grna_end}\t{grna_cut}"
                print(detail, file=output)

        # NGG or NAG reverse
        if (
            i <= chr_len - 23
            and chr_seq[i].upper() == "C"
            and (chr_seq[i + 1].upper() == "T" or chr_seq[i + 1].upper() == "C")
        ):
            sgRNA = reverse_complement(chr_seq[i + 3 : i + 23].upper())
            if sgRNA.find("TTTT") == -1:
                pam_seq = reverse_complement(chr_seq[i : i + 3].upper())
                grna_start = i + 23
                grna_end = i + 4
                grna_cut = i + 7
                detail = f"{sgRNA}\t{pam_seq}\t{chr_name}\t-\t{grna_start}\t{grna_end}\t{grna_cut}"
                print(detail, file=output)


def main() -> None:
    # chr_li = obtain_chr_li(REFPATH)
    chr_name_li = list(range(1, 23)) + ["X", "Y"]
    # chr_name_li = ['M']
    chr_li = ["chr" + str(i) for i in chr_name_li]
    async_in_iterable_structure(generate_sgRNA_table, chr_li, 20)
    # generate_sgRNA_table('chr1',ref_handle)


if __name__ == "__main__":
    main()
