from pathlib import Path
import pandas as pd
import re
from pysam import FastaFile


def reverse_complement(seq: str) -> str:
    trantab = str.maketrans("ACGTNacgtnRYMKrymkVBHDvbhd", "TGCANtgcanYRKMyrkmBVDHbvdh")
    return seq.translate(trantab)[::-1]

def ref_scan(nc_no: str, chr_name:str, ref_file: str, raw_db: str) -> None:
    ref = FastaFile(ref_file)
    output = open(
        raw_db, "w"
    )
    print(
        "gRNA_Seq\tPAM\tchr\tchr_strand\tgRNA_start(in chr)\tgRNA_end(in chr)\tgRNA_cut(in chr)",
        file=output,
    )

    ref_seq = ref.fetch(nc_no)
    ref_len = len(ref_seq)
    for i in range(ref_len):
        if ref_seq[i] == "N":
            continue

        # NGG or NAG forward
        if (
            i >= 22
            and ref_seq[i] == "G"
            and (ref_seq[i - 1] == "A" or ref_seq[i - 1] == "G")
        ):
            sgRNA = ref_seq[i - 22 : i - 2]
            # if sgRNA.find("TTTT") == -1 and ref_seq[i - 22 : i - 1].find("N") == -1:
            if sgRNA.find("TTTT") == -1 and (not re.search(r"[RYMKVBHDN]",ref_seq[i - 22 : i - 1])):
                pam_seq = ref_seq[i - 2 : i + 1]
                grna_start = i - 21
                grna_end = i - 2
                grna_cut = i - 5
                detail = f"{sgRNA}\t{pam_seq}\t{chr_name}\t+\t{grna_start}\t{grna_end}\t{grna_cut}"
                print(detail, file=output)

        # NGG or NAG reverse
        if (
            i <= ref_len - 23
            and ref_seq[i] == "C"
            and (ref_seq[i + 1] == "T" or ref_seq[i + 1] == "C")
        ):
            sgRNA = reverse_complement(ref_seq[i + 3 : i + 23])
            # if sgRNA.find("TTTT") == -1 and ref_seq[i + 2 : i + 23].find("N") == -1:
            if sgRNA.find("TTTT") == -1 and (not re.search(r"[RYMKVBHDN]",ref_seq[i + 2 : i + 23])):
                pam_seq = reverse_complement(ref_seq[i : i + 3])
                grna_start = i + 23
                grna_end = i + 4
                grna_cut = i + 7
                detail = f"{sgRNA}\t{pam_seq}\t{chr_name}\t-\t{grna_start}\t{grna_end}\t{grna_cut}"
                print(detail, file=output)


if __name__ == "__main__":
    ref_scan("NC_000075.7", "chr9", "/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/GCF/fa/NC_000075.7.fa", "/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/ref_scan/NC_000075.7.tsv")