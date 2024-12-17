from pysam import FastxFile
REFPATH = "/mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"


def generate_nc2chr_table(ref_path):
    output_file = open("./nc2chr.tsv",'w')
    for seq in FastxFile(ref_path):
        # skip unlocalized genomic scaffolds and mitochondrion
        if seq.name.startswith("NC") and seq.name != "NC_012920.1":
            print(seq.name,end="\t",file=output_file)
            chr_name = "chr" + seq.comment.split(" ")[3].replace(",","")
            print(chr_name,file=output_file)
            

def main() -> None:
    generate_nc2chr_table(REFPATH)


if __name__ == "__main__":
    main()