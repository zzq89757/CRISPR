from pysam import FastaFile
REFPATH = "/mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"

ref = FastaFile(REFPATH)

chr_seq = ref.fetch('NC_000022.11')

print(chr_seq)