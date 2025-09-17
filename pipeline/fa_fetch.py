# 切割fa文件
from pysam import FastaFile
nc_no = "NC_000083.7"

nc_seq = FastaFile("/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/GCF/GCF_000001635.27_GRCm39_genomic.fna").fetch(nc_no)
print(nc_seq)
output_fa = f"./mus/GCF/{nc_no}.fa"

with open(output_fa) as ofa:
    ofa.write(f">{nc_no}\n{nc_seq}\n")

ofa.close() 