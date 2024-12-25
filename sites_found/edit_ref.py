import pysam

# print(dir(pysam.FastxFile("/mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna")))
# exit()
for chr in pysam.FastxFile("/mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"):
    print(chr)
    
    