from Bio import SeqIO

# 输入文件和目标范围
fasta_file = "GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
target_chromosome = "NC_000003.12"
start = 129529090
end = 129529109

# 提取目标区域序列
with open(fasta_file, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        if record.id == target_chromosome:  # 找到目标染色体
            target_sequence = record.seq[start - 1:end]  # 提取目标区域（注意 Python 的索引从 0 开始）
            print(f">Sequence_{target_chromosome}_{start}_{end}")
            print(target_sequence)
            break
