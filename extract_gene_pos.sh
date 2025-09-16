# awk '$3=="gene"{print $1"\t"$4"\t"$5"\t"$7"\t"$10}' head50.gtf
awk '$1~/NC_0000/&&$3=="gene"&&$0!~/pseudogene/{print $1"\t"$4"\t"$5"\t"$7"\t"$10}' GCF_000001405.40/genomic.gtf > gene_pos.tsv
awk '$1~/NC_0000/&&$3=="gene"&&$0!~/pseudogene/&&$0!~/uncharacterized/' GCF_000001405.40/genomic.gtf > gene_pos.tsv

awk '$3=="gene"{print $NF}' GCF_000001405.40/genomic.gtf
awk -F 'gene_biotype ' '{if ($2&&$0~/\tgene\t/) {match($2, /"[^"]+"/, a); print a[0]}}' GCF_000001405.40/genomic.gtf
awk '$3=="gene"&&($0~/ncRNA/||$0~/protein_coding/){print $NF}' GCF_000001405.40/genomic.gtf|wc -l 

awk '$3=="gene"&&($0~/protein_coding/||$0~/ncRNA/)&&$0~/^NC_/{print $1"\t"$10"\t"$14}' *gtf > coding.tsv