#!/bin/bash
# output delimiter can't be \t and do not support mixed sort
# xsv sort --delimiter '\t' -s "gRNA_start(in chr)","chr_strand" Homo_chrY_KI270740v1_random.tsv
# csvsort -d "\t" -c "gRNA_start(in chr),chr_strand" --reverse "chr_strand" Homo_chrY_KI270740v1_random.csv
# csvsort  -t -k5,4r  ./Homo_chrY_KI270740v1_random.tsv| csvformat -T
# csvsort -t -c 5,6 ./Homo_chrY_KI270740v1_random.tsv | csvformat -T

# LC_ALL=C sort -t$'\t' -k5,5n -k4 ./Homo_chrY_KI270740v1_random.tsv |awk '{if(NR==1) print "gRNA-ID", $0;else print NR-1, $0}' 
# LC_ALL=C sort -t$'\t'  -k5,5n -k4  Homo_chr1.tsv
# LC_ALL=C sort

# scan genome and generate database of single chr
python /mnt/ntc_data/wayne/Repositories/CRISPR/process_ref.py

# final valid sort and merge command

# LC_ALL=C sort -t$'\t' -k5,5n -k4 ./Homo_chrY_KI270740v1_random.tsv |gawk 'NR==1 { print "gRNA-ID", $0} NR>1{ print NR-1, $0}' 

# awk FNR 多文件行号 sort完后多个文件一起加编号（NR去表头 FNR留表头）
# ls Homo*tsv |xargs -I {} echo "LC_ALL=C sort -t \$\'\\\t\' -k5,5n -k4 {}"|xargs -I {} -P 24 bash {}
# ls Homo*tsv |xargs -I {} echo "LC_ALL=C sort -t$\'\\\t\' -k5,5n -k4 {}" |xargs -I {} -P 24 bash {}
# ls Homo*tsv |xargs -I {} echo "LC_ALL=C sort -t$\'\t\' -k5,5n -k4 {} > sorted/{}" > all_chr_sort.sh
# cat all_chr_sort.sh|xargs -I {} -P 24 bash -c {}
# mkdir sorted
ls Homo*tsv |xargs -I {} echo "LC_ALL=C sort -t$'\t' -k5,5n -k4 {} > sorted/{}&" > all_chr_sort.sh
bash all_chr_sort.sh

