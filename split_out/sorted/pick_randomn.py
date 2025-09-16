from pathlib import Path
from random import randint
import pandas as pd


table_li = Path("./").glob("sp*tsv")

# 获取文件行数
line_df = pd.read_csv("line_count.csv",header=None)
line_count_dict = dict(zip(line_df[1],line_df[0]))

# 生成随机数字典
random_line_dict = {table:int(line/2) + randint(-100086,100086) for table, line in line_count_dict.items()}

sorted_li = [f"spCas9_Homo_chr{x}.tsv" for x in list(range(1,23)) + ["X", "Y"]]

new_header = ["gRNA_Sequence", "PAM", "chr", "chr_strand", "gRNA_start(in_chr)", "gRNA_end(in_chr)", "gRNA_cut(in_chr)"]
# 提取行 生成excel
with pd.ExcelWriter('ppp.xlsx') as writer:
    
    for table in sorted_li:
    # for table in ["spCas9_Homo_chrY.tsv"]:
        chr_name = str(table).split(".")[0].split("_")[-1]
        # 取行
        random_line = random_line_dict[table]
        df = pd.read_csv(table,sep="\t",skiprows=random_line,nrows=30,names=new_header,header=None)
        df.to_excel(writer, sheet_name=chr_name, index=False, header=True)
        print(chr_name)