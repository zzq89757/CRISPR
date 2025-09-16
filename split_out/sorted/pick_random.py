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
# 提取行 生成excel
with pd.ExcelWriter('789.xlsx') as writer:
    
    for table in sorted_li:
    # for table in ["spCas9_Homo_chrY.tsv"]:
        chr_name = str(table).split(".")[0].split("_")[-1]
        df = pd.read_csv(table,sep="\t")
        new_header = [x.replace(" ","_")for x in df.columns]
        df.columns = new_header
        # 取行
        random_line = random_line_dict[table]
        df[random_line: random_line + 31].to_excel(writer, sheet_name=chr_name, index=False, header=True)