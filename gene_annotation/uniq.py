import pandas as pd
from collections import Counter

# 读取两个文件
file1_df = pd.read_csv("7_1_re.tsv", header=None)[0]
file2_df = pd.read_csv("7_1.tsv", header=None)[0]
print(Counter(file1_df)-Counter(file2_df))
# print(file2_df)
exit()
# 取并集和交集
union_ids = pd.concat([file1_df, file2_df]).drop_duplicates()
intersection_ids = pd.merge(file1_df, file2_df, on="ID")

# 找到只在一个文件中的 ID
unique_ids = pd.concat([union_ids, intersection_ids]).drop_duplicates(keep=False)

# 输出结果
print(unique_ids)