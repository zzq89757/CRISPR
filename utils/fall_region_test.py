import pandas as pd
import time

# 读取gene pos 文件和nc2chr文件 记录基因信息
nc_df = pd.read_csv("nc2chr.tsv",sep="\t",header=None)
nc2chr = dict(zip(nc_df[0],nc_df[1]))

gene_df = pd.read_csv("coding_y.tsv",sep="\t",header=None)
# gene_df = gene_df[gene_df[0] == "NC_000001.11"]



def search_ordered(query_number: int) -> None:
    time_st_1 = time.time()
    # 转为 NumPy 数组
    starts = gene_df[1].to_numpy()
    ends = gene_df[2].to_numpy()
    # print(starts)
    if query_number < starts[0] or query_number > ends[-1]:
        print(f"线性查找(sk)...数字{query_number}不在任何区间内。")
        print('s cost:',time.time() - time_st_1)
        return
    mask = (starts <= query_number) & (ends >= query_number)

    if mask.any():
        result = gene_df[mask]
        print(f"线性查找...数字{query_number}在以下区间内：")
        print(result)
    else:
        print(f"线性查找...数字{query_number}不在任何区间内。")
    print('s cost:',time.time() - time_st_1)


def search_sorted(query_number: int) -> None:
    time_st_2 = time.time()
    # 转为 NumPy 数组
    starts = gene_df[1].to_numpy()
    ends = gene_df[2].to_numpy()
    if query_number < starts[0] or query_number > ends[-1]:
        print(f"二分查找(sk)...数字{query_number}不在任何区间内。")
        print('ss cost:',time.time() - time_st_2)
        return
    # 使用 searchsorted 找到可能的起点
    start_idx = starts.searchsorted(query_number, side="right") - 1

    # 向前或向后扩展搜索，找到所有可能的区间
    matches = []
    if start_idx >= 0 and ends[start_idx] >= query_number:
        matches.append(start_idx)

    # print(matches)
    # print(start_idx)

    # 检查周围的重叠区间
    for i in range(start_idx - 1, -1, -1):
        if start_idx - i > 10:
            break
        # print(i)
        if ends[i] >= query_number:
            matches.append(i)

    for i in range(start_idx + 1, len(starts)):
        if starts[i] > query_number or i - start_idx > 10:
            break
        if ends[i] >= query_number:
            matches.append(i)

    # 输出结果
    if matches:
        result = gene_df.iloc[matches]
        print(f"二分查找...数字{query_number}在以下区间内：")
        print(result)
    else:
        print(f"二分查找...数字{query_number}不在任何区间内。")
    print('ss cost:',time.time()-time_st_2)
    

def main() -> None:
    num = 164564
    
    cut_li = pd.read_csv("split_out/sorted/spCas9_Homo_chrY_nohead.tsv",header = None,sep="\t")[6]
    for cut_point in cut_li:
        
        search_ordered(cut_point)
        search_sorted(cut_point)
    

if __name__ == "__main__":
    main()
    
    