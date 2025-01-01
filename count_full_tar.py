import numpy as np
import pandas as pd
from utils.read_tsv import tsv2df


# 自定义转换函数
def fraction_to_int(value):
    try:
        if '/' in value:
            numerator, denominator = map(int, value.split('/'))
            return numerator // denominator  # 取整除
        return int(value)  # 如果是整数，直接转换
    except:
        return None  # 处理错误情况
    
    
    
def run_count(nc_no: str) -> None:
    # 读取gdb
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/tran_count/{nc_no}.tsv"
    type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category"]
    
    gdb_df = tsv2df(gdb_path, type_li)
    # gdb_df[16] = gdb_df[16].apply(lambda x:int(eval(x)))
    fraction_to_int_vec = np.vectorize(fraction_to_int)
    int_arr = fraction_to_int_vec(gdb_df[16].to_numpy())
    # print(int_arr)
    gdb_df[16] = int_arr
    gdb_df= gdb_df[(gdb_df[15]=="yes") & (gdb_df[16]==1)]
    print(f"{nc_no}\t{len(gdb_df)}")

def main() -> None:
    # 读取nc2chr_file 生成 NC -> chr 的映射字典
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    [run_count(x) for x in nc_li]
    # run_count(nc_li[0])
    
    
if __name__ == "__main__":
    main()