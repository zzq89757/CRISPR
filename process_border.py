import time
import pandas as pd
from json import dump


# 找出重叠区域的起始终止索引
def scaffold_detective(gene_pos_df: pd.DataFrame) -> list:
    common_idx_li = []
    terminal_common_idx_li = []
    for i in range(len(gene_pos_df) - 1):
        # 跳过交集区
        if i in common_idx_li:
            continue
        # 无重叠 直接append
        if gene_pos_df.iloc[i + 1][1] > gene_pos_df.iloc[i][2]:
            continue
        else:
            # 记录MAX end 并将后续区间同maxend比较以判断是否为连续交集区
            max_end = gene_pos_df.iloc[i][2]
            min_j = i
            max_j = i
            for j in range(i, len(gene_pos_df) - 1):
                if gene_pos_df.iloc[j][2] > max_end:
                    max_end = gene_pos_df.iloc[j][2]
                # tmp_res.append(df.iloc[j])
                common_idx_li.append(j)
                if j > max_j:
                    max_j = j
                # print(f"max_end is {max_end}")
                # 下一区间的起始大于maxend 退出循环
                if gene_pos_df.iloc[j + 1][1] > max_end:
                    # print(j)
                    break
                # 判断最后一行是否可并入连续交集区
                if j == len(gene_pos_df) - 2:
                    if gene_pos_df.iloc[j][1] < max_end:
                        # tmp_res.append(df.iloc[j + 1])
                        common_idx_li.append(j + 1)
                        if j + 1 > max_j:
                            max_j = j + 1
            terminal_common_idx_li.append([min_j, max_j])
    # print(common_idx_li)
    # print(terminal_common_idx_li)
    # json_handle = open("sb.json",'w')
    # json_str = dump({'ScaffoldBorder' : terminal_common_idx_li}, json_handle)
    # print(json_str)
    # return pd.DataFrame(res)
    return terminal_common_idx_li

# 找出重叠区域的起始终止索引
def scaffold_detective_numpy(gene_pos_df: pd.DataFrame) -> list:
    gene_start_array = gene_pos_df[1].to_numpy()
    gene_end_array = gene_pos_df[2].to_numpy()
    common_idx_li = []
    terminal_common_idx_li = []
    for i in range(len(gene_start_array) - 1):
        # 跳过交集区
        if i in common_idx_li:
            continue
        # 无重叠 直接append
        if gene_start_array[i + 1] > gene_end_array[i]:
            continue
        else:
            # 记录MAX end 并将后续区间同maxend比较以判断是否为连续交集区
            max_end = gene_end_array[i]
            min_j = i
            max_j = i
            for j in range(i, len(gene_pos_df) - 1):
                if gene_end_array[j] > max_end:
                    max_end = gene_end_array[j]
                # tmp_res.append(df.iloc[j])
                common_idx_li.append(j)
                if j > max_j:
                    max_j = j
                # print(f"max_end is {max_end}")
                # 下一区间的起始大于maxend 退出循环
                if gene_start_array[j + 1] > max_end:
                    # print(j)
                    break
                # 判断最后一行是否可并入连续交集区
                if j == len(gene_pos_df) - 2:
                    if gene_start_array[j] < max_end:
                        # tmp_res.append(df.iloc[j + 1])
                        common_idx_li.append(j + 1)
                        if j + 1 > max_j:
                            max_j = j + 1
            terminal_common_idx_li.append([min_j, max_j])
    return terminal_common_idx_li


def relative_pos_calc_1(gr_start: int, gr_end: int, gr_cut: int, gr_ori: str, gene_start: int, gene_end: int, gene_ori: str) -> list:
    real_gene_start = gene_start if gene_ori == "+" else gene_end
    relative_start = abs(gr_start - real_gene_start) + 1
    relative_end = abs(gr_end - real_gene_start) + 1
    relative_cut = abs(gr_cut - real_gene_start) + 1
    right_softclip = ''
    left_softclip = ''
    # 若相对起止 超出基因范围
    gene_len = gene_end - gene_start + 1
    g_left_terminal = min(gr_start,gr_end)
    g_right_terminal = max(gr_start,gr_end)

    if g_right_terminal > gene_end:
        # gr end 超出gene右端
        right_softclip = f"(+{g_right_terminal - gene_end}bp)"
        if gr_ori == "+":
            relative_end = gene_len if gene_ori == "+" else 1
        else:
            relative_start = gene_len if gene_ori == "+" else 1

    
    if g_left_terminal < gene_start:
        left_softclip = f"(-{gene_start - g_left_terminal}bp)"
        if gr_ori == "+":
            relative_start = 1 if gene_ori == "+" else gene_len
        else:
            relative_end = 1 if gene_ori == "+" else gene_len
    
    relative_loc = f"{left_softclip}{relative_start}-{relative_end}{right_softclip}" if gr_ori == "+" else f"{right_softclip}{relative_start}-{relative_end}{left_softclip}"
    return [relative_loc, relative_cut]

def relative_pos_calc_2(gr_start: int, gr_end: int, gr_cut: int, gr_ori: str, gene_start: int, gene_end: int, gene_ori: str) -> list:
    real_gene_start = gene_start if gene_ori == "+" else gene_end
    relative_start = abs(gr_start - real_gene_start) + 1
    relative_end = abs(gr_end - real_gene_start) + 1
    relative_cut = abs(gr_cut - real_gene_start) + 1
    right_softclip = ''
    left_softclip = ''
    relative_start_softclip = ''
    relative_end_softclip = ''
    # 若相对起止 超出基因范围
    gene_len = gene_end - gene_start + 1
    g_left_terminal = min(gr_start,gr_end)
    g_right_terminal = max(gr_start,gr_end)

    
    relative_start = gr_start - real_gene_start + 1 
    if relative_start > gene_len:
        relative_start = gene_len
        relative_start_softclip = f"+{relative_start - gene_len}"
    if relative_start < 1:
        relative_start -= 1
    relative_end = gr_end - real_gene_start + 1
    if relative_end > gene_len:
        relative_end = gene_len
        relative_end_softclip = f"+{relative_end - gene_len}"
    if relative_end < 1:
        relative_end -= 1
    loc_str = f"{relative_start_softclip}{relative_start}-{relative_end}{relative_end_softclip}"
    # print(loc_str)
    return [loc_str, relative_cut]


def relative_pos_calc(gr_start: int, gr_end: int, gr_cut: int, gr_ori: str, gene_terminal_left: int, gene_terminal_right: int, gene_ori: str) -> list:
    # gr start > gr end while gr ori == -,gene_terminal_left eternal less than gene_terminal_right 
    loc_str = ''
    relative_cut = 0
    relative_start = 0
    relative_end = 0
    gene_len = gene_terminal_right - gene_terminal_left + 1
    real_gene_start = gene_terminal_left if gene_ori == "+" else gene_terminal_right
    relative_cut = abs(gr_cut - real_gene_start) + 1
    
    # 根据gr位置判断gRNA是否超出基因范围
    gr_left = min(gr_start, gr_end)
    gr_right = max(gr_start, gr_end)
    if  gr_left >= gene_terminal_left and gr_right <= gene_terminal_right:
        relative_start = abs(gr_start - real_gene_start) + 1
        relative_end = abs(gr_end - real_gene_start) + 1
        loc_str = f"{relative_start}-{relative_end}"
        return [loc_str, relative_cut]
    
    # ++ and --
    if gr_ori == gene_ori:
        if gr_ori == "+": # ++
            if gr_end > gene_terminal_right:
                relative_start = gr_start - gene_terminal_left + 1
                relative_end = gene_len
                
                loc_str = f"{relative_start}-{relative_end}(+{gr_end - gene_terminal_right})"
                return [loc_str, relative_cut]
            if gr_start < gene_terminal_left:
                relative_start = gr_start - gene_terminal_left
                relative_end = gr_end - gene_terminal_left + 1
        
                loc_str = f"({relative_start})-{relative_end}"
                return [loc_str, relative_cut]
        else: # --
            if gr_start > gene_terminal_right:
                relative_end = gene_terminal_right - gr_end + 1
                relative_start = gene_terminal_right - gr_start
        
                loc_str = f"({relative_start})-{relative_end}"
                return [loc_str, relative_cut]
                
            if gr_end < gene_terminal_left:
                relative_start = gene_terminal_right - gr_start + 1
                relative_end = gene_len
                
                loc_str = f"{relative_start}-{relative_end}(+{gene_terminal_left - gr_end})"
                return [loc_str, relative_cut]
    else: # +- and -+
        if gr_ori == "+": # +-
            if gr_end > gene_terminal_right:
                relative_start = gene_terminal_right - gr_start + 1
                relative_end = gene_terminal_right - gr_end
                
                loc_str = f"{relative_start}-({relative_end})"
                return [loc_str, relative_cut]
            if gr_start < gene_terminal_left:
                relative_start = gene_len
                relative_end = gene_terminal_right - gr_end + 1
                
                loc_str = f"(+{gene_terminal_left - gr_start}){relative_start}-{relative_end}"
                return [loc_str, relative_cut]
        else: # -+
            if gr_start > gene_terminal_right:
                relative_end = gene_terminal_left - gr_start + 1
                relative_start = gene_len
                
                loc_str = f"(+{gr_start - gene_terminal_right}){relative_start}-{relative_end}"
                return [loc_str, relative_cut]
            if gr_end < gene_terminal_left:
                relative_start = gr_start - gene_terminal_left + 1
                relative_end = gr_end - gene_terminal_left
            
                loc_str = f"{relative_start}-({relative_end})"
                return [loc_str, relative_cut]
    return ["not consider", relative_cut]
        
    


def gdb_annotation(gdb: pd.DataFrame, gene_pos_df: pd.DataFrame, scaffold_pos_li: list) -> pd.DataFrame:
    # 需要分别考虑切点和起止 三个位置 超出基因末端时 loc 采用(+/- len(over))
    # 提取必要列为 NumPy 数组
    gene_start_array = gene_pos_df[1].to_numpy()
    gene_end_array = gene_pos_df[2].to_numpy()
    gRNA_ori_array = gdb[3].to_numpy()
    gRNA_pos_array = gdb[6].to_numpy()
    gdb_array = gdb.to_numpy()
    gene_pos_array = gene_pos_df.to_numpy()
    gene_pos_len = len(gene_start_array)
    scaffold_pos_len = len(scaffold_pos_li)
    current_gene_pos_idx = 0
    current_scaffold_pos_idx = 0
    # max_end = max(gene_end_array)
    # min_start = min(gene_end_array)
    ## 遍历 gRNA <考虑gRNA不同方向时切点位置不同>
    for g_idx, g_raw_pos in enumerate(gRNA_pos_array):
        g_ori = gRNA_ori_array[g_idx]
        g_pos_offset = float(g_ori + "0.5")
        g_pos = g_raw_pos + g_pos_offset
        # 不分类考虑正负链的情况 后续过滤掉即可
        gene_start = gene_start_array[current_gene_pos_idx]
        gene_end = gene_end_array[current_gene_pos_idx]
        # 跳过位于当前基因起始位置之前以及基因间的的 gRNA
        if g_pos < gene_start:
            continue            
        # 跨多个基因时的情况
        while (current_gene_pos_idx < gene_pos_len -1  and g_pos > gene_end):
            current_gene_pos_idx += 1
            # 同步更新scaffold_pos_idx
            border_right = scaffold_pos_li[current_scaffold_pos_idx][1] if current_scaffold_pos_idx < scaffold_pos_len else scaffold_pos_li[-1][1]
            if current_gene_pos_idx > border_right and current_scaffold_pos_idx <= scaffold_pos_len:
                current_scaffold_pos_idx += 1
            gene_start = gene_start_array[current_gene_pos_idx]
            gene_end = gene_end_array[current_gene_pos_idx]
        
        # 大于max end时 结束
        if current_gene_pos_idx == gene_pos_len -1 and g_pos > gene_end_array[-1]:
            break
        
        # 处于scaffold中的gRNA处理
        if current_scaffold_pos_idx < scaffold_pos_len and current_scaffold_pos_idx < scaffold_pos_li[current_scaffold_pos_idx][0] <= current_gene_pos_idx <= scaffold_pos_li[current_scaffold_pos_idx][1]:
            for j in range(current_gene_pos_idx, scaffold_pos_li[current_scaffold_pos_idx][1] + 1):
                if gene_start_array[j] < g_pos < gene_end_array[j]:
                    print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t")
                    print(gene_pos_array[j][0],end="\t") # Gene Symbol
                    print(gene_pos_array[j][4],end="\t") # Gene ID from NCBI
                    print(gene_pos_array[j][5],end="\t") # gene type
                    relative_ori = 'fwd' if g_ori == gene_pos_array[j][3] else 'rev'
                    print(relative_ori,end="\t") # relative ori
                    # # 计算gRNA relative loc和relative cut
                    relative_loc, relative_cut = relative_pos_calc(gdb_array[g_idx][4], gdb_array[g_idx][5], g_raw_pos, g_ori, gene_start_array[j], gene_end_array[j], gene_pos_array[j][3])
                    print(relative_loc,end="\t") # relative loc
                    print(relative_cut,end="\n") # relative cut
                    # print(g_pos,end="\t")
                    ...
                    # print(gene_start_array[j],end="\t")
                    # print(gene_end_array[j],end="\t")
                    # print(gene_pos_array[j])
            continue
        
        
        # 如果 gRNA 位于当前基因范围内，记录信息 chr1:1335323 跨两个或以上的位点会被跳过
        if gene_start < g_pos < gene_end:
            # print(f"{g_pos} in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}", flush=True, end="\r")
            # print(f"{g_pos} in {gene_start}-{gene_end}")
            # g_item = gdb.iloc[g_idx]
            # g_item = gdb_array[g_idx]
            ...
            # print(g_item)
        else:
            # print(f"{g_pos} not in {gene_start}-{gene_end} and {gene_start_array[current_gene_pos_idx - 1]}-{gene_end_array[current_gene_pos_idx - 1]}")
            # g_item = gdb_array[g_idx]
            # print(g_item)
            ...
    
    return pd.DataFrame([])
    




def main() -> None:
    t1 = time.time()
    nc_no = "NC_000024.10"
    chr = "chrY"

    nc_table = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    df = pd.read_csv(nc_table, sep="\t", header=None)
    chr2nc_dict = df.set_index(1)[0].to_dict()
    type_li = ["string", "int32", "int32", "category", "string", "category"]
    type_dict = dict(enumerate(type_li))
    # 读取基因位置信息文件
    gene_pos_df = pd.read_csv(
        f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{chr2nc_dict[chr]}/Gene_list.tsv",
        sep="\t",
        header=None,
        # low_memory=False,
        dtype=type_dict,
    )
    # 读取gRNA db 文件
    type_li = ["string", "string", "category", "category", "int32", "int32", "int32"]
    type_dict = dict(enumerate(type_li))
    gdb_df = pd.read_csv(
        f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/spCas9_Homo_{chr}.tsv",
        header=None,
        sep="\t",
        # low_memory=False,
        dtype=type_dict
    )
    # 找出重叠区域的起始终止索引
    # t2 = time.time()
    # sccaffold_pos_li = scaffold_detective(gene_pos_df)
    # print(f"time cost:{time.time() - t2}")
    # print(sccaffold_pos_li)
    print(f"load gdb, time cost:{time.time() - t1}")
    t2 = time.time()
    sccaffold_pos_li = scaffold_detective_numpy(gene_pos_df)
    print(f"find insertion, time cost:{time.time() - t2}")
    # print(sccaffold_pos_li)
    # 同时遍历gdb和gtf 为gdb条目添加基因注释
    gdb_annotation(gdb_df, gene_pos_df, sccaffold_pos_li)
    
    print(f"annotation,time cost:{time.time() - t2}")
    
    
    
if __name__ == "__main__":
    main()