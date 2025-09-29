import pandas as pd

from .filter_intron import scaffold_detective_numpy, gdb2df, region2df


def group_by_gene_and_count(tran_df: pd.DataFrame) -> dict:
    # 统计每个基因包含的转录本数目
    gene_tran_num_counter = dict(tran_df[0].value_counts())
    return gene_tran_num_counter


def tran_count(gdb_df: pd.DataFrame, tran_df: pd.DataFrame, output_file: str) -> None:
    # 只需额外检测没打在外显子的转录本即可 感觉不如双指针
    
    # 统计每个基因包含的转录本数目
    gene_tran_num_counter = group_by_gene_and_count(tran_df)
    
    output_handle = open(output_file,'w')
    exon_df = tran_df
    # 根据gdb 的基因和切点 找到cut的转录本、外显子以及cds并集前2/3位置
    # group_order = pd.Categorical(exon_df[1], categories=exon_df[1].unique(), ordered=True)
    # exon_df = exon_df.sort_values(by=[2, 3], key=lambda col: group_order if col.name == 2 else col)
    # exon_df[4] = exon_df.groupby(1).cumcount() + 1
    exon_df = exon_df.sort_values(2)
    scaffold_pos_li = scaffold_detective_numpy(exon_df,2,3)
    # print(scaffold_pos_li)
    # exon_df.to_csv("exp.tsv",sep="\t",header=None,index=False)
    exon_array = exon_df.to_numpy()
    # 提取所需列并矢量化计算
    gdb_ori_array = gdb_df[4].astype(str) + "0.5"
    gdb_cut_pos_array = gdb_df[7].to_numpy()
    gdb_gene_name_array = gdb_df[8].to_numpy()
    gdb_array = gdb_df.to_numpy()
    
    gene_start_array = exon_df[2].to_numpy()
    gene_end_array = exon_df[3].to_numpy()

    # 计算 g_pos 的矢量化结果
    g_pos_array = gdb_cut_pos_array + gdb_ori_array.astype(float)

    # 筛选 exon 数据
    result = []
    res_count = 0
    idx = 0
    current_gene_pos_idx = 0
    scaffold_pos_len = len(scaffold_pos_li)
    gene_pos_len = len(gene_start_array)
    current_scaffold_pos_idx = 0
    for g_idx, g_pos in enumerate(g_pos_array):
        # idx += 1
        # # 一次性筛选出满足条件的 exon 数据
        # sub_exon_df = exon_df[
        #     (exon_df[0] == gene_name) & (g_pos > exon_df[2]) & (g_pos < exon_df[3])
        # ]
        # print(f"{idx}/{len(gdb_gene_name_array)} finished !!!",end="\r",flush=True)
        # if sub_exon_df.empty:
        #     res_count += 1
        gene_start = gene_start_array[current_gene_pos_idx]
        gene_end = gene_end_array[current_gene_pos_idx]
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
        if current_scaffold_pos_idx < scaffold_pos_len and scaffold_pos_li[current_scaffold_pos_idx][0] <= current_gene_pos_idx <= scaffold_pos_li[current_scaffold_pos_idx][1]:
            tran_li = []
            for j in range(current_gene_pos_idx, scaffold_pos_li[current_scaffold_pos_idx][1] + 1):
                # if gene_start_array[j] < g_pos < gene_end_array[j]:
                if gene_start_array[j] < g_pos < gene_end_array[j]:
                    if gdb_array[g_idx][8] == exon_array[j][0]:
                        # print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                        # print(exon_array[j][3],end="\n",file=output_handle)
                        tran_li.append(exon_array[j][1])
                    
                    # print(g_pos,end="\t")
                        ...
                    # print(gene_start_array[j],end="\t")
                    # print(gene_end_array[j],end="\t")
                    # print(gene_pos_array[j])
            if tran_li:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                print(f"{len(tran_li)}/{gene_tran_num_counter[gdb_array[g_idx][8]]}",end="\n",file=output_handle)
            continue
        
        
        # 如果 gRNA 位于当前基因范围内，记录信息 chr1:1335323 跨两个或以上的位点会被跳过
        if gene_start < g_pos < gene_end:
            # print(f"{g_pos} in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}", flush=True, end="\r")
            # print(f"{g_pos} in {gene_start}-{gene_end}")
            # g_item = gdb.iloc[g_idx]
            # g_item = gdb_array[g_idx]
            if gdb_array[g_idx][8] == exon_array[current_gene_pos_idx][0]:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                print("1/1",end="\n",file=output_handle)
            #     print(gdb_array[g_idx][8],end="\t")
            #     print(exon_array[current_gene_pos_idx][0],end="\t")
            #     print(g_pos)
            # ...
            # print(g_item)
        else:
            # print(f"{g_pos} not in {gene_start}-{gene_end} and {gene_start_array[current_gene_pos_idx - 1]}-{gene_end_array[current_gene_pos_idx - 1]}")
            # g_item = gdb_array[g_idx]
            # print(g_item)
            ...
        
        
    return 
    
    


def tran_ratio_count(project_dir: str, nc_no: str) -> None:       
    # 读取gdb
    gdb_file = f"{project_dir}/cds_mark/{nc_no}.tsv"
    gdb_df = gdb2df(gdb_file)
    # 读取对应的tran表
    tran_file = f"{project_dir}/GCF/gtf/{nc_no}/TRAN.tsv"
    output_file = f"{project_dir}/tran_count/{nc_no}.tsv"
    tran_df = region2df(tran_file)
    # 注释gdb
    tran_count(gdb_df, tran_df, output_file)
    


def main() -> None:
    nc2chr_file = "/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/nc_li"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    tran_ratio_count("/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus",nc_li[-1])
    # run_count("NC_000024.10")
    

if __name__ == "__main__":
    main()