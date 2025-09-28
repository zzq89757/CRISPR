import pandas as pd
from .filter_intron import scaffold_detective_numpy, gdb2df, region2df


def common_cds_front_region(project_dir: str, cds_df: pd.DataFrame, nc_no: str) -> pd.DataFrame:
    gene_df: pd.DataFrame = pd.read_csv(f"{project_dir}/GCF/gtf/{nc_no}/GENE.tsv",sep="\t",header=None)

    gene_ori_dict = dict(zip(gene_df[0], gene_df[3]))
    # 找到所有转录本对应CDS的前2/3后取并集 注意方向!!!
    new_res = pd.DataFrame([])
    cds_df.columns = ["Gene", "Transcript", "CDS_start", "CDS_end"]
    # 按照转录本分组
    for tran, tran_df in cds_df.groupby("Transcript"):
        ori = gene_ori_dict[tran_df["Gene"].to_numpy()[0]]
        result = []
        # 计算总长
        full_cds_len = sum(tran_df["CDS_end"] - tran_df["CDS_start"])
        # 找出前2/3 full_cds_len - 每行的区间长度 当剩余数字为负时 将终点改为该行终点加上剩余数字 丢弃后续cds区间
        target_length = full_cds_len * 2 / 3  # 目标长度（前2/3）

        accumulated_length = 0  # 累积长度
        for idx, row in tran_df.iterrows():
            start, end = row["CDS_start"], row["CDS_end"]
            interval_length = end - start

            # 判断是否超过目标长度
            if accumulated_length + interval_length > target_length:
                # 修正终点，取部分区间
                remaining_length = target_length - accumulated_length
                if ori == '+':
                    new_end = int(start + remaining_length)
                    result.append([row["Gene"], tran, start, new_end])
                else:
                    new_start = int(end - remaining_length)
                    result.append([row["Gene"], tran, new_start, end])
                break
            else:
                # 累加区间
                result.append([row["Gene"], tran, start, end])
                accumulated_length += interval_length

        # 转为 DataFrame
        result_df = pd.DataFrame(result, columns=["Gene", "Transcript", "CDS_start", "CDS_end"])
        # 查看结果
        new_res = pd.concat([new_res, result_df])
    new_res.to_csv(f"{project_dir}/cds_mark/{nc_no}_region.tsv",header=None,index=False,sep="\t")
    new_res = new_res.sort_values("CDS_start")
    # print(new_res)
    return new_res

    
def mark_cds(project_dir: str, nc_no: str, gdb_df: pd.DataFrame, cds_df: pd.DataFrame, output_file: str) -> None:
    # output_file = "exon_filter/exu.tsv"
    output_handle = open(output_file,'w')
    
    # 根据gdb 的基因和切点 找到cut的转录本、外显子以及cds并集前2/3位置
    # group_order = pd.Categorical(cds_df[1], categories=cds_df[1].unique(), ordered=True)
    # cds_df = cds_df.sort_values(by=[2, 3], key=lambda col: group_order if col.name == 2 else col)
    # cds_df[4] = cds_df.groupby(1).cumcount() + 1
    # cds_df = cds_df.sort_values(2)
    
    
    cds_df = common_cds_front_region(project_dir, cds_df, nc_no)
    cds_df.columns = [0,4,1,2]
    scaffold_pos_li = scaffold_detective_numpy(cds_df)
    # print(scaffold_pos_li)
    # cds_df.to_csv("cdsp.tsv",sep="\t",header=None,index=False)
    exon_array = cds_df.to_numpy()
    # 提取所需列并矢量化计算
    gdb_ori_array = gdb_df[4].astype(str) + "0.5"
    gdb_cut_pos_array = gdb_df[7].to_numpy()
    gdb_gene_type_array = gdb_df[10].to_numpy()
    gdb_array = gdb_df.to_numpy()
    
    gene_start_array = cds_df[1].to_numpy()
    gene_end_array = cds_df[2].to_numpy()

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
        # sub_cds_df = cds_df[
        #     (cds_df[0] == gene_name) & (g_pos > cds_df[2]) & (g_pos < cds_df[3])
        # ]
        # print(f"{idx}/{len(gdb_gene_name_array)} finished !!!",end="\r",flush=True)
        # if sub_cds_df.empty:
        #     res_count += 1
        gene_start = gene_start_array[current_gene_pos_idx]
        gene_end = gene_end_array[current_gene_pos_idx]
        if  gdb_gene_type_array[g_idx] == "non_coding" or g_pos < gene_start or g_pos > gene_end_array[-1]:
            print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
            print("no",end="\n",file=output_handle)
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
        # if current_gene_pos_idx == gene_pos_len -1 and g_pos > gene_end_array[-1]:
        #     break
        
        # 处于scaffold中的gRNA处理
        if current_scaffold_pos_idx < scaffold_pos_len and scaffold_pos_li[current_scaffold_pos_idx][0] <= current_gene_pos_idx <= scaffold_pos_li[current_scaffold_pos_idx][1]:
            fall_flag = 0
            for j in range(current_gene_pos_idx, scaffold_pos_li[current_scaffold_pos_idx][1] + 1):
                # if gene_start_array[j] < g_pos < gene_end_array[j]:
                if gene_start_array[j] < g_pos < gene_end_array[j]:
                    if gdb_array[g_idx][8] == exon_array[j][0]:
                        fall_flag = 1
                        print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                        # print(exon_array[current_gene_pos_idx][2],end="\t",file=output_handle)
                        # print(exon_array[current_gene_pos_idx][3],end="\n",file=output_handle)
                        print("yes",end="\n",file=output_handle)
                        break
            
            if fall_flag == 0:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                print("no",end="\n",file=output_handle)
            continue
        
        
        # 如果 gRNA 位于当前基因范围内，记录信息 chr1:1335323 跨两个或以上的位点会被跳过
        if gene_start < g_pos < gene_end:
            
            # print(f"{g_pos} in {gene_start[current_gene_pos_idx]}-{gene_end[current_gene_pos_idx]}", flush=True, end="\r")
            # print(f"{g_pos} in {gene_start}-{gene_end}")
            # g_item = gdb.iloc[g_idx]
            # g_item = gdb_array[g_idx]
            if gdb_array[g_idx][8] == exon_array[current_gene_pos_idx][0]:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                # print(exon_array[current_gene_pos_idx][2],end="\t",file=output_handle)
                # print(exon_array[current_gene_pos_idx][3],end="\n",file=output_handle)
                print("yes",end="\n",file=output_handle)
            else:
                print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
                print("no",end="\n",file=output_handle)
            #     print(gdb_array[g_idx][8],end="\t")
            #     print(exon_array[current_gene_pos_idx][0],end="\t")
            #     print(g_pos)
            # ...
            # print(g_item)
        else:
            print("\t".join(str(x) for x in gdb_array[g_idx]),end="\t",file=output_handle)
            print("no",end="\n",file=output_handle)
            # print(f"{g_pos} not in {gene_start}-{gene_end} and {gene_start_array[current_gene_pos_idx - 1]}-{gene_end_array[current_gene_pos_idx - 1]}")
            # g_item = gdb_array[g_idx]
            # print(g_item)
            ...
        
        
    return 

    
    
def top_cds_mark(project_dir: str, nc_no:str) -> None:
    gdb_file = f"{project_dir}/intron_filtered/{nc_no}.tsv"
    gdb_df = gdb2df(gdb_file)
    # 读取对应的CDS exon 表
    cds_file = f"{project_dir}/GCF/gtf/{nc_no}/CDS.tsv"
    cds_df = region2df(cds_file)
    # cds_df.to_csv("./nc22.cds",header=None,sep="\t")
    output_file = f"{project_dir}/cds_mark/{nc_no}.tsv"
    mark_cds(project_dir, nc_no, gdb_df,cds_df,output_file)



def main() -> None:
    nc2chr_file = "/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/nc_li"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    top_cds_mark(nc_li[-1])
    
    

if __name__ == "__main__":
    main()