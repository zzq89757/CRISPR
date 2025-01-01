import os
import time
import pandas as pd
import psutil
from generate_split_ori import async_in_iterable_structure


def gene_in_li(gene_info: str, gene_li: list) -> bool:
    if str(gene_info).find("|") != -1:
        gene_info_li = [x.split(":")[0] for x in gene_info.split("|")]
        for gene in gene_info_li:
            if gene in gene_li:
                return True
    else:
        if str(gene_info) in gene_li:
            return True
    return False


# gene_in_li("LOC107987347:107987347|LOC107987346:107987346",[])

# exit()
# 读取 VCF 文件
def process_vcf(vcf_file: str, gene_file: str, output_tsv: str) -> None:
    # 读取基因位置信息
    gene_pos = pd.read_csv(gene_file,sep="\t",header=None,usecols=[0])
    gene_li = gene_pos[0].to_numpy()
    # 指定列名（VCF 文件无表头）
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    
    # 使用 pandas 读取 VCF 文件
    df = pd.read_csv(vcf_file, sep='\t', comment='#', names=columns, dtype=str)
    # 提取需要的列
    extracted = df[["CHROM", "POS", "REF", "ALT", "INFO"]]

    # 提取 INFO 列中 VC= 和 GENEINFO= 的值
    extracted[['GENEINFO', 'VC']] = extracted['INFO'].str.extract(r'GENEINFO=([^;]+):.*VC=([^;]+)')
    
    extracted = extracted[["CHROM", "POS", "REF", "ALT", "GENEINFO"]]

    # 寻找 GENEINFO 为 NaN 的索引
    nan_indices = extracted[extracted['GENEINFO'].isna()].index
    # 初始化新 DataFrame
    nan_indices_group = []
    # 确保 POS 列为数值型
    extracted['POS'] = pd.to_numeric(extracted['POS'])
    # 遍历基因区域
    current_group_idx = 0
    for idx, nan_idx in enumerate(nan_indices):
        if idx == 0:
            nan_indices_group.append([nan_idx])
            continue
        if nan_idx == nan_indices[idx - 1] + 1:
            nan_indices_group[current_group_idx].append(nan_idx)
        else:
            current_group_idx += 1
            nan_indices_group.append([nan_idx])

    # 以首尾为基准 挑选基准范围+-16bp的值
    filter_na_pos_li = []
    for group in nan_indices_group:
        filter_na_pos_li += [group[0],group[-1]]
        left_base_pos = extracted['POS'][group[0]]
        right_base_pos = extracted['POS'][group[-1]]
        
        for i in group[1:-1]:
            if extracted['POS'][i] <= left_base_pos + 16 or extracted['POS'][i] >= right_base_pos - 16:
                filter_na_pos_li.append(i)
        
    

    # 过滤掉不含 GENEINFO 和不在上下游16bp的行 (需要将以|分隔的基因名进行额外判断)
    extracted = extracted[extracted['GENEINFO'].apply(lambda x:gene_in_li(x,gene_li)) | extracted.index.isin(filter_na_pos_li)]

    # print(len(filter_na_pos_li))
    # 去掉 多余 列
    extracted = extracted[["POS", "REF", "ALT", "GENEINFO"]]
    extracted = extracted.sort_values("POS")

    # 保存为 TSV 文件
    extracted.to_csv(output_tsv, sep='\t', index=False)
    
    
def run_filter(nc_no: str) -> None:
    t1 = time.time()
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    # 使用示例
    vcf_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/vcf_split/raw/NC_000024.10.tsv"  # 输入的 VCF 文件路径
    vcf_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/vcf_split/raw/{nc_no}.tsv"  # 输入的 VCF 文件路径
    gene_pos_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/NC_000024.10/Gene_list.tsv"  # 输入的 VCF 文件路径
    gene_pos_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/Gene_list.tsv"  # 输入的 VCF 文件路径
    output_tsv = "/mnt/ntc_data/wayne/Repositories/CRISPR/vcf_split/filter/NC_000024.10.tsv"   # 输出的 TSV 文件路径
    output_tsv = f"/mnt/ntc_data/wayne/Repositories/CRISPR/vcf_split/filter/{nc_no}.tsv"   # 输出的 TSV 文件路径
    process_vcf(vcf_file, gene_pos_file, output_tsv)
    peak_memory_gb = memory_info.peak_wset / (1024**3) if hasattr(memory_info, 'peak_wset') else memory_info.rss / (1024**3)

    print(f"process {nc_no} peak memory cost: {peak_memory_gb:.2f} GB")
    print(f"process {nc_no} time cost:{time.time() - t1}")

def main() -> None:
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # run_filter(nc_li[-1])
    async_in_iterable_structure(run_filter,nc_li,24)


if __name__ == "__main__":
    main()
