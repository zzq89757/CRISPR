from collections import defaultdict
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from cds_mark import gene_ori_dict
from utils.read_tsv import tsv2df


# 根据exon（含cds和UTR）和cds 计算5'UTR和3'UTR位置
def utr_region_obtain(exon_file_path: str, cds_file_path: str, ori_dict: dict) -> defaultdict:
    # 读取exon和cds region 文件
    exon_df = tsv2df(exon_file_path,[])
    cds_df = tsv2df(cds_file_path,[])
    # 按照基因分组
    for gene, sub_exon_df in exon_df.groupby(0,sort=False):
        # 若无cds（NR） 跳过
        sub_cds_df = cds_df[cds_df[0]==gene]
        if sub_cds_df.empty:continue
        print(sub_exon_df)
        print(sub_cds_df)
        # 根据转录本分组 根据方向分别寻找UTR区域
        
        # 
        



# 根据UTR位置 添加新列 region 标记UTR和CDS
def utr_mark(nc_no: str) -> None:
    exon_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/EXON.tsv"
    cds_file = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/{nc_no}/CDS.tsv"
    ori_dict = gene_ori_dict(nc_no)
    utr_pos_dict = utr_region_obtain(exon_file, cds_file, ori_dict)
    ...
    

def main() -> None:
    nc_no = "NC_000024.10"
    
    utr_mark(nc_no)


if __name__ == "__main__":
    main()
