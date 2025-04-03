from collections import defaultdict



# 根据exon（含cds和UTR）和cds 计算5'UTR和3'UTR位置
def utr_region_obtain(exon_file_path: str, cds_file_path: str) -> defaultdict:
    # 读取exon和cds region 文件
    ...
    # 判断基因方向并寻找UTR区域



# 根据UTR位置 添加新列 region 标记UTR和CDS
def utr_mark(nc_no: str) -> None:
    ...
