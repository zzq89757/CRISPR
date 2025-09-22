from os import system
from pathlib import Path
from pysam import FastaFile

def process_gtf(project_dir: str, nc_no: str) -> None:
    ...


def data_prepare(project_dir: str, nc_no: str, chr_name: str, input_dict: dict) -> None:
    # 创建路径
    # deal_dir(project_dir)
    # 接受文件路径字典
    fna_file = input_dict.fna_file
    gtf_file = input_dict.gtf_file
    vcf_file = input_dict.vcf_file

    # 处理逻辑：
    # 考虑输入文件是否解压
    if gtf_file.endswith("gz"):
        system(f"zcat {gtf_file} | awk '$1==\"{nc_no}\"' > {project_dir}/GCF/gtf/{nc_no}.gtf")
    else:
        system(f"awk '$1==\"{nc_no}\"' {gtf_file} > {project_dir}/GCF/gtf/{nc_no}.gtf")
    # 根据nc号获取染色体
    chr_no = chr_name.replace("chr","")
    if vcf_file.endswith("gz"):
        system(f"zcat {vcf_file} | awk '$1==\"{chr_no}\"' > {project_dir}/GCF/vcf/{nc_no}.vcf")
    else:
        system(f"awk '$1==\"{chr_no}\"' {vcf_file} > {project_dir}/GCF/vcf/{nc_no}.vcf")
    
    # 对于fa 直接解压
    if fna_file.endswith("gz"):
        system(f"gzip -d {fna_file}")
    fna_file = fna_file.replace(".gz","")

    # 切割fa文件
    nc_seq = FastaFile(fna_file).fetch(nc_no).upper()

    output_fa = f"{project_dir}/GCF/fa/{nc_no}.fa"

    with open(output_fa,'w') as ofa:
        ofa.write(f">{nc_no}\n{nc_seq}\n")

    ofa.close()    
    
    # 将gtf拆分为CDS、GENE_LIST、EXON、TRAN等文件

    
    