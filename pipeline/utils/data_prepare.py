from os import system
from pathlib import Path
from pysam import FastaFile
import pandas as pd
from collections import defaultdict


def gtf2df(gtf_file: str) -> pd.DataFrame:
    """提前预设每列的数据类型并将gtf文件存入DataFrame"""
    use_col_li = [0, 2, 3, 4, 6, 8]
    type_li = ["category", "category", "int32", "int32", "category", "string"]

    type_dict = dict(zip(use_col_li, type_li))

    gtf_df = pd.read_csv(
        gtf_file,
        sep="\t",
        header=None,
        usecols=use_col_li,
        dtype=type_dict,
        # low_memory=False,
        comment="#",
    )

    return gtf_df


def append_gene_id_col(gtf_df: pd.DataFrame) -> None:
    """为gtf添加gene id 和 tran id 列"""
    gtf_df[[10 ,9]] = gtf_df[8].str.extract(
        r'transcript_id "([^"]*)".*gene "([^"]*)";'
    )
    gtf_df[9] = gtf_df[9].fillna('')
    gtf_df[10] = gtf_df[10].fillna('')


def obtain_gene_id_li(gtf_df: pd.DataFrame) -> defaultdict:
    """挑出gene行,只保留gene biotype 为protein_coding和ncRNA的gene id 列表并统计数目"""
    gene_df = gtf_df[gtf_df[2] == "gene"]
    contains_bool_li1 = gene_df[8].str.contains(r"gene_biotype \"protein_coding\";")
    contains_bool_li2 = gene_df[8].str.contains(r"gene_biotype \"ncRNA\";")
    contains_bool_li3 = gene_df[8].str.contains(r"gene_biotype \"lncRNA\";")
    contains_bool_li = [x or y or z for x, y, z in zip(contains_bool_li1, contains_bool_li2, contains_bool_li3)]
    gene_df = gene_df[contains_bool_li]
    gene_id_li = gene_df[9].to_list()
    return gene_id_li


def process_gtf(project_dir: str, nc_no: str) -> None:
    """
    将gtf拆分为CDS、EXON、TRAN以及GENE
    """
    Path(f"{project_dir}/GCF/gtf/{nc_no}").mkdir(exist_ok=True,parents=True)
    input_gtf = f"{project_dir}/GCF/gtf/{nc_no}.gtf"
    out_cds = f"{project_dir}/GCF/gtf/{nc_no}/CDS.tsv"
    out_exon = f"{project_dir}/GCF/gtf/{nc_no}/EXON.tsv"
    out_tran = f"{project_dir}/GCF/gtf/{nc_no}/TRAN.tsv"
    out_gene = f"{project_dir}/GCF/gtf/{nc_no}/GENE.tsv"

    # 读取gtf文件 生成df
    gtf_df = gtf2df(input_gtf)
    # 为gtf_df 添加gene id 和tran id列
    append_gene_id_col(gtf_df)
    # 找到coding和uncoding的gene id 列表
    gene_id_li = obtain_gene_id_li(gtf_df)

    # 挑出N开头的转录本和基因行(非X开头行)
    sub_df = gtf_df[~gtf_df[10].str.startswith("X")]
    
    all_gene_li=[]
    tran_exon_li = pd.DataFrame([])
    tran_cds_li = pd.DataFrame([])
    gene_tran_li = []
    for gene_name, sub_gene_df in sub_df.groupby(9,sort=False):
        # 跳过coding和nc之外的基因
        if gene_name not in gene_id_li:
            continue
        # 跳过没有NX or NM转录本的基因
        if len(sub_gene_df) == 1:
            print(f"{gene_name} has no NX or NM transcript found !!!")
            continue
        # 记录基因信息
        gene_item = sub_gene_df[[9, 3, 4, 6]].iloc[0]
        gene_item[10] = sub_gene_df.iloc[0][8].split('GeneID:')[1].split('"')[0] # ID
        gene_type_raw = sub_gene_df.iloc[0][8].split('gene_biotype "')[1].split('"')[0]
        gene_item[11] = "protein_coding" if gene_type_raw.find("RNA") == -1 else "non_coding" # type
        all_gene_li.append(gene_item)
        # 遍历按照转录本分类的子表（基因行的转录本id为空且只有一行）
        for tran_id, sub_tran_df in sub_gene_df.groupby(10,sort=False):
            # print(sub_tran_df)
            # 跳过基因行(基因行的转录本 id 为空 分组后仅一行)
            if len(sub_tran_df) == 1:
                continue
            # 保存转录本信息
            tran_df = sub_tran_df[[9, 10, 3, 4]].iloc[0]
            gene_tran_li.append(tran_df)
            # 保存exon 和 cds 信息
            exon_df = sub_tran_df[sub_tran_df[2]=="exon"][[9, 10, 3, 4]]
            cds_df = sub_tran_df[sub_tran_df[2]=="CDS"][[9, 10, 3, 4]]
            tran_cds_li = pd.concat([tran_cds_li, cds_df])
            tran_exon_li = pd.concat([tran_exon_li, exon_df])
    # 生成只含NM NR 的 gene表 cds表 和 exon表
    pd.DataFrame(all_gene_li).to_csv(out_gene,header=None,index=False,sep="\t")
    pd.DataFrame(gene_tran_li).to_csv(out_tran,header=None,index=False,sep="\t")
    pd.DataFrame(tran_cds_li).to_csv(out_cds,header=None,index=False,sep="\t")
    pd.DataFrame(tran_exon_li).to_csv(out_exon,header=None,index=False,sep="\t")




def data_prepare(project_dir: str, nc_no: str, chr_name: str, input_dict: dict) -> None:
    """
    根据nc号 切割reference、gtf以及vcf
    """
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
    # 将gtf拆分为CDS、GENE_LIST、EXON、TRAN等文件
    process_gtf(project_dir,nc_no)
    
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
    

    
    