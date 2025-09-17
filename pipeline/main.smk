from utils.data_prepare import data_prepare
from utils.ref_scan import ref_scan
import pandas as pd
from pathlib import Path


# 读取配置文件
configfile: "config.yaml"


project_dir = config.get("project_dir")
nc_li = (
    config.get("nc_li")
    if isinstance(config.get("nc_li"), list)
    else pd.read_csv(config.get("nc_li"), header=None)[0].to_list()
)
fna_file = config.get("fna_file")
gtf_file = config.get("gtf_file")
vcf_file = config.get("vcf_file")
print(nc_li)


rule all:
    input:
        expand(
            "{project_dir}/GCF/gtf/{sample}.gtf", project_dir=project_dir, sample=nc_li
        ),
        expand(
            "{project_dir}/GCF/fa/{sample}.fa", project_dir=project_dir, sample=nc_li
        ),
        expand(
            "{project_dir}/GCF/vcf/{sample}.vcf", project_dir=project_dir, sample=nc_li
        ),


rule data_prepare:
    input:
        fna_file=fna_file,
        gtf_file=gtf_file,
        vcf_file=vcf_file,
    output:
        gtf="{project_dir}/GCF/gtf/{sample}.gtf",
        fa="{project_dir}/GCF/fa/{sample}.fa",
        vcf="{project_dir}/GCF/vcf/{sample}.vcf",
    run:
        Path(f"{project_dir}/GCF/fa").mkdir(exist_ok=True,parents=True)
        Path(f"{project_dir}/GCF/vcf").mkdir(exist_ok=True,parents=True)
        Path(f"{project_dir}/GCF/gtf").mkdir(exist_ok=True,parents=True)
        # 这里调用单个样本的处理逻辑
        data_prepare(project_dir, wildcards.sample, input)


rule search_ref:
    input: 
        ref = "{project_dir}/GCF/fa/{sample}.fa"
    output: 
        raw_db = "{project_dir}/ref_scan/{sample}.tsv"
    run: 
        Path(f"{project_dir}/ref_scan").mkdir(exist_ok=True,parents=True)
        ref_scan(project_dir, input.ref, output.raw_db)