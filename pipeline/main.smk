from utils.data_prepare import data_prepare
from utils.ref_scan import ref_scan
from utils.merge_raw import merge_with_order
from utils.add_raw_id import add_raw_id
from utils.gene_annotate import gene_annotate
import pandas as pd
from pathlib import Path


# 读取配置文件信息
# configfile: "config.yaml"


project_dir = config.get("project_dir")
nc_li = pd.read_csv(config.get("nc_li"), header=None, sep="\t")[0].to_list()
chr_li = pd.read_csv(config.get("nc_li"), header=None, sep="\t")[1].to_list()
nc2chr_dict = dict(zip(nc_li, chr_li))
fna_file = config.get("fna_file")
gtf_file = config.get("gtf_file")
vcf_file = config.get("vcf_file")


rule all:
    input:
        # expand(
        #     "{project_dir}/GCF/gtf/{sample}.gtf", project_dir=project_dir, sample=nc_li
        # ),
        # expand(
        #     "{project_dir}/GCF/fa/{sample}.fa", project_dir=project_dir, sample=nc_li
        # ),
        # expand(
        #     "{project_dir}/GCF/vcf/{sample}.vcf", project_dir=project_dir, sample=nc_li
        # ),
        # expand(
        #     "{project_dir}/ref_scan/{sample}.tsv",
        #     project_dir=project_dir,
        #     sample=nc_li,
        # ),
        # f"{project_dir}/ref_scan/all.tsv"
        # f"{project_dir}/ref_scan/lc_all",
        # f"{project_dir}/raw_with_id/{sample}.tsv"
        # expand(
        #     "{project_dir}/raw_with_id/{sample}.tsv",
        #     project_dir=project_dir,
        #     sample=nc_li,
        # ),
        expand(
            "{project_dir}/gene_annotated/{sample}.tsv",
            project_dir=project_dir,
            sample=nc_li,
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
    params:
        chr_name=lambda wildcards: nc2chr_dict[wildcards.sample],
    run:
        Path(f"{project_dir}/GCF/fa").mkdir(exist_ok=True, parents=True)
        Path(f"{project_dir}/GCF/vcf").mkdir(exist_ok=True, parents=True)
        Path(f"{project_dir}/GCF/gtf").mkdir(exist_ok=True, parents=True)
        # 这里调用单个样本的处理逻辑
        data_prepare(project_dir, wildcards.sample, params.chr_name, input)


rule search_ref:
    input:
        ref="{project_dir}/GCF/fa/{sample}.fa",
    output:
        raw_db="{project_dir}/ref_scan/{sample}.tsv",
    params:
        chr_name=lambda wildcards: nc2chr_dict[wildcards.sample],
    run:
        Path(f"{project_dir}/ref_scan").mkdir(exist_ok=True, parents=True)
        ref_scan(wildcards.sample, params.chr_name, input.ref, output.raw_db)


rule generate_line_count_file:
    input:
        # raw_db="{project_dir}/ref_scan/{sample}.tsv",
        tsvs=expand(
            "{project_dir}/ref_scan/{sample}.tsv",
            project_dir=project_dir,
            sample=nc_li,
        ),
        # tsvs = [f"{project_dir}/ref_scan/{s}.tsv" for s in nc_li]
    output:
        tmp_lc_file="{project_dir}/ref_scan/lc_tmp",
        lc_file="{project_dir}/ref_scan/lc_all",
    shell:
        r"""
        wc -l {input.tsvs} > {output.tmp_lc_file}
        sed -i 's/^[[:space:]]*//' {output.tmp_lc_file}
        cut -d " " -f1 {output.tmp_lc_file} > {output.lc_file}
        """


rule add_raw_id:
    input:
        raw_db="{project_dir}/ref_scan/{sample}.tsv",
        lc_file="{project_dir}/ref_scan/lc_all",
    output:
        raw_db_with_id="{project_dir}/raw_with_id/{sample}.tsv",
    params:
        nc_idx=lambda wildcards: nc_li.index(wildcards.sample),
    run:
        Path(f"{project_dir}/raw_with_id").mkdir(exist_ok=True, parents=True)
        add_raw_id(params.nc_idx, input.lc_file, input.raw_db, output.raw_db_with_id)


rule gene_annotate:
    input:
        raw_db="{project_dir}/raw_with_id/{sample}.tsv",
    output:
        gene_annotated_db="{project_dir}/gene_annotated/{sample}.tsv",
    run:
        Path(f"{project_dir}/gene_annotated").mkdir(exist_ok=True, parents=True)
        gene_annotate(project_dir,wildcards.sample)
