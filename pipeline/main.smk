from utils.data_prepare import data_prepare
from utils.ref_scan import ref_scan
from utils.merge_raw import merge_with_order
from utils.add_raw_id import add_raw_id
from utils.gene_annotate import gene_annotate
from utils.filter_intron import exon_annotate
from utils.cds_mark import top_cds_mark
from utils.tran_count import tran_ratio_count
from utils.ag_end import ag_mark
from utils.rm_n0 import remove_n0
from utils.flashfry_seq_construct import construct_seq
from utils.filter_ot import filter_ot
from utils.cfd_score_append import append_cfd_score
from utils.flank_fill import flank_fill
import pandas as pd
from pathlib import Path
from os import system


# 读取配置文件信息
# configfile: "config.yaml"


project_dir = config.get("project_dir")
flashfry_bin = config.get("flashfry_bin")
nc_li = pd.read_csv(config.get("nc_li"), header=None, sep="\t")[0].to_list()
chr_li = pd.read_csv(config.get("nc_li"), header=None, sep="\t")[1].to_list()
nc2chr_dict = dict(zip(nc_li, chr_li))
chr2nc_dict = dict(zip(chr_li, nc_li))
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
        # f"{project_dir}/GCF/fa/NCA.fasta",
        expand(
            "{project_dir}/flank_fill/{sample}.tsv",
            project_dir=project_dir,
            sample=nc_li,
        ),
        f"{project_dir}/GCF/flashfry_db/NCA_cas9_db",


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

rule merge_fa:
    input: 
        # fa="{project_dir}/GCF/fa/{sample}.fa",
        expand(
            "{project_dir}/GCF/fa/{sample}.fa", project_dir=project_dir, sample=nc_li
        ),
    output: 
        fasta="{project_dir}/GCF/fa/NCA.fasta",
    run: 
        system(f"cat {project_dir}/GCF/fa/NC* > {project_dir}/GCF/fa/NCA.fasta")


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
        gene_annotate(project_dir, wildcards.sample)


rule filter_intron:
    input:
        gene_annotated_db="{project_dir}/gene_annotated/{sample}.tsv",
    output:
        intron_filtered_db="{project_dir}/intron_filtered/{sample}.tsv",
    run:
        Path(f"{project_dir}/intron_filtered").mkdir(exist_ok=True, parents=True)
        exon_annotate(project_dir, wildcards.sample)


rule top_cds_mark:
    input:
        intron_filtered_db="{project_dir}/intron_filtered/{sample}.tsv",
    output:
        cds_marked_db="{project_dir}/cds_mark/{sample}.tsv",
    run:
        Path(f"{project_dir}/cds_mark").mkdir(exist_ok=True, parents=True)
        top_cds_mark(project_dir, wildcards.sample)


rule tran_count:
    input:
        cds_marked_db="{project_dir}/cds_mark/{sample}.tsv",
    output:
        tran_counted_db="{project_dir}/tran_count/{sample}.tsv",
    run:
        Path(f"{project_dir}/tran_count").mkdir(exist_ok=True, parents=True)
        tran_ratio_count(project_dir, wildcards.sample)


rule ag_mark:
    input:
        tran_counted_db="{project_dir}/tran_count/{sample}.tsv",
    output:
        ag_marked_db="{project_dir}/ag_mark/{sample}.tsv",
    run:
        Path(f"{project_dir}/ag_mark").mkdir(exist_ok=True, parents=True)
        ag_mark(project_dir, wildcards.sample)


rule flashfry_index_build:
    input:
        nca_fa_file="{project_dir}/GCF/fa/NCA.fasta",
        # ag_marked_db="{project_dir}/ag_mark/{sample}.tsv",
        ag_marked_db = expand(
            "{project_dir}/ag_mark/{sample}.tsv", project_dir=project_dir, sample=nc_li
        ),
    output:
        flashfry_index="{project_dir}/GCF/flashfry_db/NCA_cas9_db",
    shell:
        # 若已经构建该物种flashfry index 直接使用 否则进行构建
        r"""
        mkdir -p {project_dir}/GCF/flashfry_db/tmp
        java -Xmx200g -jar {flashfry_bin} index --tmpLocation {project_dir}/GCF/flashfry_db/tmp --database {project_dir}/GCF/flashfry_db/NCA_cas9_db --reference {input.nca_fa_file} --enzyme spcas9
        """


rule remove_n0_seq:
    input: 
        flashfry_index=f"{project_dir}/GCF/flashfry_db/NCA_cas9_db",
    output: 
        ag_marked_rm0_db = expand(
            "{project_dir}/ag_mark_n0_rm/{sample}.tsv", project_dir=project_dir, sample=nc_li
        ),
    run: 
        Path(f"{project_dir}/ag_mark_n0_rm").mkdir(exist_ok=True, parents=True)
        remove_n0(project_dir,chr2nc_dict)


rule flashfry_input_seq_construct:
    input: 
        ag_marked_rm0_db = expand(
            "{project_dir}/ag_mark_n0_rm/{sample}.tsv", project_dir=project_dir, sample=nc_li
        ),
    output: 
        flashfry_seq_input="{project_dir}/flashfry_out/fa/{sample}.fa"
    run: 
        Path(f"{project_dir}/flashfry_out/fa").mkdir(exist_ok=True, parents=True)
        construct_seq(project_dir, wildcards.sample)

rule flashfry_search_ot:
    input: 
        flashfry_seq_input="{project_dir}/flashfry_out/fa/{sample}.fa",
    output: 
        flashfry_raw_out="{project_dir}/flashfry_out/raw_out/{sample}.tsv"
    shell:
        r"""
        mkdir -p {project_dir}/flashfry_out/raw_out
        java -Xmx200g -jar {flashfry_bin} discover --database {project_dir}/GCF/flashfry_db/NCA_cas9_db --fasta {input.flashfry_seq_input} --output {output.flashfry_raw_out} --maxMismatch 3 --maximumOffTargets 100 --forceLinear
        """

rule flashfry_offtarget_filter:
    input: 
        flashfry_raw_out="{project_dir}/flashfry_out/raw_out/{sample}.tsv"
    output: 
        flashfry_filter_out="{project_dir}/flashfry_out/filter_out/{sample}.tsv"
    run: 
        Path(f"{project_dir}/flashfry_out/filter_out").mkdir(exist_ok=True, parents=True)
        filter_ot(project_dir, wildcards.sample)


rule flashfry_score:
    input: 
        flashfry_filter_out="{project_dir}/flashfry_out/filter_out/{sample}.tsv"
    output: 
        flashfry_score_out="{project_dir}/flashfry_out/score_out/{sample}.tsv"
    shell:
        r"""
        mkdir -p {project_dir}/flashfry_out/score_out
        java -Xmx200g -jar {flashfry_bin} score --input {input.flashfry_filter_out} --output {output.flashfry_score_out} --scoringMetrics doench2016cfd --database {project_dir}/GCF/flashfry_db/NCA_cas9_db
        """

rule cfd_score:
    input: 
        ag_marked_db="{project_dir}/ag_mark_n0_rm/{sample}.tsv",
        flashfry_score_out="{project_dir}/flashfry_out/score_out/{sample}.tsv"
    output: 
        cfd_score_db="{project_dir}/cfd_score/{sample}.tsv",
    run: 
        Path(f"{project_dir}/cfd_score").mkdir(exist_ok=True, parents=True)
        append_cfd_score(project_dir, wildcards.sample)

rule flank_fill:
    input: 
        cfd_score_db="{project_dir}/cfd_score/{sample}.tsv"
    output: 
        flank_fill_db="{project_dir}/flank_fill/{sample}.tsv"
    run: 
        Path(f"{project_dir}/flank_fill").mkdir(exist_ok=True, parents=True)
        flank_fill(project_dir, wildcards.sample)

# rule rs2_score:
    
#     input: 
#     output: 
#     run: 


# rule snp_mark:
#     input: 
#     output: 
#     run: 


# rule utr_mark:
#     input: 
#     output: 
#     run: 


# rule low_score_mark:
#     input: 
#     output: 
#     run: 


# rule final_filter:
#     input: 
#     output: 
#     run: 