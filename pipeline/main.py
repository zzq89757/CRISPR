# encoding=utf-8
import click
import pandas as pd
from collections import defaultdict
from pathlib import Path
from sys import path
path.append("/mnt_data/Wayne/Repositories/CRISPR/")
from generate_split_ori import async_in_iterable_structure



# 添加命令行提示信息
@click.command()
@click.option('-l', '--nclist', required=False, type=str, help="nc list")
@click.option('-r', '--ref', required=False, type=str, help="raw Ref")
@click.option('-g', '--gtf', required=False, type=str, help="raw gtf")
@click.option('-g', '--vcf', required=False, type=str, help="raw gtf")
@click.option('-t', '--threads', default=18, type=int, help="threads of all programs")
@click.option('-o', '--outdir', required=False, type=str, help="directory of result")
class Base:
    def __init__(self,r1,r2,addref,addgtf,ref,gtf,threads,outdir) -> None:
        self.r1, self.r2, self.addref,self.addgtf, self.ref, self.gtf, self.threads ,self.outdir=  Path(r1).absolute(), Path(r2).absolute(), Path(addref).absolute(),Path(addgtf).absolute(), Path(ref).absolute(),Path(gtf).absolute(), threads, Path(outdir).absolute()
        # 按照核酸号拆分gtf、fna、vcf文件
        self.data_prepare()
    

    
    def data_prepare() -> None:
        """
        按照核酸号拆分gtf、fna、vcf文件
        """
        ...


if __name__=="__main__":
    Base()
