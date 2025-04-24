from pathlib import Path
import time
import pandas as pd
from collections import defaultdict
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from cds_mark import gene_ori_dict
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure



