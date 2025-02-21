from pathlib import Path
import pandas as pd
from sys import path
from os import system
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


def dual() -> pd.DataFrame:
    ...


def main() -> None:
    dual()


if __name__ == "__main__":
    main()