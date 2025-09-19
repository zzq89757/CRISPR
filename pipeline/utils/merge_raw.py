import pandas as pd
from pathlib import Path

def merge_with_order(project_dir: str, nc_li: list[str], out_file: str):
    """
    用 pandas 合并 {project_dir}/ref_scan/{sample}.tsv，
    保持 nc_li 的顺序，并加上全局行号作为索引列。
    """
    df_li = []
    for sample in nc_li:
        fpath = Path(project_dir) / "ref_scan" / f"{sample}.tsv"
        if not fpath.exists():
            print(f"[WARN] {fpath} 不存在，跳过")
            continue
        df = pd.read_csv(fpath, sep="\t", header=None, dtype=str)  # 保证纯文本
        df_li.append(df)

    if not df_li:
        print("[ERROR] 没有有效的文件被读取")
        return

    merged = pd.concat(df_li, ignore_index=True)
    merged.index = merged.index + 1  # 行号从 1 开始
    merged.reset_index(names="ID", inplace=True)

    merged.to_csv(out_file, sep="\t", index=False, header=False)



if __name__ == "__main__":
    project_dir = "your_project"
    nc_li = ["nc_001", "nc_002", "nc_003"]  # 按顺序
    out_file = Path(project_dir) / "ref_scan" / "all.tsv"
    merge_with_order(project_dir, nc_li, out_file)
