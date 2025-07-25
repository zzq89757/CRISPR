## gRNA数据库(spCas9.a/i)构建流程

### 0.数据预处理
总览：需要寻找每个基因对应的转录本，并将每个转录本的起始位置作为TSS，而后将基因对应的每个(TSS - 1000,TSS)区域取并集，筛出落在并集区的gRNA

#### 0.1 根据gtf寻找TSS及其对应转录本
脚本路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/obtain_tss.py](./obtain_tss.py)

结果路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/tss_tran](./tss_tran)

#### 0.2 寻找TSS上游1kb 并获取交集区域
脚本路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/insertion_region_obtain.py](./insertion_region_obtain.py)

结果路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/tss_regions](./tss_regions)


#### 0.3 去除PAM为NAG的gRNA
脚本路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/remove_nag_async.py](./remove_nag_async.py)

结果路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/nag_remove](./nag_remove)

### 1.CFD Scoring
计算每条gRNA的CFD分数

#### 1.0 off target sites search
脚本路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/search_ot.py](./search_ot.py)

结果路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/flashfry_out/](./cfd_score/flashfry_out/)

#### 1.1 过滤OVERFLOW和0 mismatch的结果并算分
脚本路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/filter_ot.py](./filter_ot.py)

结果路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/cfd_score/flashfry_filter_out/](./cfd_score/flashfry_filter_out)


### 2.rs2 Scoring