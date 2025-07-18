## gRNA数据库(spCas9.a/i)构建流程

### 0.数据预处理
总览：需要寻找每个基因对应的转录本，并将每个转录本的起始位置作为TSS，而后将基因对应的每个(TSS - 1000,TSS)区域取并集，筛出落在并集区的gRNA

#### 0.1 根据gtf寻找TSS及其对应转录本
脚本路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/obtain_tss.py](./obtain_tss.py)

结果路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/tss_tran](./tss_tran)

#### 0.2 寻找TSS上游1kb 并获取交集区域 记录落在交集区能影响的转录本
脚本路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/insertion_region_obtain.py](./insertion_region_obtain.py)

结果路径：[/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/tss_regions](./tss_regions)

