# cas-offinder 试运行 单条约120s
# /mnt/ntc_data/wayne/Software/cas/cas-offinder input.txt C output.txt
# ./cas-offinder input.txt C output.txt
# cas-offinder 优化策略 1.删除参考中多余的scaffold和不含gRNA的区域

# 不使用cas-offinder 直接扫原始库 速度慢 读取数据库耗费时间远超120s

# 先统计原始库各个gRNA的数目 随后聚类分析？

# 切kmer 若最多三错配 20bp可切成5mer 即 16 个 数字的numpy数组 此方法不如直接将序列转为numpy数组

# 在id中记录gRNA序列的多个位置信息 构建blast数据库或编辑参考后使用bwa

# flashfry 试运行 chr22 index cost time ~10min individual gRNA search and score cost time ~3s

    # 构建索引 
#    java -Xmx4g -jar FlashFry-assembly-1.15.jar  index  --tmpLocation ./tmp  --database chr22_cas9ngg_database  --reference chr22.fa.gz
    # 搜索脱靶位点序列
#    java -Xmx4g -jar FlashFry-assembly-1.15.jar  discover  --database chr22_cas9ngg_database  --fasta EMX1_GAGTCCGAGCAGAAGAAGAAGGG.fasta  --output EMX1_ot.tsv
    # 打分
#    java -Xmx4g -jar FlashFry-assembly-1.15.jar  score  --input EMX1_ot.tsv  --output EMX1_ot_scored.tsv  --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot  --database chr22_cas9ngg_database 
#    java -Xmx4g -jar FlashFry-assembly-1.15.jar  score  --input EMX1_ot.tsv  --output EMX1_ot_scored.tsv  --scoringMetrics doench2016cfd  --database chr22_cas9ngg_database 


# y test <chr1 建库 120min y扫chr1 10min>
 # 构建y 对应的索引
#  ln -s /mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/split_fa/NC_000024.10.fa ./
 java -Xmx8g -jar FlashFry-assembly-1.15.jar  index  --tmpLocation ./tmp  --database /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/NC_000024.10_cas9_db --reference NC_000024.10.fa --enzyme spcas9
 # 提取y中的序列 生成fasta文件
 python dbseq2fasta.py
 # 搜索脱靶位点序列
 java -Xmx8g -jar FlashFry-assembly-1.15.jar  discover  --database /mnt/ntc_data/wayne/Repositories/CRISPR/sites_found/flashfry/db/NC_000024.10_cas9_db  --fasta y_sg.fasta  --output y_sg.tsv --maxMismatch 3 --maximumOffTargets 100 --forceLinear
 java -Xmx8g -jar FlashFry-assembly-1.15.jar  discover  --database NCA_cas9_db  --fasta y_sg.fasta  --output ya_sg.tsv --maxMismatch 3 --maximumOffTargets 100 --forceLinear
 java -Xmx8g -jar FlashFry-assembly-1.15.jar  discover  --database NCA_cas9_db  --fasta NC_000001.11.fa  --output 1a_sg.tsv --maxMismatch 3 --maximumOffTargets 100 --forceLinear
 