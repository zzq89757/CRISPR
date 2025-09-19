# zcat /mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/GCF/GCF_000001635.27_GRCm39_genomic.fna.gz | grep ">NC_0000" | cut -d" " -f1 | sed  's/>//' > mus/nc_li
# zcat /mnt_data/Wayne/Repositories/CRISPR/pipeline/rat/GCF/GCF_036323735.1_GRCr8_genomic.fna.gz | grep ">NC_0860" | cut -d" " -f1 | sed  's/>//' > rat/nc_li
snakemake -s main.smk -c 48 --configfile /mnt_data/Wayne/Repositories/CRISPR/pipeline/config_mus.yaml
snakemake -s main.smk -c 48 --configfile /mnt_data/Wayne/Repositories/CRISPR/pipeline/config_rat.yaml