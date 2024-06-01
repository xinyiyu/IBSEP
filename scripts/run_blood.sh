for chrom in {1..22}
do
    python ibsep/prioritize_batch.py \
        --data-dir /data2/xiaojiashun/eQTL/data/blood_data/1M/merged_sumstats/by_chrom_gene/chr${chrom} \
        --out-dir /data2/xiaojiashun/eQTL/data/blood_data/1M/merged_sumstats/by_chrom_gene/chr${chrom} \
        --avg-props /data2/xiaojiashun/eQTL/data/xcell/xcell_6meanprops.csv \
        --trun-corr
done
