for chrom in {1..22}
do
    python ibsep/prioritize_batch.py \
        --data-dir /data2/xiaojiashun/eQTL/data/blood_data/1M_neutro2015/merged_sumstats/by_chrom_gene_input/chr${chrom} \
        --out-dir /data2/xiaojiashun/eQTL/data/blood_data/1M_neutro2015/merged_sumstats/by_chrom_gene_output/chr${chrom} \
        --avg-props /data2/xiaojiashun/eQTL/data/cibersortx/blood_meanprops7.csv \
        --trun-corr
done
