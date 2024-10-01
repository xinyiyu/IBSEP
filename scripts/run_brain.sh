for chrom in {1..22}
do
    python ibsep/prioritize_batch.py \
        --data-dir /data2/xiaojiashun/eQTL/data/brain_data/merged_sumstats/by_chrom_gene_input/chr${chrom} \
        --out-dir /data2/xiaojiashun/eQTL/data/brain_data/merged_sumstats/by_chrom_gene_output/chr${chrom} \
        --cell-type-props /data2/xiaojiashun/eQTL/data/cibersortx/brain_meanprops8.csv \
        --trun-corr
done
