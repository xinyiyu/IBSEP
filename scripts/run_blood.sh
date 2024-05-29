method=tram
shrink=0.9
hsq=0.8
for chrom in {22..22}
do
    python /home/xiaojiashun/eQTL/IBSEP/scripts/real_data_debug.py \
        --data-dir /data2/xiaojiashun/eQTL/data/blood_data/1M/merged_sumstats/by_chrom_gene/chr${chrom} \
        --out-dir /data2/xiaojiashun/eQTL/data/blood_data/1M/merged_sumstats/omega_shrink_${shrink}_max_hsq_${hsq}/chr${chrom} \
        --cell-type-props /data2/xiaojiashun/eQTL/data/xcell/xcell_6meanprops_20231105.csv \
        --omega-shrink ${shrink} --max-hsq ${hsq} --method ${method} --trun-corr --pval-thres 1e-10 --reg-int-ident --tau 1.0
done
