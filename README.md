# IBSEP: Integrating bulk RNA-seq and scRNA-seq data for cell-type-specific eQTLs prioritization
This repository holds the code for the paper "A unified framework for cell-type-specific eQTLs prioritization by integrating bulk and scRNA-seq data".

## Description
Traditional eQTL analyses using bulk RNA sequencing (bulk RNA-seq) provide tissue-level insights but suffer from signal loss and distortion due to unaddressed cellular heterogeneity. Recently, single-cell RNA sequencing (scRNA-seq) has provided higher resolution enabling cell-type-specific eQTL (ct-eQTL) analyses. However, these studies are limited by their smaller sample sizes and technical constraints. We present a novel statistical framework, IBSEP, which integrates bulk RNA-seq and scRNA-seq data for enhanced ct-eQTLs prioritization. IBSEP demonstrates superior performance in identifying ct-eQTLs and unveils new transcriptional regulatory mechanisms specific to cell types, offering deeper insights into the genetic basis of complex diseases at a cellular resolution.

## Installation

Git clone the repository and install the package:

```
conda env create -f environment.yml
conda activate ibsep
python setup.py develop
```

Check if it is successfully installed:

```
python ibsep/prioritize.py -h
```

## Quick start

### Data preparation

IBSEP performs ct-eQTLs prioritization one gene at a time. For each gene, IBSEP requires the following input files:

* eQTLs summary statistics for cell type 1, ..., $K$
* Tissue-level eQTLs summary statistics
* Average proportions of cell type 1, ..., $K$ in tissue samples

Suppose we apply IBSEP on blood ct-eQTLs prioritization using public eQTLs summary statistics for B cells, CD4 cells, CD8 cells, natural killer (NK) cells, monocytes and dendritic cells (DC), and public blood tissue eQTLs summary statistics (e.g from GTEx). For one gene (e.g ENSG00000104894), its eQTLs file for B cell looks like:

``` shell
$ head ENSG00000104894_B.csv
CHR     BP      A1      A2      SNP     GENEID  GENENAME        BETA    SE      PVAL
19      49742632        T       C       rs4802592       ENSG00000104894 CD37    0.0130254061242457      0.098058067569092       0.8941934
19      49743782        A       C       rs2878364       ENSG00000104894 CD37    0.0360282990701006      0.098058067569092       0.7133207
19      49744538        A       C       rs2335186       ENSG00000104894 CD37    0.0346677139651592      0.098058067569092       0.7237886
19      49744995        A       T       rs1716276       ENSG00000104894 CD37    -0.0020510021528886     0.098058067569092       0.9834052
19      49745289        T       G       rs4318314       ENSG00000104894 CD37    0.033686819503652       0.098058067569092       0.7312976
```

The eQTLs files for other cell types are of the same form and are named as `<gene_id>_<cell_type>.csv`. Columns `CHR, BP, A1, A2, SNP, BETA, SE` are required. The blood tisue eQTL summary statistics also has a similar format and is named as `<gene_id>_gtex_withld.csv` if GTEx data is used:

``` shell
$ head ENSG00000104894_gtex_withld.csv
CHR     BP      A1      A2      SNP     GENEID  GENENAME        BETA    SE      PVAL    LD
19      49742632        T       C       rs4802592       ENSG00000104894 CD37    0.0225815       0.0183049       0.217823        50.2383728027
19      49743782        A       C       rs2878364       ENSG00000104894 CD37    -0.0114515      0.0179085       0.522777        49.6510276794
19      49744538        A       C       rs2335186       ENSG00000104894 CD37    -0.00790063     0.0179737       0.66041 49.65102005
19      49744995        A       T       rs1716276       ENSG00000104894 CD37    0.0881959       0.0212669       3.85289e-05     54.7710914612
19      49745289        T       G       rs4318314       ENSG00000104894 CD37    -0.00775093     0.0180235       0.667316        49.65102005
```

Here, an additional column `LD` is present as IBSEP requires LD scores. LD scores can be calculated using the 1000 Genome reference panel (for European population). Before running IBSEP, we need to estimate proportions of cell types of interest in tissue samples using some cell type deconvolution method (e.g xCell/Cibersortx). Then we obtain the average cell type proportions:

``` shell
$ cat avg_props.csv
cell_type       prop
B       0.20722903601366718
CD4T    0.21487472191709173
CD8T    0.17785976564054523
DC      0.016122039530804287
monocyte        0.1519457097589802
NK      0.07374181989335338
```

### Usage
IBSEP integrates cell-type-level eQTLs summary statistics and tissue-level eQTLs summary statistics to prioritize ct-eQTLs. Specify `<data_dir>` containing summary statistics, `<out_dir>` for saving IBSEP output files and `gene_id` (e.g ENSG00000104894). As mentioned in our paper, IBSEP performs test on genetic correlations between cell types and set a small $p$-value threshold so that insignificant correlations are truncated to zero. We also recommand restricting the intercept matrix of LDSC to be an identity matrix if no strong evidence on sample overlap. The command for running IBSEP on one gene is as follow:

``` shell
python ibsep/prioritize.py \
    --data-dir <data_dir> \
    --out-dir <out_dir> \
    --avg-props avg_props.csv \
    --gene-id gene_id \
    --trun-corr
```

### Output
The major output of IBSEP contains improved ct-eQTLs summary statistics of cell types 1, ..., $K$. Here is an example output summary statistics file:

``` shell
$ head ENSG00000104894_B_IBSEP_pval1e-10_truncorr.csv
CHR     BP      A1      A2      SNP     GENEID  GENENAME        BETA    SE      PVAL     BETA_BLUE       SE_BLUE PVAL_BLUE
19      49742632        T       C       rs4802592       ENSG00000104894 CD37    0.0130254061242457      0.098058067569092       0.8941934       0.06599009919719703  0.0656300056917588      0.3146625406113458
19      49743782        A       C       rs2878364       ENSG00000104894 CD37    0.0360282990701006      0.098058067569092       0.7133207       -0.015352725024721278        0.06483390165999987     0.8128112681380122
19      49744538        A       C       rs2335186       ENSG00000104894 CD37    0.0346677139651592      0.098058067569092       0.7237886       -0.006172963972433082        0.06496643986209406     0.9243007400910752
19      49744995        A       T       rs1716276       ENSG00000104894 CD37    -0.0020510021528886     0.098058067569092       0.9834052       0.20204551346764305  0.07089711105082191     0.0043741051405087735
19      49745289        T       G       rs4318314       ENSG00000104894 CD37    0.033686819503652       0.098058067569092       0.7312976       -0.006101310669275541        0.0650672474336624      0.9252924140723515
```

### Example
The [example scripts](https://github.com/xinyiyu/IBSEP/tree/main/scripts) and [example data](https://github.com/xinyiyu/IBSEP/tree/main/examples) demonstrate how IBSEP works on one gene and genes in one chromosome.

## Reproducibility
We provide [scripts](https://github.com/xinyiyu/IBSEP/tree/main/scripts) and [notebooks](https://github.com/xinyiyu/IBSEP/tree/main/notebooks) for reproducing the real data analysis results of IBSEP. 

## Reference
Xinyi Yu, Xianghong HU, Xiaomeng Wan, Zhiyong Zhang, Mingxuan Cai, Tianwei Yu, Jiashun Xiao. A unified framework for cell-type-specific eQTLs prioritization by integrating bulk and scRNA-seq data. 2024. doi: https://doi.org/10.1101/2024.05.27.24307972.

## Contact information
Please contact Xinyi Yu (xyyu98@gmail.com) and Dr. Jiashun Xiao (xiaojiashun@sribd.cn) if any enquiry.
