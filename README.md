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

Suppose we apply IBSEP on brain ct-eQTLs prioritization using public eQTLs summary statistics for excitatory neurons, inhibitory neurons, astrocytes, oligodendrocytes, microglia, OPCs, endothelial, pericytes, and GTEx brain cortex eQTL summary statistics. For one gene (e.g ENSG00000117500), its eQTLs file for excitatory neurons looks like:

``` shell
$ head ENSG00000117500_Excitatory_neurons.csv
CHR     BP      A1      A2      SNP     GENEID  GENENAME        BETA    SE      PVAL
1       92647223        A       G       rs560389        ENSG00000117500 TMED5   0.0708231       0.0465664584890307      0.128284
1       92647636        G       T       rs17578364      ENSG00000117500 TMED5   0.088615        0.0752085231859952      0.238694
1       92649292        A       T       rs570368        ENSG00000117500 TMED5   0.0708231       0.0465664584890307      0.128284
1       92650349        C       T       rs539288        ENSG00000117500 TMED5   0.0583558       0.0464225571977956      0.208733
1       92652691        T       C       rs112517106     ENSG00000117500 TMED5   -0.105728       0.0948710715648838      0.265091
```

The eQTLs files for other cell types are of the same form and are named as `<gene_id>_<cell_type>.csv`. Columns `CHR, BP, A1, A2, SNP, BETA, SE` are required. The tissue-level eQTL summary statistics also has a similar format and is named as `<gene_id>_gtex_withld.csv` if GTEx data is used:

``` shell
$ head ENSG00000117500_gtex_withld.csv
CHR     BP      A1      A2      SNP     GENEID  GENENAME        BETA    SE      PVAL    LD
1       92647223        A       G       rs560389        ENSG00000117500 TMED5   -0.0938936      0.0453816       0.0401076       140.784896851
1       92647636        G       T       rs17578364      ENSG00000117500 TMED5   -0.100632       0.0749145       0.181022        27.6296691895
1       92649292        A       T       rs570368        ENSG00000117500 TMED5   -0.0818084      0.0448585       0.070007        140.784896851
1       92650349        C       T       rs539288        ENSG00000117500 TMED5   -0.0729709      0.0454674       0.110427        157.327850342
1       92652691        T       C       rs112517106     ENSG00000117500 TMED5   -0.0994571      0.072868        0.174145        151.522323608
```

Here, an additional column `LD` is present as IBSEP requires LD scores. LD scores can be calculated using the 1000 Genome reference panel (for European population). Before running IBSEP, we need to estimate proportions of cell types of interest in tissue samples using some cell type deconvolution method (e.g xCell/Cibersortx). Then we obtain the average cell type proportions:

``` shell
$ cat brain_meanprops8.csv
cell_type       prop
Inhibitory neurons      0.102130737487027
Oligodendrocytes        0.169931676536728
Excitatory neurons      0.308488897051598
Astrocytes      0.290533778997567
OPCs    0.00379891008352746
Microglia       0.0490726269125998
Pericytes       0.024231706135807
Endothelial     0.0518116667951463
```

### Usage
IBSEP integrates cell-type-level eQTLs summary statistics and tissue-level eQTLs summary statistics to prioritize ct-eQTLs. Specify `<data_dir>` containing summary statistics, `<out_dir>` for saving IBSEP output files, a table `<avg_props.csv>` providing average cell type proportions and `<gene_id>` (e.g ENSG00000117500). As mentioned in our paper, IBSEP performs test on genetic correlations between cell types and set a small $p$-value threshold so that insignificant correlations are truncated to zero. We also recommand restricting the intercept matrix of LDSC to be an identity matrix if no strong evidence on sample overlap. The command for running IBSEP on one gene is as follow:

``` shell
python ibsep/prioritize.py \
    --data-dir <data_dir> \
    --out-dir <out_dir> \
    --avg-props <avg_props.csv> \
    --gene-id <gene_id> \
    --trun-corr
```

### Output
The major output of IBSEP contains improved ct-eQTLs summary statistics of cell types 1, ..., $K$. Here is an example output summary statistics file:

``` shell
$ head ENSG00000117500_Excitatory_neurons_IBSEP_truncorr_pval1e-10.csv
CHR     BP      A1      A2      SNP     GENEID  GENENAME        BETA    SE      PVAL    BETA_IBSEP      SE_IBSEP        PVAL_IBSEP
1       92647223        A       G       rs560389        ENSG00000117500 TMED5   0.0708231       0.0465664584890307      0.128284        0.015446450635914216    0.0260622271503747      0.5533981577249139
1       92647636        G       T       rs17578364      ENSG00000117500 TMED5   0.088615        0.0752085231859952      0.238694        -0.035877049952112265   0.03472936391585988     0.30158210979172717
1       92649292        A       T       rs570368        ENSG00000117500 TMED5   0.0708231       0.0465664584890307      0.128284        0.0175173325550378      0.026022131792201596    0.5008387999823593
1       92650349        C       T       rs539288        ENSG00000117500 TMED5   0.0583558       0.0464225571977956      0.208733        0.014585469956598996    0.0263835006060214      0.5803829562933673
1       92652691        T       C       rs112517106     ENSG00000117500 TMED5   -0.105728       0.0948710715648838      0.265091        -0.041378004994756364   0.043257948643890196    0.33879893928307037
```

The input, output and script for this example can be found at [example scripts](https://github.com/xinyiyu/IBSEP/tree/main/scripts) and [example data](https://github.com/xinyiyu/IBSEP/tree/main/examples).

## Reproducibility
We provide [scripts](https://github.com/xinyiyu/IBSEP/tree/main/scripts) and [notebooks](https://github.com/xinyiyu/IBSEP/tree/main/notebooks) for reproducing the real data analysis results of IBSEP. 

## Reference
Xinyi Yu, Xianghong HU, Xiaomeng Wan, Zhiyong Zhang, Mingxuan Cai, Tianwei Yu, Jiashun Xiao. A unified framework for cell-type-specific eQTLs prioritization by integrating bulk and scRNA-seq data, _The American Journal of Human Genetics_ (2025), https://doi.org/10.1016/j.ajhg.2024.12.018

## Contact information
Please contact Xinyi Yu (xyyu98@gmail.com) and Dr. Jiashun Xiao (xiaojiashun@sribd.cn) if any enquiry.
