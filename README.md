# Code for CHARM data analysis

<img src="logo.png" width="50%" height="50%">

## Introduction 

A single cell assay for chromatin <ins>**C**</ins>onformation, <ins>**H**</ins>istone modification, chromatin <ins>**A**</ins>ccessibility, and <ins>**R**</ins>NA expression <ins>**M**</ins>ulti-omics profiling (CHARM).

## Related Publications

CHARM paper:

Gene regulatory landscape dissected by single-cell four-omics sequencing (In submission, hope to be published soon!)

Yujie Chen*, Zhiyuan Liu*, Heming Xu, Jiayu Liu, Mengxuan Wang, Yi Chi, Boyuan Liang, Menghan Liu, Yongli Peng, Hao Ge, Dong Xing✉

+ Raw data: PRJNA1284811
+ Processed data: GSE303006

Our previous multi-omics profiling method, HiRES:

Linking genome structures to functions by simultaneous single-cell Hi-C and RNA-seq, Science 380, 1070–1076 (2023)

Zhiyuan Liu*; Yujie Chen*; Qimin Xia*; Menghan Liu; Heming Xu; Yi Chi; Yujing Deng; Dong Xing✉

+ Raw data: PRJNA907173
+ Processed data: GSE223917

This folder project contains two folders, and I will describe their respective uses below.

## CHARM_preprocess_pipeline

The CHARM_preprocess_pipeline contains a snakemake pipeline for pre-processing CHARM data.

### Usage
1. Place the Rawdata folder and the CHARM_preprocess_pipeline folder in the same directory.
```bash
tree -h -L 2
# .
# ├── [ 252]  CHARM_preprocess_pipeline
# │   ├── [3.0K]  config.yaml
# │   ├── [  35]  envs
# │   ├── [ 248]  CHARM_scripts
# │   ├── [3.1K]  CHARM.smk
# │   ├── [ 152]  README.md
# │   ├── [ 224]  rules
# │   ├── [ 361]  runCHARM.sh
# │   └── [ 20K]  stat.ipynb
# |   ...
# └── [  33]  Rawdata
#     └── [  40]  P10001 -> ../../CHARM_test/Rawdata/P10001
```
2. Installation of environment

```bash
mamba create -n charm -c conda-forge -c bioconda python=3.8 snakemake=5.20.1 
mamba activate charm
mamba env update --file CHARM_preprocess_pipeline/envs/charm.yaml
```
3. Prepare files

CHARM_preprocess_pipeline relies on a modified version of hickit (original link: https://github.com/lh3/hickit/ , but no need to download it) and CHARMtools(https://github.com/skelviper/CHARMtools) in addition to softwares that can be installed automatically. Also you need to build index for RNA/DNA seperately on your version of reference genome.

```bash
cd CHARM_preprocess_pipeline
vim config.yaml
```

4. Run the pipeline
```bash
cd CHARM_preprocess_pipeline; ./runCHARM.sh
```

5. generate statistics

    see analysis/stat.ipynb

## analysis_and_plot_notebooks

analysis_and_plot_notebooks holds the python/R notebooks used in the analysis of the project.