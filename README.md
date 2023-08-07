# CHARM
CHARM project data pre-process pipeline in XingLab

## What is CHARM

Single cell **C**hromatin structure\\**H**istone modification\\**A**ccessibility and **R**NA expression **M**ulti-omics measurement(**CHARM**)

## Build enviroment for CHARM

```
# run this for the first time 
mamba create -n charm -c conda-forge -c bioconda python=3.7.16 snakemake=5.20.1 
mamba env update -n charm --file charm.yaml 
```
