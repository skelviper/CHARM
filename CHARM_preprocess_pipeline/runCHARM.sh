#!/bin/bash

#usage: ./runCHARM.sh

cd ../
mkdir -p slurm_log
snakemake --use-conda --cluster 'sbatch --qos=high -w node03 --output=slurm_log/slurm-%j.out --cpus-per-task={threads} -t 7-00:00 -J CHARM!' --jobs 1024 --resources nodes=1024 --rerun-incomplete -s ./CHARM_preprocess_pipeline/CHARM.smk --keep-going

mkdir -p ./analysis
cp CHARM_preprocess_pipeline/stat.ipynb ./analysis/stat.ipynb
