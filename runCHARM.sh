#!/bin/bash
cd ../
mkdir -p slurm_log
snakemake --cluster 'sbatch --output=slurm_log/slurm-%j.out --cpus-per-task={threads} -t 7-00:00 -J CHARM!' --jobs 188 --resources nodes=188 --rerun-incomplete  -s ./CHARM/CHARM.smk
