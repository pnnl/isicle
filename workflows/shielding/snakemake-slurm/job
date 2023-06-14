#!/bin/bash
#SBATCH -J isicle_nmr
#SBATCH -N 1
#SBATCH -t 96:00:00
#SBATCH -n 1
#SBATCH -p shared,slurm,short

source /etc/bashrc
module purge
module load intel/18.0.0
module load gcc/8.1.0
#module load java/1.8.0_31
# python
module load python/miniconda3.7
source /share/apps/python/miniconda3.7/etc/profile.d/conda.sh
conda activate isicle

snakemake --unlock
snakemake --cluster 'sbatch --job-name {resources.name} -t {resources.runtime} -p {resources.partition} -N {resources.nodes}' -j 5000 --latency-wait 100 --keep-going --rerun-incomplete

