#!/bin/bash
#SBATCH --account=nrsaamr
#SBATCH --partition=ord
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --job-name=snakemake_driver
#SBATCH --output=/work/NRSAAMR/Projects/SWAM/WGS/snakemake.%j

snakemake --profile config/slurm/ --config in_dir="input/" out_dir="output/" --configfile config/config.yaml --use-conda --conda-frontend conda
