#!/bin/bash
#SBATCH --account=nrsaamr
#SBATCH --partition=ord
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --job-name=snakemake_driver
#SBATCH --output=/work/NRSAAMR/Projects/SWAM/WGS/slurm/snakemake.%j

cd /work/NRSAAMR/Projects/SWAM/WGS

input="/work/NRSAAMR/Projects/SWAM/WGS/input"
output="/work/NRSAAMR/Projects/SWAM/WGS/output"

snakemake --profile config/slurm/ --config in_dir="$input" out_dir="$output" --use-conda --conda-frontend conda -j 10 --local-cores 1

#snakemake --report report.html --config in_dir="$input" out_dir="$output" --configfile config/config.yaml
