#!/bin/bash
#SBATCH --account=nrsaamr
#SBATCH --partition=ord
#SBATCH --time=0-10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --job-name=snakemake_driver
#SBATCH --output=/work/NRSAAMR/Projects/SWAM/WGS/slurm/snakemake.%j

cd /work/NRSAAMR/Projects/SWAM/WGS

snakemake --profile config/slurm/ --config in_dir="input/" out_dir="output/" --configfile config/config.yaml --use-conda --conda-frontend conda \
	  --group-components group1=100 group2=100 group3=100 group4=100 -j 10

#snakemake --report report.html --config in_dir="input/" out_dir="output/" --configfile config/config.yaml
