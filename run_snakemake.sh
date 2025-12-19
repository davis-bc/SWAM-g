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

snakemake --profile config/slurm/ \
          --config in_dir="$input" out_dir="$output" \
          --use-conda \
          --conda-frontend conda \
          -j 10 \
          --local-cores 1 \
          --group-components group1=3 group2=3 group3=3 \
	  --keep-incomplete --rerun-incomplete
#	  --cleanup-metadata /work/NRSAAMR/Projects/SWAM/WGS/output/data/gtdb-tk/gtdbtk.bac120.summary.tsv /work/NRSAAMR/Projects/SWAM/WGS/output/data/checkm/genome.stats.tsv /work/NRSAAMR/Projects/SWAM/WGS/output/data/checkm/storage/bin_stats.analyze.tsv

#snakemake --report report.html --config in_dir="$input" out_dir="$output" --configfile config/config.yaml
