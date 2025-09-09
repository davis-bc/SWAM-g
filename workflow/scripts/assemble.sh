#!/bin/bash
#SBATCH --account=nrsaamr
#SBATCH --nodes=1
#SBATCH --partition=ord
#SBATCH --time=1-00:00:00 #days-hours:min:sec
#SBATCH --output=/work/NRSAAMR/Projects/SWAM/WGS/slurm/test.%j
#SBATCH --ntasks-per-node=8
#SBATCH --mem=20g

source activate swam-wgs

### Set env variables
home="/work/NRSAAMR/Projects/SWAM/WGS"
input="/work/NRSAAMR/Projects/SWAM/WGS/input"
out="/work/NRSAAMR/Projects/SWAM/WGS/output"

export TMPDIR="/work/NRSAAMR/Projects/SWAM/WGS/work"

### Import sample
base="SRR34965641"

### Clean reads with fastp
fastp -i "$input"/${base}_1.fastq -I "$input"/${base}_2.fastq -o "$out"/"$base"/${base}_R1.clean.fastq.gz -O "$out"/"$base"/${base}_R2.clean.fastq.gz --html /dev/null/ --json /dev/null/

### Assemble chromosome
mkdir -p "$out"/"$base"/chromosome
spades.py --isolate -1 "$input"/SRR34965641_1.fastq -2 "$input"/SRR34965641_2.fastq -o "$out"/"$base"/chromosome/"$base"

### Assemble plasmid(s)
mkdir -p "$out"/"$base"/plasmid
spades.py --plasmid -1 "$input"/SRR34965641_1.fastq -2 "$input"/SRR34965641_2.fastq -o "$out"/"$base"/plasmid/"$base"
