# SWAM-g
This pipeline is designed to process Illumina paired-end whole-genome sequencing (WGS) data for the U.S. EPA's Surface Water AMR Monitoring (SWAM) system. It was constructed to be a robust, species-agnostic pipeline for assembling and annotating WGS data for antimicrobial resistance (AMR) monitoring. SWAM-g (SWAM-genome) calls several canonical pipelines, including [NCBI's AMRFinderPlus](https://github.com/ncbi/amr/tree/master), [CGE's Resfinder](https://github.com/genomicepidemiology/resfinder) and [MobileElementFinder](https://pypi.org/project/MobileElementFinder/), [PHAC NML's MOB-suite](https://github.com/phac-nml/mob-suite) and [ECTyper](https://github.com/phac-nml/ecoli_serotyping), [TXXScan](https://github.com/macsy-models/TXSScan) for bacterial secretion systems, as well as [MLST](https://github.com/tseemann/mlst) for comprehensive assessment of AMR genotype, predicted phenotype, sequence type, and the presence and functionality of plasmids and assocaited mobile geneteic elements (MGEs). It leverages the circularity flags from the [SPAdes](https://github.com/ablab/spades/tree/main) wrapper, [Unicylcer](https://github.com/rrwick/Unicycler?tab=readme-ov-file#installation) for accurate plasmid reconstruction and typing using [MOB-recon](https://github.com/phac-nml/mob-suite). SWAM-g additionally calls [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) and [CheckM](https://github.com/Ecogenomics/CheckM) for accurate speciation and assembly QAQC. The taxonomic assignments from GTDB-tk are fed to both [AMRFinderPlus](https://github.com/ncbi/amr/tree/master) and [Resfinder](https://github.com/genomicepidemiology/resfinder) models of species-specific AMR profiling. If Salmonella or E. coli genomes are detected, SWAM-g will conditionally run [SeqSero2](https://github.com/denglab/SeqSero2/tree/master) for serotyping Salmonella or [ECTyper](https://github.com/phac-nml/ecoli_serotyping) for serotyping + pathotyping E. coli. All outputs are then collated into easy-to-use summary tables for downstream analysis.


<img width="1045" height="825" alt="SWAM-g_diagram drawio" src="https://github.com/user-attachments/assets/60767ac2-14a0-4cd3-b6f8-82c0c1f64d64" />


## Setup and Configuration
### Step 0. Install conda/mamba manager, install Snakemake and SRA-tools
`SWAM-g` was constructed using [Snakemake](https://github.com/snakemake/snakemake) and relies entirely on conda environments. 
Users should install either [miniforge](https://github.com/conda-forge/miniforge) or [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2).

Code snippet for installing miniforge3 in your home directory and reconfiguring your shell:
```bash
cd ~
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
rm Miniforge3-Linux-x86_64.sh
source ~/miniforge3/bin/activate

### Close and reopen terminal, then run
conda init --all

conda install snakemake sra-tools
```
### Step 1. Clone repo
```bash
cd /pathto/workspace
git clone https://github.com/davis-bc/SWAM-g
```
### Step 2. Configure slurm profile
`SWAM-g` was constructed and tested on an HPC cluster managed by a Slurm scheduler.
Minimum memory requirement of ~150GB due to GTDB-tk may restrict its use. 
To take advantage of Slurm job management, the following [slurm profile](https://github.com/davis-bc/EPA-SWAM-WGS/blob/main/config/slurm/config.yaml) must be edited to be configured to your particular HPC environment. 

Replace "###" with your credentials.

```bash
cluster: "sbatch --account=### --partition=### --time={resources.time} --output=/pathto/slurm/{rule}.%j --cpus-per-task={threads} --mem={resources.mem_mb}M"
jobs: 50
default-resources:
  - mem_mb=20000
  - threads=8
  - time='1-00:00:00'
```

More information on Snakemake configurations for different compute environments can be found in their documentation [Snakemake docs](https://snakemake.readthedocs.io/en/stable/)

### Step 3. Download test data, run the pipeline
To test the pipeline, first download example AMR-laden E.coli, S. enterica, and E. faecalis genomes:

```bash
mkdir -p /pathto/directory/input
cd /pathto/directory/input
fasterq-dump SRR30768419 SRR34965641 SRR7839461
```
`SWAM-g` takes a directory of paired-end fastqs as input.

The following is an example bash script for executing `SWAM-g` using a Snakemake driver. 

The `-j` flag dictates the max number of jobs (isolates) to be run in parallel.

```bash
#!/bin/bash
#SBATCH --account=###
#SBATCH --partition=###
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --job-name=snakemake_driver
#SBATCH --output=/pathto/workspace/slurm/snakemake.%j

cd /pathto/SWAM-g

input="/pathto/directory/input"
output="/pathto/directory/output"

snakemake --profile config/slurm/ \
          --config in_dir="$input" out_dir="$output" \
          --use-conda \
          --conda-frontend conda \
          -j 10 \
          --local-cores 1 

```






