# SWAM-g
This pipeline is designed to process Illumina paired-end whole-genome sequencing (WGS) data for the U.S. EPA's Surface Water AMR Monitoring (SWAM) system. It was constructed to be a robust, species-agnostic pipeline for assembling and annotating WGS data for antimicrobial resistance (AMR) monitoring. SWAM-g (SWAM-genome) calls several canonical pipelines, including [NCBI's AMRFinderPlus](https://github.com/ncbi/amr/tree/master), [CGE's Resfinder](https://github.com/genomicepidemiology/resfinder) and [MobileElementFinder](https://pypi.org/project/MobileElementFinder/), [PHAC NML's MOB-suite](https://github.com/phac-nml/mob-suite) and [ECTyper](https://github.com/phac-nml/ecoli_serotyping), as well as [MLST](https://github.com/tseemann/mlst) for comprehensive assessment of AMR genotype, predicted phenotype, sequence type, and the presence and functionality of plasmids and assocaited mobile geneteic elements (MGEs). It leverages the circularity flags from the [SPAdes](https://github.com/ablab/spades/tree/main) wrapper, [Unicylcer](https://github.com/rrwick/Unicycler?tab=readme-ov-file#installation) for accurate plasmid reconstruction and typing using [MOB-recon](https://github.com/phac-nml/mob-suite). Albeit cumbersome, SWAM-g additionally calls [GTDB-tk](https://ecogenomics.github.io/GTDBTk/) and [CheckM](https://github.com/Ecogenomics/CheckM) for accurate speciation and assembly QAQC. If Salmonella or E. coli genomes are detected, SWAM-g will conditionally run [SeqSero2](https://github.com/denglab/SeqSero2/tree/master) for serotyping Salmonella or [ECTyper](https://github.com/phac-nml/ecoli_serotyping) for serotyping + pathotyping E. coli. All outputs are then collated into easy-to-use summary tables for downstream analysis.

<img width="952" height="402" alt="EPA-SWAM-WGS drawio" src="https://github.com/user-attachments/assets/e4fc8fcc-2b67-4ec9-b420-e26670774c20" />

## Setup and Configuration
SWAM-g was constructed using [Snakemake](https://github.com/snakemake/snakemake) and relies entirely on conda environments. Because SWAM-g strings several existing pipelines together, it will require some effort to prepare the appropriate package manager, databases, and config files.

### Step 0. Install conda/mamba manager, install Snakemake and SRA-tools
Recomend using [miniforge](https://github.com/conda-forge/miniforge) over [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2)  

Code snippet for installing miniforge3 in your home directory:
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
git clone https://github.com/davis-bc/EPA-SWAM-WGS
```
### Step 2. Download and configure databases
A suite of databases are used, some that will require downloading beforehand including ([Resfinder](https://github.com/genomicepidemiology/resfinder), [GTDB](https://ecogenomics.github.io/GTDBTk/), [CheckM](https://github.com/Ecogenomics/CheckM)). 

[AMRFinderPlus](https://github.com/ncbi/amr/tree/master) and [MOB-suite](https://github.com/phac-nml/mob-suite) are configured to auto-update within the pipeline.

Careful, this will take several hours to download and best to be run as a batch job.
```bash
cd /pathto/database/storage

git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/
wget https://data.ace.uq.edu.au/public/CheckM_databases
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
```
Once downloaded (or if already downloaded elsewhere), provide the absolute paths in the [config.yaml](https://github.com/davis-bc/EPA-SWAM-WGS/blob/main/config/config.yaml)

```bash
gtdbtk_db: /pathto/release226
checkm_db: /pathto/checkm_data
res_db:    /pathto/resfinder_db
pt_db:     /pathto/pointfinder_db
dis_db:    /pathto/disinfinder_db
```
### Step 3. Configure slurm profile
SWAM-g minimizes runtimes through job parallelization on an HPC cluster managed by a Slurm scheduler.
The following [slurm profile](https://github.com/davis-bc/EPA-SWAM-WGS/blob/main/config/slurm/config.yaml) must be edited to be configured to your particular HPC environment.

```bash
cluster: "sbatch --account=### --partition=### --time={resources.time} --output=/pathto/slurm/{rule}.%j --cpus-per-task={threads} --mem={resources.mem_mb}M"
jobs: 50
default-resources:
  - mem_mb=20000
  - threads=8
  - time='1-00:00:00'
```
More information on Snakemake configurations for different compute environments can be found in their documentation [Snakemake docs](https://snakemake.readthedocs.io/en/stable/)

### Step 4. Download test data, run the pipeline
To test the pipeline, first download example AMR-laden E.coli, S. enterica, and E. faecalis genomes:

```bash
mkdir -p /pathto/workspace/input
cd /pathto/workspace/input
fasterq-dump SRR30768419 SRR34965641 SRR7839461
```

The following is an example bash script for executing the pipeline as a batch job. SWAM-g uses process grouping for more efficient slurm management.

```bash
#!/bin/bash
#SBATCH --account=###
#SBATCH --partition=###
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --job-name=snakemake_driver
#SBATCH --output=/pathto/workspace/slurm/snakemake.%j

cd /pathto/EPA-SWAM-WGS

input="/pathto/workspace/input"
output="/pathto/workspace/output"

snakemake --profile config/slurm/ \
          --config in_dir="$input" out_dir="$output" \
          --use-conda \
          --conda-frontend conda \
          -j 10 \
          --local-cores 1 \
          --group-components group1=3 group2=3 group3=3

### Optionally generate a snakemake report upon completion
#snakemake --report report.html --config in_dir="$input" out_dir="$output" --configfile config/config.yaml
```
The `-j` flag denotes the total number of jobs that can be run in parallel.

The `--group-components` flag denotes the number samples to be grouped into singular jobs. For example, with the 3 test samples, these values are set to 3.

|   Group |    Process                    |
|---------|-------------------------------|
| group1  | calculate coverage (minimap2) |
| group2  | mob-suite                     |
| group3  | AMRFinderPlus                 |





