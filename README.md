# EPA SWAM WGS
This pipeline is designed for processing whole-genome sequencing (WGS) data for the U.S. EPA's Surface Water AMR Monitoring (SWAM) system. It was constructed to be a robust, species-agnostic pipeline for assembling and annotating WGS data for antimicrobial resistance (AMR) monitoring. EPA-SWAM-WGS calls several canonical pipelines, including [NCBI's AMRFinderPlus](https://github.com/ncbi/amr/tree/master), [CGE's Resfinder](https://github.com/genomicepidemiology/resfinder), [PHAC NML's MOB-suite and ECTyper](https://github.com/phac-nml), and [MLST](https://github.com/tseemann/mlst) for comprehensive assessment of AMR genotype, predicted phenotype, sequence type, and the presence and functionality of plasmids. It leverages the [SPAdes](https://github.com/ablab/spades/tree/main) wrapper, [Unicylcer](https://github.com/rrwick/Unicycler?tab=readme-ov-file#installation)'s circularity flags for accurate plasmid reconstruction and typing using [MOB-recon](https://github.com/phac-nml/mob-suite). If Salmonella or E. coli genomes are detected, EPA-SWAM-WGS will conditionally run [SeqSero2](https://github.com/denglab/SeqSero2/tree/master) for serotyping Salmonella or [ECTyper](https://github.com/phac-nml/ecoli_serotyping) for serotyping+pathotyping E. coli.

<img width="952" height="402" alt="EPA-SWAM-WGS drawio" src="https://github.com/user-attachments/assets/e4fc8fcc-2b67-4ec9-b420-e26670774c20" />

## Setup and Configuration
This pipeline was constructed using [Snakemake](https://github.com/snakemake/snakemake) and relies entirely on conda environments. Because EPA-SWAM-WGS strings several existing pipelines together, it will require some effort to prepare the appropriate databases and config files. This setup and configuration guide assumes that [miniforge](https://github.com/conda-forge/miniforge) or an equivalent conda/mamba manager is already installed and configured on an HPC system.

### Step 1. Clone repo
```bash
cd /homebase
git clone https://github.com/davis-bc/EPA-SWAM-WGS
```
### Step 2. Download and configure databases
A large collection of databases are used, some that will require downloading beforehand ([Resfinder](https://github.com/genomicepidemiology/resfinder), [GTDB](https://ecogenomics.github.io/GTDBTk/), [CheckM](https://github.com/Ecogenomics/CheckM)). 

[AMRFinderPlus](https://github.com/ncbi/amr/tree/master) and [MOB-suite](https://github.com/phac-nml/mob-suite) are configured to auto-update within the pipeline.

Careful, this will take several hours to download and best to be run as a batch job.
```bash
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
This pipeline was constructed to minimize runtimes through job parallelization on an HPC cluster managed by a Slurm scheduler.
The following [slurm profile](https://github.com/davis-bc/EPA-SWAM-WGS/blob/main/config/slurm/config.yaml) must be edited to be configured to your particular HPC environment.
```bash
cluster: "sbatch --account=nrsaamr --partition=ord --time={resources.time} --output=/work/NRSAAMR/Projects/SWAM/WGS/slurm/{rule}.%j --ntasks-per-node={threads} --mem={resources.mem_mb}M"
jobs: 50
default-resources:
  - mem_mb=20000
  - threads=8
  - time='1-00:00:00'
```
### Step 4. Run the pipeline
EPA-SWAM-WGS was designed to be executed on a directory of paired-end fastq's and is agnostic to common filename syntax and gzipping.
For assembly, all detected samples will be submitted as individual jobs to minimize runtimes. 
Because GTDB-tk and CheckM are executed, the more genomes that are provided, the more efficient the runtimes will be.




