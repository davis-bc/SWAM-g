# SWAM-g: A Comprehensive Whole-Genome Sequencing Workflow for Environmental AMR Surveillance
This pipeline is designed to process Illumina paired-end whole-genome sequencing (WGS) data for the U.S. EPA's Surface Water AMR Monitoring (SWAM) system. It was constructed to be a robust, species-agnostic pipeline for assembling and annotating WGS data for antimicrobial resistance (AMR) monitoring. `SWAM-g` (SWAM-genome) calls several canonical pipelines, including NCBI's [AMRFinderPlus](https://github.com/ncbi/amr/tree/master), CGE's [ResFinder](https://github.com/genomicepidemiology/resfinder) and [MobileElementFinder](https://pypi.org/project/MobileElementFinder/), PHAC NML's [MOB-suite](https://github.com/phac-nml/mob-suite) and [ECTyper](https://github.com/phac-nml/ecoli_serotyping), [TXSScan](https://github.com/macsy-models/TXSScan), as well as [MLST](https://github.com/tseemann/mlst) for comprehensive assessment of AMR genotype, predicted phenotype, sequence type, and the presence and functionality of plasmids, protein secretion systems, and associated mobile genetic elements (MGEs). It leverages the circularity flags from the [SPAdes](https://github.com/ablab/spades/tree/main) wrapper, [Unicycler](https://github.com/rrwick/Unicycler?tab=readme-ov-file#installation), for accurate plasmid reconstruction and typing using [MOB-recon](https://github.com/phac-nml/mob-suite). `SWAM-g` uses [MASH screen](https://github.com/marbl/Mash) against a curated GTDB reference sketch for rapid speciation and [CheckM2](https://github.com/chklovski/CheckM2) for assembly quality assessment. The taxonomic assignments from MASH are fed to both [AMRFinderPlus](https://github.com/ncbi/amr/tree/master) and [ResFinder](https://github.com/genomicepidemiology/resfinder) to enable species-specific AMR profiling. If *Salmonella* or *E. coli* genomes are detected, `SWAM-g` will conditionally run [SeqSero2](https://github.com/denglab/SeqSero2/tree/master) for serotyping *Salmonella* or [ECTyper](https://github.com/phac-nml/ecoli_serotyping) for serotyping + pathotyping *E. coli*. All outputs are then collated into easy-to-use summary tables for downstream analysis.


![SWAM-g pipeline](docs/pipeline_dag.png)


## Setup and Configuration
**Note:** The pipeline is currently designed for Linux systems. macOS users can attempt adaptation but might face Snakemake system-specific dependencies. Windows is not officially supported. If you encounter issues, refer to the [Snakemake issues page](https://github.com/snakemake/snakemake/issues).

The peak memory requirement is approximately **20 GB** for MASH taxonomy classification; all other rules use substantially less.

### Step 0. Install conda/mamba manager, Snakemake, and SRA-tools
`SWAM-g` was constructed using [Snakemake](https://github.com/snakemake/snakemake) and relies entirely on conda environments.
Users should install either [miniforge](https://github.com/conda-forge/miniforge) or [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2) if their shells are not already configured for Conda.

Code snippet for installing miniforge3 in your home directory and reconfiguring your shell:
```bash
cd ~
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -O installer.sh
bash installer.sh -b -p $HOME/miniforge3 && rm installer.sh
source ~/miniforge3/bin/activate
conda init && exec $SHELL
conda install "snakemake>=8" snakemake-executor-plugin-slurm sra-tools
```

> **HPC users:** `snakemake-executor-plugin-slurm` is required for the Slurm profile. If you are running on a local workstation only, it can be omitted.

### Step 1. Clone repo
Replace "/pathto/workspace" with where you want `SWAM-g` to live.

```bash
cd /pathto/workspace
git clone https://github.com/davis-bc/SWAM-g
```
### Step 2. Configure slurm profile and rule resources
`SWAM-g` was constructed and tested on an HPC cluster managed by a Slurm scheduler. Two Slurm profiles are provided — choose based on the number of samples in your run:

| Profile | Recommended for | Behavior |
|---------|----------------|----------|
| `config/slurm/large-batch` | **≥50 samples** | Annotation rules are batched 32 samples per Slurm job; Slurm receives the summed resources across all samples in a batch |
| `config/slurm/small-batch` | **<50 samples** | Each sample runs as its own individual Slurm job with per-job memory headroom |

Edit the `slurm_account` and `slurm_partition` fields in **both** profile configs to match your HPC environment. All resource allocation (`mem_mb`, `runtime`, `threads`) is controlled exclusively in the profile YAML — no edits to rule files are needed.

```yaml
default-resources:
  slurm_account: "your_account"
  slurm_partition: "your_partition"
```

More information on Snakemake configurations for different compute environments can be found in the [Snakemake docs](https://snakemake.readthedocs.io/en/stable/).

### Step 3. Tune resources (optional)
All per-rule resource allocations are defined in the profile config under `set-resources` and `set-threads`. Edit the appropriate profile file to adjust memory or runtime for any rule — no changes to `.smk` files are needed.

For example, to give `unicylcer` more RAM in the large-batch profile:
```yaml
# config/slurm/large-batch/config.yaml
set-resources:
  unicylcer:
    mem_mb: 200000
    runtime: "6-00:00:00"
```

### Step 4. Download test data, run the pipeline on HPC
To test the pipeline, first download example AMR-laden *E. coli*, *S. enterica*, and *E. faecalis* genomes using `fasterq-dump`. This will produce paired-end R1/R2 FASTQ files automatically:

```bash
mkdir -p /pathto/directory/input
cd /pathto/directory/input
fasterq-dump SRR30768419 SRR34965641 SRR7839461
```
`SWAM-g` takes a directory of paired-end FASTQs as input and a target directory for output. It currently cannot handle single-end Illumina or long-read datatypes.

The following is an example `run_swam-g.sh` driver script for executing `SWAM-g` via Slurm. Choose `large-batch` or `small-batch` based on the number of samples (see Step 2).

Replace "/pathto/" placeholders with appropriate paths.

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

# Use large-batch for ≥50 samples, small-batch for <50 samples
snakemake --profile config/slurm/large-batch \
          --config in_dir="$input" out_dir="$output" \
          --quiet

```
Then submit as:
```bash
sbatch run_swam-g.sh
```

### Step 5. Running on a local workstation
For running without Slurm (e.g., a laptop or single server), use the included `run_swam-g_local.sh` script. It uses the `config/local/` profile which caps CPU usage at 8 cores and RAM at 30 GB, and serializes the memory-intensive MASH classify step automatically.

```bash
bash run_swam-g_local.sh /path/to/input /path/to/output
```

If no arguments are provided, it defaults to `./input` and `./output` relative to the repository root.

---

## Outputs

### Filesystem overview

```
output/
├── SWAM-g_results.xlsx         # Primary deliverable — multi-sheet workbook (see below)
├── contig_map.csv              # Per-contig annotation map (see below)
└── data/
    ├── amrfinderplus/          # Per-sample AMR, stress, and virulence hits ({sample}.afp.tsv)
    ├── benchmarks/             # Snakemake per-rule wall-clock timing files
    ├── checkm2/                # Per-sample assembly quality reports (quality_report.tsv)
    ├── clean_reads/            # fastp-trimmed paired FASTQs ({sample}_R1/R2.clean.fastq.gz)
    ├── mash/                   # Per-sample MASH screen results + aggregated taxonomy table
    ├── mlst/                   # Multi-locus sequence typing across all samples (mlst.tsv)
    ├── mob-suite/              # Per-sample plasmid reconstruction and typing (mobtyper_results.txt)
    ├── mobileelementfinder/    # Per-sample mobile element annotations ({sample}.csv)
    ├── resfinder/              # Per-sample resistance gene hits and predicted phenotypes
    ├── serotype/
    │   ├── E.coli/             # ECTyper serotype + pathotype output (E. coli samples only)
    │   └── Salmonella/         # SeqSero2 antigenic profile + serotype (Salmonella samples only)
    ├── txsscan/                # Per-sample secretion system predictions (all_systems.tsv)
    └── unicycler/              # Assemblies (assembly.fasta), coverage TSVs, and protein FAAs
```

### `SWAM-g_results.xlsx`

A multi-sheet workbook collating all tool outputs. Each sheet can be used independently for downstream analysis.

| Sheet | Contents |
|-------|----------|
| `summary_out` | One row per sample: MASH species assignment, MLST sequence type, predicted serotype, AMRFinderPlus AMR genotype, ResFinder AMR genotype and predicted phenotype, PointFinder mutations, plasmid count / rep types / relaxase types, and TXSScan secretion systems present |
| `AMRFinderPlus` | Full per-hit AMRFinderPlus output including AMR, stress, and virulence elements with contig ID, coordinates, and strand |
| `assembly_QA` | CheckM2 completeness and contamination estimates plus mean read coverage; includes a `QA` pass/fail flag (N50 > 20 kb, total contigs < 500, coverage ≥ 30×) |
| `MOBrecon_summary` | MOB-suite plasmid typing: rep type, relaxase type, MPF type, predicted mobility, and predicted host range per plasmid cluster |
| `salmonella_serotype` | SeqSero2 antigenic profile and predicted serotype (populated for *Salmonella* samples only; empty otherwise) |
| `ecoli_serotype` | ECTyper serotype and pathotype predictions (populated for *E. coli* samples only; empty otherwise) |

### `contig_map.csv`

A flat per-row annotation table linking every detected genetic element to its genomic context. One row per gene/element hit per contig. Columns:

| Column | Description |
|--------|-------------|
| `Sample` | Sample identifier |
| `contig_id` | Full Unicycler contig header (includes length and depth) |
| `circularity_status` | `complete` (circular) or `incomplete` (linear) as called by Unicycler |
| `molecule_type` | `chromosome` or `plasmid` as assigned by MOB-recon |
| `primary_cluster_id` | MOB-suite plasmid cluster ID (`-` for chromosomal contigs) |
| `predicted_mobility` | MOB-suite mobility prediction for the plasmid cluster |
| `predicted_host_range_overall_name` | Predicted host range taxon from MOB-suite |
| `gene` | Gene or element name from AMRFinderPlus, MobileElementFinder, or TXSScan |
| `type` | Element category (e.g. `AMR`, `VIRULENCE`, `T4SS`, `miniature inverted repeat`) |
| `start` / `end` | Genomic coordinates on the contig (bp) |
| `strand` | `positive` or `negative` |

---

## Software Versions

Key tools are pinned to specific versions in the conda environment files (`workflow/envs/*.yaml`) to ensure reproducibility. Tools installed via `pip` or marked with `>=` use the latest compatible release at environment creation time.

| Tool | Version | Notes |
|------|---------|-------|
| AMRFinderPlus | 4.2.7 (latest) | Exact pin |
| ECTyper | 2.0.0 (latest) | Exact pin |
| MacSyFinder (TXSScan) | 2.1.6 (latest) | Exact pin |
| MLST | 2.25.0 | Exact pin |
| CheckM2 | 1.1.0 (latest) | Minimum version |
| Unicycler | ≥0.5.0 | Minimum version |
| SPAdes | ≥3.14.0 | Minimum version |
| MASH | ≥2.3 | Minimum version |
| fastp | latest | Unpinned |
| Prodigal | latest | Unpinned |
| ResFinder | latest | Unpinned |
| MOB-suite | latest (pip) | Unpinned |
| MobileElementFinder | latest (pip) | Unpinned |
| SeqSero2 | git HEAD | Installed from source |






