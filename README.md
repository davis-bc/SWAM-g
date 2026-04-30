```                                                                                                                                                                                                                                                           
                                         _______          __     __  __             
                                        / ____\ \        / /\   |  \/  |            
                                       | (___  \ \  /\  / /  \  | \  / |______ __ _ 
                                        \___ \  \ \/  \/ / /\ \ | |\/| |______/ _` |
                                        ____) |  \  /\  / ____ \| |  | |     | (_| |
                                       |_____/    \/  \/_/    \_\_|  |_|      \__, |
                                                                               __/ |
                                                                              |___/ 
```
`SWAM-g` is a Snakemake pipeline for Illumina bacterial whole-genome sequencing data. It combines read QC, assembly, taxonomy, AMR profiling, plasmid and mobile-element screening, MLST, conditional serotyping, and optional Pathogen Detection metadata enrichment into one reporting workflow.

`SWAM-g` takes a directory of paired-end FASTQs as input and a target directory for output. It currently cannot handle single-end Illumina or long-read datatypes.



<img width="976" height="755" alt="SWAM-g_diagram-V3 drawio" src="https://github.com/user-attachments/assets/1ad338d2-5b63-4f1b-8fe3-227384eb3831" />



## Setup and Configuration

**Platform support:**
- **Linux:** fully tested and recommended.
- **macOS:** supported via conda; Apple Silicon (M-series) users may need [Rosetta 2](https://support.apple.com/en-us/102527) for some tools.
- **Windows:** use [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install) with Ubuntu — install conda inside WSL and follow the Linux instructions below.

> **Databases:** Reference databases and model bundles are bootstrapped by Snakemake local init rules into the repo-level `dbs/` directory. On HPC, the first bootstrap or any later refresh still needs to run from an internet-connected login/head node; downstream compute-node jobs then reuse the staged `dbs/` contents fully offline.

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
### Step 2. Configure your Slurm credentials
Edit the `slurm_account` and `slurm_partition` fields in `config/slurm/config.yaml`.

```yaml
default-resources:
  slurm_account: "your_account"
  slurm_partition: "your_partition"
```

More information on Snakemake configurations for different compute environments can be found in the [Snakemake docs](https://snakemake.readthedocs.io/en/stable/).

### Step 3. Download test data
To test the pipeline, download example AMR-laden *E. coli*, *S. enterica*, and *E. faecalis* genomes using `fasterq-dump`. This will produce paired-end R1/R2 FASTQ files automatically:

```bash
mkdir -p /pathto/directory/input
cd /pathto/directory/input
fasterq-dump SRR30768419 SRR34965641 SRR7839461
```

### Step 4. Run the pipeline on an HPC cluster
Executing `bash run_swam-g.sh input output` on a login node is not recommended, especially for busy clusters.

Best to execute within a driver script `submit-swam-g.sh`:

```bash
#!/bin/bash
#SBATCH --account=###
#SBATCH --partition=###
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --job-name=swam-g_driver
#SBATCH --output=/pathto/workspace/slurm/swam-g.%j

swam_g="/pathto/SWAM-g/run_swam-g.sh"
input="/pathto/directory/input"
output="/pathto/directory/output"

bash "$swam_g" "$input" "$output"

```
Then submit as:
```bash
sbatch submit-swam-g.sh
```

## Notes

### Run outputs and troubleshooting
The HPC entrypoint `run_swam-g.sh` writes:

- `output/logs/run_status/run_report.txt`
- `output/logs/run_status/sample_rule_status.tsv`

If startup logs show a warning like `Mismatch in number of R1 (...) and R2 (...) files`, `SWAM-g` will continue with the paired subset only. Fix or remove unmatched FASTQs before a production run.

### Resource tuning
The default Slurm profile in `config/slurm/config.yaml` has been tuned from benchmarking data.

For `unicylcer`, assembly retries scale automatically in `workflow/rules/assemble.smk`:

1. attempt 1: `40000 MB`, `4h`
2. attempt 2: `60000 MB`, `6h`
3. attempt 3: `80000 MB`, `8h`

If assembly still fails, raise retry tiers before broadly increasing other rule allocations.

### Salmonella serotyping behavior
Salmonella serotyping uses two complementary methods:

1. **SeqSero2 allele mode** on the `fastp`-cleaned paired-end reads
2. **SISTR** on the draft assembly

Both methods are gated by the run-level MASH taxonomy table. If a sample is not classified as `g__Salmonella`, SWAM-g writes the expected placeholder TSV and skips the serotyping runtime for that sample.

### Running on a local workstation
For running without Slurm (e.g., a workstation, laptop, or VM), use the included `run_swam-g_local.sh` script. It uses the `config/local/` profile, which caps CPU usage at 8 cores and total RAM at 30 GB, and serializes the memory-intensive MASH classify step automatically. This is a good way to validate the pipeline with the test data before running a full HPC batch.

```bash
bash run_swam-g_local.sh /path/to/input /path/to/output
```

If no arguments are provided, it defaults to `./input` and `./output` relative to the repository root.

**macOS:** Run as above — conda handles all dependencies natively.

**Windows (WSL2):** Open a WSL2 Ubuntu terminal, navigate to the cloned repository, and run the same command. Ensure conda is installed inside WSL (not the Windows host) before proceeding.

### Pathogen Detection SNP-cluster enrichment
`SWAM-g` appends NCBI Pathogen Detection (PD) isolate metadata and SNP-cluster assignments to the final workbook by default. This is a **post-processing enrichment step** only — it does **not** run local SNP clustering and it does **not** block the main pipeline if PD lookups are unavailable.

By default, `SWAM-g` now queries the **current public PD FTP `latest_snps` release** at runtime (`pd_backend=ftp`). The lookup flow is:

1. resolve `SRR -> BioSample` through NCBI SRA metadata when needed
2. map the sample to a PD organism group using MASH species (or SRA species as fallback)
3. find the current PD isolate record and SNP cluster from the live FTP release
   - some taxgroups publish `Clusters/` and `Exceptions/` without a `Metadata/` directory; SWAM-g now falls back to `cluster_list.tsv` where possible instead of failing the run
4. return up to `pd_comparator_limit` cluster comparators
   - exact SNP-ranked neighbors when the public `SNP_distances.tsv` is small enough to scan
   - otherwise 10 same-cluster comparators from the live cluster membership table

Optional inputs:

- a per-sample metadata file (`pd_sample_metadata_tsv`) if your sample names are not already SRR or BioSample accessions
- `pd_backend=table` plus a PD isolates export (`pd_isolates_tsv`) and optional exceptions export (`pd_exceptions_tsv`) if you want to use local tables instead of the live FTP backend

Supported metadata columns are matched case-insensitively and may include:

- `Sample` (required if `pd_sample_metadata_tsv` is provided)
- `SRR` / `srr_acc`
- `BioSample` / `biosample_acc`
- `target_acc`
- `asm_acc`
- `scientific_name` / `species`

If `Sample` itself looks like an SRR accession (for example `SRR30768419`), `SWAM-g` will try to resolve the corresponding BioSample automatically through NCBI SRA metadata before joining against the PD table.

For the bundled driver scripts, use explicit SWAM-g flags instead of appending raw `pd_lookup=true` tokens after the script arguments:

```bash
bash run_swam-g.sh "$input" "$output" --pd-lookup=false

bash run_swam-g.sh "$input" "$output" \
    --pd-backend=table \
    --pd-isolates-tsv="/path/to/pd_isolates.tsv" \
    --pd-exceptions-tsv="/path/to/pd_isolate_exceptions.tsv"

bash run_swam-g_local.sh /path/to/input /path/to/output \
    --pd-sample-metadata-tsv="/path/to/sample_metadata.tsv"
```

If you are calling `snakemake` directly, the equivalent config looks like:

```bash
snakemake --profile config/local/ \
          --config \
            in_dir="$INPUT" \
            out_dir="$OUTPUT" \
            pd_backend=ftp \
            pd_comparator_limit=10 \
            pd_sample_metadata_tsv="/path/to/sample_metadata.tsv"
```

If you prefer local table lookups instead:

```bash
snakemake --profile config/local/ \
          --config \
            in_dir="$INPUT" \
            out_dir="$OUTPUT" \
            pd_backend=table \
            pd_isolates_tsv="/path/to/pd_isolates.tsv" \
            pd_exceptions_tsv="/path/to/pd_isolate_exceptions.tsv"
```

The lookup outputs are written to:

- `output/data/pd/pd_isolate_metadata.tsv`
- `output/data/pd/pd_cluster_comparators.tsv`

Reported statuses include `FOUND`, `FOUND_NO_CLUSTER`, `NOT_FOUND`, `QC_EXCEPTION`, `NO_ACCESSION`, `LOOKUP_ERROR`, `CONFIG_ERROR`, `UNSUPPORTED_ORGANISM`, and `DISABLED`.

Because the public PD distance tables can be extremely large for some taxgroups, exact pairwise comparator ranking is not always practical. In those cases, `SWAM-g` falls back automatically to same-cluster comparators from the current live cluster membership table and records that choice in `pd_lookup_note` / `PD_Comparator_Mode`.

---

## Outputs

### Filesystem overview

```
output/
├── mashtree.nwk               # Quick assembly relationship sketch from mashtree
├── SWAM-g_results.xlsx         # Primary deliverable — multi-sheet workbook (see below)
├── pd_isolate_metadata.xlsx    # Optional PD workbook containing isolate metadata + comparators
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
    ├── pd/                     # Optional Pathogen Detection metadata + comparator tables
    ├── resfinder/              # Per-sample resistance gene hits and predicted phenotypes
    ├── serotype/
    │   ├── E.coli/             # ECTyper serotype + pathotype output (E. coli samples only)
    │   └── Salmonella/         # SeqSero2 allele-mode + SISTR serotyping outputs (MASH-called Salmonella only)
    │       └── {sample}/
    │           ├── SeqSero_result.tsv
    │           └── sistr.tsv
    ├── txsscan/                # Per-sample secretion system predictions (all_systems.tsv)
    └── unicycler/              # Assemblies (assembly.fasta), coverage TSVs, and protein FAAs
```

### `SWAM-g_results.xlsx`

A multi-sheet workbook collating all tool outputs. Each sheet can be used independently for downstream analysis.

| Sheet | Contents |
|-------|----------|
| `summary_out` | One row per sample: MASH species assignment, MLST scheme, MLST sequence type, raw Salmonella serotype calls (`SeqSero2_Serotype`, `SISTR_Serovar`), `Ecoli_serotype` from ECTyper when applicable, AMRFinderPlus AMR genotype, ResFinder AMR genotype and predicted phenotype, PointFinder mutations, plasmid count / rep types / relaxase types, and deduplicated TXSScan secretion-system models present |
| `AMRFinderPlus` | Full per-hit AMRFinderPlus output including AMR, stress, and virulence elements with contig ID, coordinates, and strand |
| `assembly_QA` | CheckM2 completeness and contamination estimates plus mean read coverage; includes a `QA` pass/fail flag (N50 > 20 kb, total contigs < 500, coverage ≥ 30×) |
| `MOBrecon_summary` | MOB-suite plasmid typing: rep type, relaxase type, MPF type, predicted mobility, and predicted host range per plasmid cluster |
| `salmonella_serotype` | Combined Salmonella serotyping table retaining the raw SeqSero2 and SISTR outputs side-by-side for MASH-called *Salmonella* samples only |
| `ecoli_serotype` | ECTyper serotype and pathotype predictions for samples classified as *Escherichia coli* by the pipeline MASH taxonomy |

### `pd_isolate_metadata.xlsx`

Optional PD workbook written separately from the main results workbook.

| Sheet | Contents |
|-------|----------|
| `pd_isolate_metadata` | Optional PD enrichment table containing resolved SRR/BioSample accessions, PD target/assembly accessions, SNP cluster, source metadata, comparator mode/count, lookup source, and status for each sample |
| `pd_cluster_comparators` | Up to `pd_comparator_limit` live PD comparators per queried isolate, including rank, SNP distance when available, accession fields, and source environment metadata for rough source-tracking review |

### `mashtree.nwk`

A run-level Newick tree sketched from the symlinked assembly FASTAs in `output/data/unicycler/batch/{sample}.fasta` using `mashtree`. Tip labels follow the pipeline sample names.

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
| MacSyFinder (TXSScan) | 2.1.6 | Installed via `pip` inside a Python 3.12 conda env with `hmmer`/`prodigal` from conda; required because TXSScan 1.1.4 models use 2.1 grammar that MacSyFinder 2.1.4 rejects; models are staged under `dbs/macsyfinder/models` |
| MLST | 2.25.0 | Exact pin |
| CheckM2 | 1.1.0 (latest) | Minimum version |
| Unicycler | ≥0.5.0 | Minimum version |
| SPAdes | ≥3.14.0 | Minimum version |
| MASH | ≥2.3 | Minimum version |
| fastp | latest | Unpinned |
| Prodigal | latest | Unpinned |
| ResFinder | latest | Unpinned |
| MOB-suite | latest (pip) | Unpinned; database is staged under `dbs/mob_suite` and passed explicitly to `mob_recon` for offline HPC execution |
| MobileElementFinder | 1.1.2 (pip) | Exact pin with `mgedb==1.1.1`; MEF inputs are sanitized to strip UTF-8 BOMs and trim FASTA headers before BLAST |
| SeqSero2 | 1.3.1 | Exact pin from Bioconda; environment currently solved with `python=3.7`, `biopython=1.73` |
| SISTR | ≥1.1.3 | Minimum version; environment currently pinned to `python=3.11`, `setuptools<81` for CLI compatibility |
