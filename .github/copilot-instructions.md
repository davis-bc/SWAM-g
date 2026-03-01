# SWAM-g Copilot Instructions

## What this project is

**SWAM-g** is a Snakemake-based bioinformatics pipeline for whole-genome sequencing (WGS) data, built for the U.S. EPA's Surface Water AMR Monitoring (SWAM) system. It assembles and annotates paired-end Illumina reads for AMR surveillance. It runs entirely on HPC (Slurm) using conda environments — there is no local dev server and no test suite.

## Running the pipeline

```bash
# Standard execution via Slurm profile
snakemake --profile config/slurm/ \
          --config in_dir="$input" out_dir="$output" \
          --use-conda \
          --conda-frontend conda \
          -j 10 \
          --local-cores 1 \
          --quiet

# Dry-run to validate rules/DAG without executing
snakemake --profile config/slurm/ \
          --config in_dir="$input" out_dir="$output" \
          --use-conda -n

# Force re-run of a specific rule (e.g., after editing summarize.smk)
snakemake --profile config/slurm/ \
          --config in_dir="$input" out_dir="$output" \
          --use-conda \
          --forcerun summarize_results
```

The pipeline is designed to be submitted as an SBATCH driver job (`sbatch run_swam-g.sh`), which then dispatches per-sample jobs through Slurm.

## Architecture

### DAG flow

```
Input: paired-end FASTQs (R1/R2)
  └─ fastp (QC trim) → Unicycler (assembly)
       └─ MOB-recon → chromosome + plasmid FASTAs
            ├─ GTDB-tk (taxonomy) → feeds species into AMRFinder + Resfinder
            ├─ CheckM (assembly QC)
            ├─ AMRFinderPlus (AMR/stress/virulence — chromosome + plasmids)
            ├─ Resfinder + PointFinder + DisinFinder (phenotype prediction)
            ├─ MobileElementFinder (MGEs)
            ├─ MLST (sequence typing)
            ├─ TXSScan/MacSyFinder (secretion systems, via Prodigal proteins)
            ├─ SeqSero2 (Salmonella only — conditional)
            ├─ ECTyper (E. coli only — conditional)
            └─ FastANI (pairwise ANI matrix)
                 └─ data_summary.R → SWAM-g_results.xlsx + contig_map.csv
```

### Key design patterns

- **`in_dir` / `out_dir` are passed at runtime** via `--config`; there is no hardcoded `configfile:` directive in `Snakefile` (it is commented out). All paths are constructed with `os.path.join(output_dir, ...)`.
- **Sample discovery** happens in `workflow/rules/common.smk`: R1 files are glob-matched from `in_dir`, and R2 pairs are resolved by name substitution. Supported naming: `*R1*`, `*_1.fastq*`. The regex in `extract_sample_name()` strips read suffixes to produce the sample name.
- **Database init rules** (`gtdb_init`, `checkm_init`, `mob_init`, `afp_init`, `res_init`, `txsscan_init`, `ectyper_init`) are declared as `localrules` — they run on the head node and create sentinel `.touch` files. Downstream rules depend on these sentinel files as inputs. Databases live in a `dbs/` directory relative to the working directory where Snakemake is invoked.
- **Species-aware annotation**: GTDB-tk output (`gtdbtk.bac120.summary.tsv`) is parsed by two Python scripts (`get_amrfinder_organism.py`, `get_resfinder_species.py`) to produce per-sample organism/species maps. These are loaded into Python dicts in `common.smk` and injected as `params` into AMRFinder and Resfinder rules via `lambda wildcards:` functions.
- **Plasmid handling uses marker files**: Because plasmid count per sample is variable (0–N), rules use hidden `.done` marker files (e.g., `.{sample}.plasmids.done`) to signal completion when no plasmid FASTAs are produced. Shell blocks always `touch {output}` at the end to satisfy Snakemake.
- **MOB-recon exploits Unicycler circularity flags**: `mob_recon --unicycler_contigs` is intentional — Unicycler annotates circular contigs in FASTA headers, and MOB-suite uses this for more accurate plasmid reconstruction.
- **Final summarization is in R**: `workflow/scripts/data_summary.R` runs under the `Renv.yaml` conda env. It uses the `snakemake@input` / `snakemake@output` S4 object to access named inputs. It joins all tool outputs and writes a multi-sheet `.xlsx` (via `writexl`) plus a `contig_map.csv`.

## Configuration

- `config/config.yaml` — hardcoded DB paths for the EPA HPC environment (not used at runtime; `configfile:` is commented out in `Snakefile`)
- `config/slurm/config.yaml` — Slurm profile: account, partition, log paths, per-job resource defaults. **Must be edited** before first use.

## Resource tuning

Each `.smk` rule has a `resources:` block (`mem_mb`, `threads`, `time`). Adjust these directly in the relevant rule file. Notable defaults:
- `fastp_and_unicylcer`: 50GB RAM, 32 threads, 1 day
- `gtdbtk`: 150GB RAM (minimum system requirement)
- `checkm`: 100GB RAM

## Conda environments

Each rule specifies its own conda env via `conda: "../envs/<name>.yaml"`. All envs are in `workflow/envs/`. Rules without a `conda:` directive (like `get_resfinder_species`, `get_amrfinder_organism`) run in the base environment — these only use the Python stdlib.

## Input requirements

- Paired-end Illumina FASTQs only (no single-end, no long-read)
- Files must match glob patterns: `*R1*.fastq*` / `*R2*.fastq*` or `*_1.fastq*` / `*_2.fastq*`
- All FASTQs for a run must be in a single flat directory (`in_dir`)
