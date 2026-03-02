# SWAM-g Copilot Instructions

## Session continuity

**At the start of every session:** read `.copilot/session-log.md` for context on
recent work, known bugs, and pending next steps. This file is local to the
workstation (gitignored) — do not commit it.

**When wrapping up** (or when asked to "log progress" / "log our progress"):
append a dated `## YYYY-MM-DD` entry to `.copilot/session-log.md` summarising
what was done, what files changed, and what remains outstanding.

---

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
            ├─ MASH screen (taxonomy) → mash_taxonomy.tsv
            │    └─ feeds species into AMRFinder + ResFinder via
            │       get_amrfinder_organism / get_resfinder_species
            ├─ CheckM2 (assembly QC)
            ├─ AMRFinderPlus (AMR/stress/virulence)
            ├─ ResFinder + PointFinder + DisinFinder (phenotype prediction)
            ├─ MobileElementFinder (MGEs)
            ├─ MLST (sequence typing)
            ├─ TXSScan/MacSyFinder (secretion systems, via Prodigal proteins)
            ├─ SeqSero2 (Salmonella only)
            ├─ ECTyper (E. coli only — rule currently commented out in serotype.smk)
            └─ data_summary.R → SWAM-g_results.xlsx + contig_map.csv
```

### Key design patterns

- **`in_dir` / `out_dir` are passed at runtime** via `--config`; there is no hardcoded `configfile:` directive in `Snakefile` (it is commented out). All paths are constructed with `os.path.join(output_dir, ...)`.
- **Sample discovery** happens in `workflow/rules/common.smk`: R1 files are glob-matched from `in_dir`, and R2 pairs are resolved by name substitution. Supported naming: `*R1*`, `*_1.fastq*`. The regex in `extract_sample_name()` strips read suffixes to produce the sample name.
- **Database init rules** (`checkm_init`, `mob_init`, `afp_init`, `res_init`, `txsscan_init`, `ectyper_init`) are declared as `localrules` — they run on the head node and create sentinel `.touch` files. Downstream rules depend on these sentinel files as inputs. Databases live in a `dbs/` directory relative to the working directory where Snakemake is invoked. `ectyper_init` also downloads the `EnteroRef_GTDBSketch_20231003_V2.msh` MASH sketch used by `mash_classify`.
- **Species-aware annotation**: MASH screen output (`mash_taxonomy.tsv`) is parsed by two Python scripts (`get_amrfinder_organism.py`, `get_resfinder_species.py`) to produce per-sample organism/species maps. These are loaded into Python dicts in `common.smk` and injected as `params` into AMRFinder and ResFinder rules via `lambda wildcards:` functions. Note: `common.smk` loads `amrfinder_species.tsv` but the rule outputs `amrfinder_organism.tsv` — a known mismatch that causes all samples to default to `Escherichia`.
- **Plasmid handling uses marker files**: Because plasmid count per sample is variable (0–N), rules use hidden `.done` marker files (e.g., `.{sample}.plasmids.done`) to signal completion when no plasmid FASTAs are produced. Shell blocks always `touch {output}` at the end to satisfy Snakemake.
- **MOB-recon exploits Unicycler circularity flags**: `mob_recon --unicycler_contigs` is intentional — Unicycler annotates circular contigs in FASTA headers, and MOB-suite uses this for more accurate plasmid reconstruction.
- **Final summarization is in R**: `workflow/scripts/data_summary.R` runs under the `Renv.yaml` conda env. It uses the `snakemake@input` / `snakemake@output` S4 object to access named inputs. It joins all tool outputs and writes a multi-sheet `.xlsx` (via `writexl`) plus a `contig_map.csv`.

## Configuration

- `config/config.yaml` — hardcoded DB paths for the EPA HPC environment (not used at runtime; `configfile:` is commented out in `Snakefile`)
- `config/slurm/config.yaml` — Slurm profile: account, partition, log paths, per-job resource defaults. **Must be edited** before first use.
- `config/local/config.yaml` — Local workstation profile (8 cores, 32 GB RAM, greedy scheduler). Run via `bash run_swam-g_local.sh`.

## Resource tuning

Each `.smk` rule has a `resources:` block (`mem_mb`, `threads`, `time`). Adjust these directly in the relevant rule file. Notable defaults:
- `fastp` / `unicylcer`: 10 GB / 150 GB RAM, unicylcer uses 32 threads (capped to 8 via local profile)
- `mash_classify`: declared 8 GB but overridden to 20 GB in both profiles to prevent OOM when loading the 944 MB EnteroRef sketch
- `checkm2`: 2 GB per sample (runs sequentially with `threads: 1`)

## Conda environments

Each rule specifies its own conda env via `conda: "../envs/<name>.yaml"`. All envs are in `workflow/envs/`. Rules without a `conda:` directive (like `get_resfinder_species`, `get_amrfinder_organism`) run in the base environment — these only use the Python stdlib.

## Input requirements

- Paired-end Illumina FASTQs only (no single-end, no long-read)
- Files must match glob patterns: `*R1*.fastq*` / `*R2*.fastq*` or `*_1.fastq*` / `*_2.fastq*`
- All FASTQs for a run must be in a single flat directory (`in_dir`)
