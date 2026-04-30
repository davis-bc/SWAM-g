import os
import glob
import pandas as pd
import re
import csv


# ----------------------------------------------
#   Generalized Sample Parsing from FASTQ Files
# ----------------------------------------------

# Load config variables
input_dir  = config["in_dir"]
output_dir = config["out_dir"]
mash_taxonomy_file = os.path.join(output_dir, "data", "mash", "mash_taxonomy.tsv")

def config_bool(key, default=False):
    raw = config.get(key, default)
    return raw if isinstance(raw, bool) else str(raw).lower() in ("true", "1", "yes")


debug_mode = config_bool("debug", False)
pd_lookup_enabled = config_bool("pd_lookup", True)
pd_lookup_backend = str(config.get("pd_backend", "ftp")).strip() or "ftp"
pd_comparator_limit = int(config.get("pd_comparator_limit", 10))
pd_sample_metadata_tsv = config.get("pd_sample_metadata_tsv")
pd_isolates_tsv = config.get("pd_isolates_tsv")
pd_exceptions_tsv = config.get("pd_exceptions_tsv")

# Helper: extract sample name from fastq filename
def extract_sample_name(filename):
    """
    Extracts the sample name from a FASTQ file, supporting common paired-end conventions.
    """
    base = os.path.basename(filename)
    if base.endswith('.gz'):
        base = base[:-3]  # Remove ".gz"
    # Match common paired-end naming conventions
    m = re.match(r"(.+?)(?:_R?[12](?:_001)?|_[12])(?:\.fastq)?$", base)
    if m:
        return m.group(1)  # Return sample name without paired-end suffix
    else:
        raise ValueError(f"Cannot extract sample name from filename: {filename}")


# Find R1 and R2 files
r1_files = glob.glob(os.path.join(input_dir, "*R1*.fastq*")) + glob.glob(os.path.join(input_dir, "*_1.fastq*"))
r2_files = glob.glob(os.path.join(input_dir, "*R2*.fastq*")) + glob.glob(os.path.join(input_dir, "*_2.fastq*"))

# Map sample names to paired files
sample_to_files = {}
for r1 in r1_files:
    try:
        sample = extract_sample_name(r1)
        r2_candidates = [r1.replace("_R1", "_R2"), r1.replace("_1", "_2")]  # Replace common patterns for R2
        found_r2 = None
        for r2 in r2_candidates:
            if r2 in r2_files:
                found_r2 = r2
                break
        if not found_r2:
            # Use a fallback check for unmatched R2 files
            found_r2 = next((r for r in r2_files if extract_sample_name(r) == sample), None)
        if found_r2:
            sample_to_files[sample] = {'r1': r1, 'r2': found_r2}
        else:
            print(f"Warning: No matching R2 file found for R1: {r1}")
    except ValueError as ve:
        print(f"Skipping invalid file: {r1} ({ve})")

samples = sorted(sample_to_files.keys())  # All valid sample names


# ----------------------------------------------
#   Validation: Check Sample Consistency
# ----------------------------------------------

def validate_sample_counts():
    """
    Ensures the total number of valid paired samples matches the input directory.
    """
    total_r1_files = len(r1_files)
    total_r2_files = len(r2_files)
    total_samples = len(samples)

    # Check if paired-end files align with the number of samples
    if total_r1_files != total_r2_files:
        print(
            f"Validation Warning: Mismatch in number of R1 ({total_r1_files}) and R2 ({total_r2_files}) files. "
            "Check input directory for unpaired FASTQ files."
        )
    
    # Ensure there are no unexpected duplicates
    if total_samples > total_r1_files or total_samples > total_r2_files:
        print(
            f"Validation Error: Detected more paired samples ({total_samples}) than available R1/R2 files. "
            "Check for duplicate or ambiguous sample names."
        )
        raise ValueError("Sample count validation failed: More paired samples than available R1/R2 files.")

    print(f"Total paired samples to process after validation = {total_samples}")

# Perform sample validation
validate_sample_counts()


# ----------------------------------------------
#   Define Wildcard Functions for Snakemake
# ----------------------------------------------

def get_r1(wildcards):
    """
    Retrieve the R1 file path for a given sample.
    """
    return sample_to_files[wildcards.sample]['r1']

def get_r2(wildcards):
    """
    Retrieve the R2 file path for a given sample.
    """
    return sample_to_files[wildcards.sample]['r2']


# -------------------------------------------------------------
#   Define function for inserting species names into Resfinder
# -------------------------------------------------------------

def load_resfinder_species_map(path):
    species_map = {}
    if os.path.exists(path):
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                species_map[row['sample']] = row['species']
    return species_map

resfinder_species_map = load_resfinder_species_map(os.path.join(output_dir, "data", "resfinder", "resfinder_species.tsv"))


# -------------------------------------------------------------
#   Define function for inserting organism names into AMRFinder
# -------------------------------------------------------------

def load_amrfinder_organism_map(path):
    organism_map = {}
    if os.path.exists(path):
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                organism_map[row['sample']] = row['amrfinderplus_organism']
    return organism_map

amrfinder_organism_map = load_amrfinder_organism_map(os.path.join(output_dir, "data", "amrfinderplus", "amrfinder_species.tsv"))
