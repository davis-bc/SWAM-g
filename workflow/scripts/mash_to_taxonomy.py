#!/usr/bin/env python3
"""
mash_to_taxonomy.py  –  Aggregate per-sample MASH screen results into a
taxonomy summary TSV compatible with the existing GTDB-tk parsers.

Usage:
    python mash_to_taxonomy.py <output_tsv> <sample1.mash_screen.tsv> [sample2 ...]

Output columns:
    user_genome              – {sample}.chromosome
    classification           – g__Genus;s__Genus species   (GTDB-style minimal)
    closest_genome_reference – RefSeq accession of the top hit
"""

import sys
import os
import re


def parse_organism(comment: str):
    """
    Extract (genus, species_binomial) from a MASH screen query-comment field.

    Comments look like:
      "Escherichia coli str. K-12 substr. MG1655, complete genome."
      "[Salmonella enterica subsp. enterica serovar Typhimurium str. LT2], complete genome."
    Returns ("Escherichia", "Escherichia coli") or (None, None).
    """
    # Strip leading/trailing brackets, dots, brackets
    name = re.sub(r"^\[|\].*$", "", comment).strip()
    # Take first two whitespace-separated tokens as genus + specific epithet
    tokens = name.split()
    if len(tokens) >= 2:
        genus = tokens[0]
        epithet = tokens[1]
        return genus, f"{genus} {epithet}"
    if len(tokens) == 1:
        return tokens[0], tokens[0]
    return None, None


def best_hit(mash_tsv: str):
    """
    Return (accession, genus, species_binomial) for the top hit in a
    sorted (descending identity) MASH screen result file.
    Returns ("unknown", None, None) if the file is empty or unreadable.
    """
    try:
        with open(mash_tsv) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                accession = parts[4].split()[0]          # strip extra tokens after whitespace
                comment = "\t".join(parts[5:]) if len(parts) > 5 else ""
                genus, species = parse_organism(comment)
                return accession, genus, species
    except OSError:
        pass
    return "unknown", None, None


def sample_name_from_path(path: str) -> str:
    """Derive sample name from e.g. 'SRR123.mash_screen.tsv' → 'SRR123'."""
    base = os.path.basename(path)
    # Strip known suffixes
    for suffix in (".mash_screen.tsv", ".tsv"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return base


def main(output_tsv: str, mash_files: list):
    with open(output_tsv, "w") as out:
        out.write("user_genome\tclassification\tclosest_genome_reference\n")
        for mf in mash_files:
            sample = sample_name_from_path(mf)
            accession, genus, species = best_hit(mf)
            if genus and species:
                classification = f"g__{genus};s__{species}"
            elif genus:
                classification = f"g__{genus};s__{genus}"
            else:
                classification = "g__unknown;s__unknown"
            out.write(f"{sample}.chromosome\t{classification}\t{accession}\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.stderr.write(
            "Usage: mash_to_taxonomy.py <output_tsv> <sample1.mash_screen.tsv> [sample2 ...]\n"
        )
        sys.exit(1)
    main(sys.argv[1], sys.argv[2:])
