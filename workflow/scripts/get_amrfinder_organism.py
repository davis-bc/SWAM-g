#!/usr/bin/env python3

import sys
import csv
import re
import subprocess


def get_amrfinderplus_organisms():
    """
    Retrieve available --organism options from the output of `amrfinder -l`.
    Returns a dictionary with valid genus and genus_species keys.
    """
    try:
        process = subprocess.Popen(
            ["amrfinder", "-l"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            sys.stderr.write(f"Error running `amrfinder -l`: {stderr}\n")
            sys.exit(process.returncode)

        organism_line = next(
            (line for line in stdout.splitlines() if line.startswith("Available --organism options")),
            "",
        )

        # Extract the organism options, clean them up, and store in a dictionary
        organisms = organism_line.split(":", 1)[1].strip().split(", ")
        valid_flags = {flag.strip(): True for flag in organisms}  # Normalize any flag spacing
        return valid_flags
    except StopIteration:
        sys.stderr.write("Failed to parse organism options from `amrfinder -l` output.\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"An unexpected error occurred: {str(e)}\n")
        sys.exit(1)


def extract_genus_species(taxonomy):
    """
    Extract genus and species from the GTDB taxonomy string.
    Returns genus and species as separate strings.
    """
    match = re.match(r".*;g__([A-Za-z0-9_]+);s__([A-Za-z0-9_\. ]+)", taxonomy)
    if match:
        genus = match.group(1).strip()
        species = match.group(2).strip()
        return genus, species
    return None, None


def sanitize_species(genus, species):
    """
    Sanitize the species string to ensure it matches the AMRFinderPlus format.
    If the species includes the genus (e.g., 'Enterococcus faecalis'), exclude the genus.
    """
    sanitized = species.replace(" ", "_").replace(".", "_")
    if sanitized.startswith(genus + "_"):  # Remove redundant genus from species
        sanitized = sanitized[len(genus) + 1 :]
    return sanitized


def get_organism_flag(genus, species, valid_flags, default="Escherichia"):
    """
    Determine the AMRFinderPlus organism flag given genus, species, and valid flags.
    Prioritizes full genus_species match before falling back to genus.
    """
    if genus and species:
        # Sanitize species and construct genus_species flag
        sanitized_species = sanitize_species(genus, species)
        genus_species_flag = f"{genus}_{sanitized_species}"

        if genus_species_flag in valid_flags:
            return genus_species_flag
        elif genus in valid_flags:
            return genus
    elif genus and genus in valid_flags:
        return genus
    return default


def main(gtdb_summary, output_map):
    """
    Map GTDB taxonomy to AMRFinderPlus organism flags and write the mapping to a TSV file.
    """
    # Retrieve valid --organism flags and store them in a dictionary
    valid_flags = get_amrfinderplus_organisms()

    with open(gtdb_summary, newline="") as infile, open(output_map, "w", newline="") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["sample", "amrfinderplus_organism"])

        for row in reader:
            sample = row.get("user_genome", "").strip()
            classification = row.get("classification", "").strip()

            if not sample or not classification:
                writer.writerow([sample, "Escherichia"])  # Default organism for missing or invalid data
                continue

            # Parse genus and species from classification
            genus, species = extract_genus_species(classification)

            if not genus:
                writer.writerow([sample, "Escherichia"])  # Default organism for missing genus
                continue

            # Get organism flag using genus and species
            organism_flag = get_organism_flag(genus, species, valid_flags)

            # Write the result to the output file
            writer.writerow([sample, organism_flag])


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: get_amrfinderplus_species.py <gtdb_summary> <output_map>\n")
        sys.exit(1)

    gtdb_summary = sys.argv[1]
    output_map = sys.argv[2]
    main(gtdb_summary, output_map)
