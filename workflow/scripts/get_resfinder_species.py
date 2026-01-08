import pandas as pd
import sys

gtdbtk_summary = sys.argv[1]
out_file = sys.argv[2]

# List of supported species/models
RESFINDER_SPECIES = [
    "Campylobacter", "Campylobacter jejuni", "Campylobacter coli", 
    "Enterococcus faecalis", "Enterococcus faecium", "Escherichia coli",
    "Helicobacter pylori", "Klebsiella", "Mycobacterium tuberculosis",
    "Neisseria gonorrhoeae", "Plasmodium falciparum", "Salmonella", 
    "Salmonella enterica", "Staphylococcus aureus"
]

def assign_species(classification):
    species = None
    genus = None
    for part in str(classification).split(';'):
        part = part.strip()
        if part.startswith("s__"):
            species = part[3:]
        if part.startswith("g__"):
            genus = part[3:]
    # Try full species match
    if species and species in RESFINDER_SPECIES:
        return species
    # Try genus-level match
    if genus:
        if genus in RESFINDER_SPECIES:
            return genus
        if species:
            combo = f"{genus} {species}"
            for s in RESFINDER_SPECIES:
                if s.lower() == combo.lower():
                    return s
    if species:
        for s in RESFINDER_SPECIES:
            if s.lower() in species.lower():
                return s
    return "Other"

df = pd.read_csv(gtdbtk_summary, sep='\t')
with open(out_file, 'w') as out:
    out.write("sample\tspecies\n")
    for _, row in df.iterrows():
        sample = row['user_genome'].replace('.chromosome', '')
        species = assign_species(row['classification'])
        out.write(f"{sample}\t{species}\n")