import pandas as pd
import sys

gtdb_tk_summary = sys.argv[1]
out_file = sys.argv[2]

df = pd.read_csv(gtdb_tk_summary, sep='\t')
salmonella_samples = []
for _, row in df.iterrows():
    sample = row['user_genome'].replace('.chromosome', '')
    classification = str(row['classification'])
    if "Salmonella" in classification:
        salmonella_samples.append(sample)
with open(out_file, 'w') as f:
    for sample in salmonella_samples:
        f.write(sample + '\n')