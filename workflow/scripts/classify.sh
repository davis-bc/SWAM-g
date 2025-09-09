source ativate wgs-qc

### Set env variables
home="/work/NRSAAMR/Projects/SWAM/WGS"
input="/work/NRSAAMR/Projects/SWAM/WGS/input"
out="/work/NRSAAMR/Projects/SWAM/WGS/output"

export TMPDIR="/work/NRSAAMR/Projects/SWAM/WGS/work"

### Assign taxonomy
#gtdbtk classify_wf --genome_dir "$out"/chromosome --out_dir "$mags"/gtdb.out  -x fasta --cpus 32 --skip_ani_screen

### Check completeness and contamination
#checkm lineage_wf "$mags" "$mags"/checkm -x fa --tmpdir "$tmp" -t 32
#checkm qa "$mags"/checkm/lineage.ms "$mags"/checkm -f "$mags"/checkm/binstats.tsv

### Run through AMRFinderPlus
#amrfinder -n "$mag" --plus --name $base --threads 32 -o "$home"/amrfinderplus/$base --organism Salmonella

### Run through MobileElementsFinder
#mefinder find --contig /path/to/genome.fna output_name

### If Salmonella, run chromosome through SeqSero2
#SeqSero2_package.py -m k -t 4 -i assembly.fasta
