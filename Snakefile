import os
import glob

input_dir = config["in_dir"]
output_dir = config["out_dir"]

### Find all sample names by looking for *_1.fastq files
samples = [os.path.basename(f).replace("_1.fastq", "") for f in glob.glob(os.path.join(input_dir, "*_1.fastq"))]

print(samples)

### Define final targets for pipeline
rule all:
    input:
        expand(os.path.join(output_dir, "samples", "{sample}", "chromosome", "scaffolds.fasta"), sample=samples),
        expand(os.path.join(output_dir, "samples", "{sample}", "plasmid", "scaffolds.fasta"), sample=samples)

# ---------------------------
# 1. Assemble (fastp + SPAdes)
# ---------------------------
rule fastp_and_spades:
    input:
        r1 = os.path.join(input_dir, "{sample}_1.fastq"),
        r2 = os.path.join(input_dir, "{sample}_2.fastq")
    output:
        r1_clean = os.path.join(output_dir, "samples", "{sample}", "{sample}_R1.clean.fastq.gz"),
        r2_clean = os.path.join(output_dir, "samples", "{sample}", "{sample}_R2.clean.fastq.gz"),
        chrom_scaffolds = os.path.join(output_dir, "samples", "{sample}", "chromosome", "scaffolds.fasta"),
        plasmid_scaffolds = os.path.join(output_dir, "samples", "{sample}", "plasmid", "scaffolds.fasta"),
	chromosomes = directory(os.path.join(output_dir, "chromosomes")
    params:
        html = "/dev/null/",
        json = "/dev/null/"
    threads: 8
    resources:
        mem_mb=20000
    conda: "envs/assemble.yaml"
    shell:
	"""
        # fastp
        mkdir -p $(dirname {output.r1_clean})
        fastp -i {input.r1} -I {input.r2} -o {output.r1_clean} -O {output.r2_clean} --html {params.html} --json {params.json}

        # chromosome assembly
        mkdir -p $(dirname {output.chrom_scaffolds})
        spades.py --isolate -1 {output.r1_clean} -2 {output.r2_clean} -o $(dirname {output.chrom_scaffolds})
	
        # plasmid assembly
        mkdir -p $(dirname {output.plasmid_scaffolds})
        spades.py --plasmid -1 {output.r1_clean} -2 {output.r2_clean} -o $(dirname {output.plasmid_scaffolds})
        
	# copy chromosome scaffolds to new directory
	mkdir -p $(dirname {output.chromosomes})
	cp {output.chrom_scaffolds} {output.chromosomes}/{wildcards.sample}.scaffold.fasta
	"""

group: "assembly"


# ------------------------------------------------------------
# 2. Classify and Annotate (GTDB-tk + AMRFinderPlus + CheckM)
# ------------------------------------------------------------

rule classify_annotate:
    input:
        chromosomes = os.path.join(output_dir, "chromosomes", "{sample}.scaffold.fasta")
    output:
        gtdbtk = os.path.join(output_dir, "GTDB-tk", "gtdbtk.bac120.summary.tsv"),
        amrf = directory(os.path.join(output_dir, "AMRFinderPlus")),
        checkm = os.path.join(output_dir, "CheckM", "genome.stats.tsv")
    threads: 32
    resources:
        mem_mb = 100000,
        time = "1-00:00:00"
    conda: "envs/classify.yaml"
    shell:
        """
        mkdir -p $(dirname {output.gtdbtk})
        gtdbtk classify_wf --genome_dir $(dirname {input.chromosomes}) \
            --out_dir $(dirname {output.gtdbtk}) -x fasta --cpus {threads} --skip_ani_screen

        mkdir -p {output.amrf}
        amrfinder -n {input.chromosomes} --plus --name {wildcards.sample} \
            --threads {threads} \
            -o {output.amrf}/{wildcards.sample}.amrfinderplus.txt

        mkdir -p $(dirname {output.checkm})
        checkm lineage_wf $(dirname {input.chromosomes}) $(dirname {output.checkm}) -x fasta --tmpdir "$TMPDIR" -t {threads}
        checkm qa $(dirname {output.checkm})/lineage.ms $(dirname {output.checkm}) -f {output.checkm}
        """

group: "classify"
