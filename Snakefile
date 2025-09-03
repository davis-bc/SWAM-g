import glob
import re
import os

configfile: "config.yaml"

# Auto-detect sample names from input directory
fastqs = glob.glob("input/*_1.fastq")
SAMPLES = [re.sub(r"_1.fastq$", "", os.path.basename(f)) for f in fastqs]

rule all:
    input:
        expand("output/{sample}/plasmids/mob_recon.done", sample=SAMPLES)

# ---------------------------
# 1. Assembly (fastp + SPAdes)
# ---------------------------
rule assemble:
    input:
        r1 = "input/{sample}_1.fastq",
        r2 = "input/{sample}_2.fastq"
    output:
        clean_r1 = "output/{sample}/clean/{sample}_R1.clean.fastq.gz",
        clean_r2 = "output/{sample}/clean/{sample}_R2.clean.fastq.gz",
        chromo_dir = directory("output/{sample}/chromosome"),
        plasmid_dir = directory("output/{sample}/plasmid")
    threads: 8
    resources:
        mem_mb = 20000,
        time = "1-00:00:00"
    conda: "envs/assemble.yaml"
    shell:
        """
        mkdir -p output/{wildcards.sample}/clean
        fastp -i {input.r1} -I {input.r2} \
              -o {output.clean_r1} -O {output.clean_r2} \
              --html /dev/null --json /dev/null

        spades.py --isolate -1 {input.r1} -2 {input.r2} -o {output.chromo_dir}
        spades.py --plasmid -1 {input.r1} -2 {input.r2} -o {output.plasmid_dir}
        """

# ---------------------------
# 2. Classification (taxonomy, QC, ARGs)
# ---------------------------
rule classify:
    input:
        chromo = "output/{sample}/chromosome/contigs.fasta"
    output:
        gtdbtk = "output/{sample}/classify/gtdbtk.done",
        amrf = "output/{sample}/classify/amrfinderplus.done"
    threads: 32
    resources:
        mem_mb = 64000,
        time = "2-00:00:00"
    conda: "envs/classify.yaml"
    shell:
        """
        mkdir -p output/{wildcards.sample}/classify

        gtdbtk classify_wf \
            --genome_dir output/{wildcards.sample}/chromosome \
            --out_dir output/{wildcards.sample}/classify/gtdbtk \
            -x fasta --cpus {threads} --skip_ani_screen
        touch {output.gtdbtk}

        amrfinder -n {input.chromo} --plus --name {wildcards.sample} \
                  --threads {threads} \
                  -o output/{wildcards.sample}/classify/amrfinderplus.txt
        touch {output.amrf}
        """

# ---------------------------
# 3. Plasmid reconstruction
# ---------------------------
rule plasmids:
    input:
        asm = "output/{sample}/plasmid/contigs.fasta"
    output:
        "output/{sample}/plasmids/mob_recon.done"
    threads: 8
    resources:
        mem_mb = 16000,
        time = "12:00:00"
    conda: "envs/plasmids.yaml"
    shell:
        """
        mkdir -p output/{wildcards.sample}/plasmids
        mob_recon --infile {input.asm} --outdir output/{wildcards.sample}/plasmids
        touch {output}
        """
