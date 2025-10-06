# ------------------------------------------
#   Clean and Assemble (fastp + Unicylcer)
# ------------------------------------------

rule fastp_and_unicylcer:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1_clean = protected(os.path.join(output_dir, "samples", "{sample}", "{sample}_R1.clean.fastq.gz")),
        r2_clean = protected(os.path.join(output_dir, "samples", "{sample}", "{sample}_R2.clean.fastq.gz")),
        unicycler = protected(os.path.join(output_dir, "samples", "{sample}", "unicycler", "assembly.fasta"))
    params:
        html = "/dev/null/",
        json = "/dev/null/"
    resources:
        mem_mb = 20000,
        threads = 16,
        time = "0-10:00:00"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.fastp_and_unicycler.txt")
    conda: "../envs/assemble.yaml"
    shell:
        """
        # run fastp
        mkdir -p $(dirname {output.r1_clean})
        fastp -i {input.r1} -I {input.r2} -o {output.r1_clean} -O {output.r2_clean} --html {params.html} --json {params.json}

        # run unicycler, keep only pertinent files
        mkdir -p $(dirname {output.unicycler})
        unicycler -1 {output.r1_clean} -2 {output.r2_clean} -o $(dirname {output.unicycler}) -t {resources.threads} --keep 0
        """

# ---------------------------------
#    Calculate genomic coverage
# ---------------------------------

rule map_and_calculate_coverage:
    input:
        fasta = os.path.join(output_dir, "samples", "{sample}", "unicycler", "assembly.fasta"),
        r1 = os.path.join(output_dir, "samples", "{sample}", "{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(output_dir, "samples", "{sample}", "{sample}_R2.clean.fastq.gz")
    output:
        bam = os.path.join(output_dir, "data", "bams", "{sample}.bam"),
        coverage = os.path.join(output_dir, "data", "bams", "{sample}_coverage.tsv")
    resources:
        mem_mb = 10000,
        threads = 1,
        time = "0-01:00:00"
    group: "group1"
    shell:
        r"""
        # Map reads using minimap2
        minimap2 -ax sr {input.fasta} {input.r1} {input.r2} | \
        samtools view -@ {resources.threads} -bS - | \
        samtools sort -@ {resources.threads} -o {output.bam}

        # Index BAM for coverage calculation
        samtools index {output.bam}

        # Calculate genome coverage with samtools depth
        cov=$(samtools depth {output.bam} | awk '{{sum+=$3; cnt++}} END {{print (cnt > 0 ? sum/cnt : 0)}}')

        # Write result as TSV: sample | coverage
        echo -e "{wildcards.sample}\t$cov" > {output.coverage}
        """





