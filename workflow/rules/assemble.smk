# ------------------------------------------
#        Clean input data (fastp)
# ------------------------------------------

rule fastp:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")
    params:
        html = "/dev/null/",
        json = "/dev/null/"
    resources:
        threads = 1,
        mem_mb = 10000,
        time = "1h"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.fastp.txt")
    conda: "../envs/assemble.yaml"
    shell:
        """
        # run fastp
        fastp -i {input.r1} \
              -I {input.r2} \
              -o {output.r1_clean} \
              -O {output.r2_clean} \
              --html {params.html} \
              --json {params.json}
        
        """

# ------------------------------------------
#         Assemble (Unicylcer)
# ------------------------------------------

rule unicylcer:
    input:
        r1_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")
    output:
        unicycler = protected(os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"))
    resources:
        threads = 32, 
        mem_mb = 150000,
        time = "6d"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.fastp_and_unicycler.txt")
    conda: "../envs/assemble.yaml"
    shell:
        """
    
        # run unicycler, keep only pertinent files
        unicycler -1 {input.r1_clean} \
                  -2 {input.r2_clean} \
                  -o $(dirname {output.unicycler}) \
                  -t {resources.threads} \
                  --keep 0
        
        """

# -----------------------------------------
#    Calculate genomic coverage (minimap2)
# -----------------------------------------

rule coverage:
    input:
        fasta = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        r1 = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")
    output:
        bam = temp(os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}.bam")),
        coverage = os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}_coverage.tsv")
    resources:
        threads = 1,
        mem_mb = 10000,
        time = "1h"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.coverage.txt")
    conda: "../envs/assemble.yaml"
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
    


