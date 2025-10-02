# ------------------------------------------
#   Clean and Assemble (fastp + Spades)
# ------------------------------------------

rule fastp_and_spades:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1_clean = protected(os.path.join(output_dir, "samples", "{sample}", "{sample}_R1.clean.fastq.gz")),
        r2_clean = protected(os.path.join(output_dir, "samples", "{sample}", "{sample}_R2.clean.fastq.gz")),
        spades = protected(os.path.join(output_dir, "samples", "{sample}", "spades", "contigs.fasta"))
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
        mkdir -p $(dirname {output.spades})
        spades.py --isolate -1 {output.r1_clean} -2 {output.r2_clean} -o $(dirname {output.spades}) -t {resources.threads}
        """
