# ---------------------------------------------------------
#        Check completeness and contamination (CheckM)
# ---------------------------------------------------------

rule checkm:
    input:
        assembly = expand(os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"), sample=samples)
    output:
        checkm = protected(os.path.join(output_dir, "data", "checkm", "genome.stats.tsv")),
        checkm_stats = os.path.join(output_dir, "data", "checkm", "storage", "bin_stats.analyze.tsv")
    resources:
        mem_mb = 100000,
        time = "1-00:00:00",
        threads = 32
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "checkm.txt")
    conda: "../envs/checkm.yaml"
    params:
            checkm_db=config["checkm_db"]
    shell:
        """
        # CheckM
        export CHECKM_DATA_PATH={params.checkm_db}
        mkdir -p $(dirname {output.checkm})
        checkm lineage_wf $(dirname {input.assembly[0]}) $(dirname {output.checkm}) -x fasta -t {resources.threads} --tmpdir $(dirname {output.checkm})
        checkm qa $(dirname {output.checkm})/lineage.ms $(dirname {output.checkm}) -f {output.checkm}
        """


# ------------------------------------------------
#            Classify taxonomy (GTDB-tk)
# ------------------------------------------------

rule gtdbtk:
    input:
        assembly = expand(os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"), sample=samples)
    output:
        gtdbtk = protected(os.path.join(output_dir, "data", "gtdb-tk", "gtdbtk.bac120.summary.tsv"))
    resources:
        mem_mb = 150000,
        time = "1-00:00:00",
        threads = 32
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "gtdbtk.txt")
    conda: "../envs/gtdb.yaml"
    params:
            gtdbtk_db=config["gtdbtk_db"]
    shell:
        """
        # GTDB-tk
        export GTDBTK_DATA_PATH={params.gtdbtk_db}
        mkdir -p $(dirname {output.gtdbtk})
        gtdbtk classify_wf --genome_dir $(dirname {input.assembly[0]}) --out_dir $(dirname {output.gtdbtk}) -x fasta --cpus {resources.threads} --skip_ani_screen
        """


