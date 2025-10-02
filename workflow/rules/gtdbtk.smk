# ------------------------------------------------
#   Taxonomically classify chromosomes (GTDB-tk)
# ------------------------------------------------

rule gtdbtk:
    input:
        assembly = expand(os.path.join(output_dir, "assemblies", "chromosomes", "{sample}.chromosome.fasta"), sample=samples)
    output:
        gtdbtk = os.path.join(output_dir, "data", "gtdb-tk", "gtdbtk.bac120.summary.tsv")
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
