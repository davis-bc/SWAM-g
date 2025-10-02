# ---------------------------------------------------------
#        Check completeness and contamination (CheckM)
# ---------------------------------------------------------

rule checkm:
    input:
        assembly = expand(os.path.join(output_dir, "assemblies", "chromosomes", "{sample}.chromosome.fasta"), sample=samples)
    output:
        checkm = os.path.join(output_dir, "data", "checkm", "genome.stats.tsv")
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
        checkm lineage_wf $(dirname {input.assembly[0]}) $(dirname {output.checkm}) -x fasta -t {resources.threads}
        checkm qa $(dirname {output.checkm})/lineage.ms $(dirname {output.checkm}) -f {output.checkm}
        """
