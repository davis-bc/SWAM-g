# -----------------------------------
#    Initialize MOB-suite database
# -----------------------------------

rule mob_init:
    output:
        touch(os.path.join(output_dir,".mob_suite_db_initialized"))
    conda:
        "../envs/mob_suite.yaml"
    shell:
        """
        MOB_DB=$(find .snakemake/conda -type d -path "*/mob_suite/databases" | head -n 1)
        if [ -z "$MOB_DB" ]; then
        echo "mob_suite database does not exist in conda environment, initializing..."
            mob_init
        fi
        touch {output}
        """

# --------------------------------------------------------------------------
#   Separate chromosome from plasmid(s), type plasmids and MGEs (MOB-suite)
# --------------------------------------------------------------------------

rule mobsuite:
    input:
        assembly = os.path.join(output_dir, "samples", "{sample}", "spades", "contigs.fasta"),
        db_init = os.path.join(output_dir,".mob_suite_db_initialized")
    output:
        mob = os.path.join(output_dir, "samples", "{sample}", "mob-suite", "chromosome.fasta"),
        chrom = os.path.join(output_dir, "assemblies", "chromosomes", "{sample}.chromosome.fasta"),
        plas = os.path.join(output_dir, "assemblies", "plasmids", ".{sample}.plasmids.done"),
        mobtyper = os.path.join(output_dir, "samples", "{sample}", "mob-suite", "mobtyper_results.txt")
    resources:
        mem_mb = 1000,
        time = "0-10:00:00",
        threads = 1
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mobsuite.txt")
    conda: "../envs/mob_suite.yaml"
    group: "group1"
    shell:
        """
        # run mob_recon to reconstruct and type plasmids
        mob_recon --infile {input.assembly} --outdir $(dirname {output.mob}) --force -n {resources.threads}

        # copy chromosome assembly to new directory
        mkdir -p $(dirname {output.chrom})
        cp {output.mob} {output.chrom}

        # copy plasmid(s) to new directory
        for f in $(dirname {output.mob})/plasmid*; do
            if [ -f $f ]; then
                base=$(basename $f)
                cp $f $(dirname {output.plas})/{wildcards.sample}.$base
            fi
        done

        # Always create mobtyper_results.txt in the absence of plasmid contigs
        if [ ! -f {output.mobtyper} ]; then
        touch {output.mobtyper}
        fi

        # create marker file to satisfy output
        touch {output.plas}
        """
