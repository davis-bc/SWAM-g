# --------------------------------------
#   Initialize AMRFinderPlus database
# --------------------------------------

rule afp_init:
    output:
        touch(os.path.join(output_dir,".afp_db_initialized"))
    conda:
        "../envs/amrfinder.yaml"
    shell:
        """
        amrfinder -u
        touch {output}
        """

# -------------------------------------------------------------------------------
#   Detect AMR + Stress + Virulence in chromosomes and plasmids (AMRFinderPlus)
# -------------------------------------------------------------------------------

rule annotate_amr:
    input:
        assembly = os.path.join(output_dir, "assemblies", "chromosomes", "{sample}.chromosome.fasta"),
        plasmid = os.path.join(output_dir, "assemblies", "plasmids", ".{sample}.plasmids.done"),
        db_init = os.path.join(output_dir,".afp_db_initialized")
    output:
        afp = os.path.join(output_dir, "data", "amrfinderplus", "chromosomes", "{sample}.afp.tsv"),
        plas = os.path.join(output_dir, "data", "amrfinderplus", "plasmids", ".{sample}.plasmids.afp.done")
    resources:
        mem_mb = 1000,
        time = "0-10:00:00",
        threads = 1
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.afp.txt")
    conda: "../envs/amrfinder.yaml"
    group: "group2"
    shell:
        """
        # annotate chromosomes
        mkdir -p $(dirname {output.afp})
        amrfinder -n {input.assembly} --plus --name {wildcards.sample} --threads {resources.threads} -o {output.afp} -i 0.8

        # annotate plasmids
        mkdir -p $(dirname {output.plas})
        for f in $(dirname {input.plasmid})/{wildcards.sample}*plasmid*; do
            if [ -f $f ]; then
                base=$(basename $f .fasta)
                amrfinder -n $f --plus --name $base --threads {resources.threads} -o $(dirname {output.plas})/$base.afp.tsv -i 0.8
            fi
        done

        # create marker file to satisfy output in the absence of plasmids
        touch {output.plas}
        """

# ------------------------------------
#     Pull species for Resfinder
# ------------------------------------

rule get_resfinder_species:
    input:
        summary = os.path.join(output_dir, "data", "gtdb-tk", "gtdbtk.bac120.summary.tsv")
    output:
        species_map = os.path.join(output_dir, "data", "resfinder_species.tsv")
    shell:
        """
        python workflow/scripts/get_resfinder_species.py {input.summary} {output.species_map}

        """

# ------------------------------------------------------
#     Predict phenotype from genotype (Resfinder)
# ------------------------------------------------------

rule resfinder:
    input:
        assembly = os.path.join(output_dir, "samples", "{sample}", "spades", "contigs.fasta"),
        species_map = os.path.join(output_dir, "data", "resfinder_species.tsv")
    output:
        results_dir = os.path.join(output_dir, "data", "resfinder", "{sample}", "ResFinder_results_tab.txt")
    params:
        species = lambda wildcards: resfinder_species_map.get(wildcards.sample, "Other"),
        res_db = config["res_db"],
        pt_db = config["pt_db"],
        dis_db = config["dis_db"]
    conda: "../envs/resfinder.yaml"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.resfinder.txt")
    shell:
        """
        python -m resfinder -ifa {input.assembly} -o $(dirname {output.results_dir}) -s "{params.species}" -t 0.9 -l 0.6 --acquired --point --disinfectant -db_res {params.res_db} -db_point {params.pt_db} -db_disinf {params.dis_db}
        """
