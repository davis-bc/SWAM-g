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

# -----------------------------------------------------------------
#   Separate chromosome from plasmid(s), type plasmids (MOB-suite)
# -----------------------------------------------------------------

rule mobsuite:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        db_init = os.path.join(output_dir, ".mob_suite_db_initialized")
    output:
        mob = os.path.join(output_dir, "data", "mob-suite", "{sample}", "chromosome.fasta"),
        chrom = os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"),
        plas = os.path.join(output_dir, "data", "assemblies", "plasmids", ".{sample}.plasmids.done"),
        mobtyper = os.path.join(output_dir, "data", "mob-suite", "{sample}", "mobtyper_results.txt")
    resources:
        mem_mb = 1000,
        time = "0-10:00:00",
        threads = 1
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mobsuite.txt")
    conda: "../envs/mob_suite.yaml"
    group: "group2"
    shell:
        """
        # run mob_recon to reconstruct and type plasmids leveraging unicycler circularity flags
        mob_recon --infile {input.assembly} --outdir $(dirname {output.mob}) --force --unicycler_contigs -n {resources.threads}

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
        assembly = os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"),
        plasmid = os.path.join(output_dir, "data", "assemblies", "plasmids", ".{sample}.plasmids.done"),
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
    group: "group3"
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
        species_map = os.path.join(output_dir, "data", "resfinder", "resfinder_species.tsv")
    shell:
        """
        python workflow/scripts/get_resfinder_species.py {input.summary} {output.species_map}

        """

# ------------------------------------------------------
#     Predict phenotype from genotype (Resfinder)
# ------------------------------------------------------

rule resfinder:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        species_map = os.path.join(output_dir, "data", "resfinder", "resfinder_species.tsv")
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

# --------------------------------
#      MobileElementFinder setup
# --------------------------------

rule mef_init:
    output:
        mef_init = os.path.join(output_dir, "bin", ".mef_initialized")
    shell:
        """
        mkdir -p $(dirname {output.mef_init}) && cd $(dirname {output.mef_init})
        
        # Clone MobileElementFinder if not already present
        if [ ! -d MobileElementFinder-1.1.2 ]; then
            wget https://files.pythonhosted.org/packages/42/09/7709dfe81fc6b695159c3cdc3d341b279af81e4e9986e545379642f07117/MobileElementFinder-1.1.2.tar.gz
            tar -xvf MobileElementFinder-1.1.2.tar.gz
            rm MobileElementFinder-1.1.2.tar.gz
            touch {output.mef_init}
        else
            echo "MobileElementFinder already exists, skipping"
        fi
        """

# ------------------------------------------------------
#     Comprehensively screen MGEs (MobileElementFinder)
# ------------------------------------------------------

rule mef:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        mef_init = os.path.join(output_dir, "bin", ".mef_initialized")
    output:
        mef = os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.csv")
    shell:
        """
        mkdir -p $(dirname {output.mef})
        
        python $(dirname {input.mef_init})/MobileElementFinder-1.1.2/mobileElementFinder.py find -c {input.assembly} $(dirname {output.mef})/{wildcards.sample} --temp-dir $(dirname {output.mef})
        
        """






