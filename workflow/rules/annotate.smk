# -----------------------------------------------------------------
#   Separate chromosome from plasmid(s), type plasmids (MOB-suite)
# -----------------------------------------------------------------

rule mobsuite:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        db_init = os.path.join(output_dir, "data", "mob-suite", ".mob_suite_db_initialized")
    output:
        mob      = os.path.join(output_dir, "data", "mob-suite", "{sample}", "chromosome.fasta"),
        chrom    = os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"),
        plas     = os.path.join(output_dir, "data", "assemblies", "plasmids", ".{sample}.plasmids.done"),
        mobtyper = os.path.join(output_dir, "data", "mob-suite", "{sample}", "mobtyper_results.txt"),
        blast    = os.path.join(output_dir, "data", "mob-suite", "{sample}", "biomarkers.blast.txt"),
        c_report = os.path.join(output_dir, "data", "mob-suite", "{sample}", "contig_report.txt")
    resources:
        mem_mb = 10000,
        time = "0-10:00:00",
        threads = 32
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mobsuite.txt")
    conda: "../envs/mob_suite.yaml"
    group: "mobsuite"
    shell:
        """
        # run mob_recon to reconstruct and type plasmids leveraging unicycler circularity flags
        mob_recon --infile {input.assembly} --outdir $(dirname {output.mob}) --force --unicycler_contigs -n {resources.threads}

        # copy chromosome assembly to new directory
        cp {output.mob} {output.chrom}

        # copy plasmid(s) to new directory
        for f in $(dirname {output.mob})/plasmid*; do
            if [ -f $f ]; then
                base=$(basename $f)
                cp $f $(dirname {output.plas})/{wildcards.sample}.$base
            fi
        done

        # Create empty files if no plasmid files or other outputs are generated
        if [ ! -f {output.mobtyper} ]; then
            echo "Creating empty mobtyper_results.txt."
            touch {output.mobtyper}
        fi

        if [ ! -f {output.blast} ]; then
            echo "Creating empty biomarkers.blast.txt."
            touch {output.blast}
        fi

        if [ ! -f {output.c_report} ]; then
            echo "Creating empty contig_report.txt."
            touch {output.c_report}
        fi
        
        touch {output.plas}
        """

# ----------------------------------------------------------------------
#     Generate organism mapping for AMRFinderPlus (get_amrfinder_organism)
# ----------------------------------------------------------------------

rule get_amrfinder_organism:
    input:
        summary = os.path.join(output_dir, "data", "gtdb-tk", "gtdbtk.bac120.summary.tsv")
    output:
        organism_map = os.path.join(output_dir, "data", "amrfinderplus", "amrfinder_organism.tsv")
    conda: "../envs/amrfinder.yaml"
    shell:
        """
        
        python workflow/scripts/get_amrfinder_organism.py {input.summary} {output.organism_map}
        
        """
# -------------------------------------------------------------------------------
#   Detect AMR + Stress + Virulence in chromosomes and plasmids (AMRFinderPlus)
# -------------------------------------------------------------------------------

rule amrfinderplus:
    input:
        assembly = os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"),
        plasmid = os.path.join(output_dir, "data", "assemblies", "plasmids", ".{sample}.plasmids.done"),
        db_init = os.path.join(output_dir, "data", "amrfinderplus", ".afp_db_initialized"),
        organism_map = os.path.join(output_dir, "data", "amrfinderplus", "amrfinder_organism.tsv")
    output:
        afp = os.path.join(output_dir, "data", "amrfinderplus", "chromosomes", "{sample}.afp.tsv"),
        plas = os.path.join(output_dir, "data", "amrfinderplus", "plasmids", ".{sample}.plasmids.afp.done")
    resources:
        mem_mb = 10000,
        time = "0-10:00:00",
        threads = 1
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.afp.txt")
    conda: "../envs/amrfinder.yaml"
    group: "amrfinderplus"
    params:
        organism = lambda wildcards: amrfinder_organism_map.get(wildcards.sample, "Escherichia")
    shell:
        """
        # annotate chromosomes
        amrfinder -n {input.assembly} --plus --name {wildcards.sample} --threads {resources.threads} --organism {params.organism} -o {output.afp} -i 0.8 > /dev/null 2>&1

        # annotate plasmids
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
        assembly    = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        species_map = os.path.join(output_dir, "data", "resfinder", "resfinder_species.tsv"),
        res_db      = os.path.join(output_dir, "data", "resfinder", ".res_db_initialized")
    output:
        results_dir = os.path.join(output_dir, "data", "resfinder", "{sample}", "ResFinder_results_tab.txt"),
        pf_result   = os.path.join(output_dir, "data", "resfinder", "{sample}", "PointFinder_results.txt")
    params:
        species = lambda wildcards: resfinder_species_map.get(wildcards.sample, "Other")
    conda: "../envs/resfinder.yaml"
    group: "resfinder"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.resfinder.txt")
    shell:
        """
        
        db_res="dbs/resfinder_db"
        db_pt="dbs/pointfinder_db"
        db_dis="dbs/disinfinder_db"
        
        python -m resfinder -ifa {input.assembly} -o $(dirname {output.results_dir}) -s "{params.species}" -t 0.9 -l 0.6 --acquired --point --disinfectant -db_res $db_res -db_point $db_pt -db_disinf $db_dis
        
        """

# ------------------------------------------------------
#     Comprehensively screen MGEs (MobileElementFinder)
# ------------------------------------------------------

rule mef:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta")
    output:
        mef        = os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.csv"),
        temp_fasta = temp(os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.tmp.fasta"))
    resources:
        mem_mb = 1000,
        time = "0-10:00:00",
        threads = 1
    conda: "../envs/mobileelementfinder.yaml"
    group: "mef"
    params:
        temp_dir = os.path.join(output_dir, "data", "mobileelementfinder", "{sample}")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mef.txt")
    shell:
        """
        
        # Preprocess FASTA headers into a temporary file
        sed -E 's/^(>[^ ]+).*/\\1/' {input.assembly} > {output.temp_fasta}

        # Run mefinder on the temporary FASTA file
        mefinder find -c {output.temp_fasta} $(dirname {output.mef})/{wildcards.sample} --temp-dir {params.temp_dir}
        
        """






