# -----------------------------------------------------------------
#    Separate chromosome from plasmid(s), type plasmids (MOB-suite)
# -----------------------------------------------------------------

rule mobsuite:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        db_init = os.path.join(output_dir, "data", "mob-suite", ".mob_suite_db_initialized")
    output:
        mobtyper = os.path.join(output_dir, "data", "mob-suite", "{sample}", "mobtyper_results.txt"),
        blast    = os.path.join(output_dir, "data", "mob-suite", "{sample}", "biomarkers.blast.txt"),
        c_report = os.path.join(output_dir, "data", "mob-suite", "{sample}", "contig_report.txt")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mobsuite.txt")
    conda: "../envs/mob_suite.yaml"
    params:
        db_dir = "dbs/mob_suite"
    shell:
        """
        # run mob_recon to reconstruct and type plasmids leveraging unicycler circularity flags
        mob_recon --infile {input.assembly} \
        --outdir "$(dirname {output.mobtyper})" \
        --database_directory "{params.db_dir}" \
        --unicycler_contigs \
        -n {threads} \
        --force
        
        # Create dummy files if no plasmids were detected
        if [ ! -f {output.mobtyper} ]; then
            echo "# No plasmids detected - writing dummy file" > {output.mobtyper}
        fi
        
        if [ ! -f {output.blast} ]; then
            echo "# No plasmids detected - writing dummy file" > {output.blast}
        fi
        
        """

# -------------------------------------------------------------------------
#     Generate organism mapping for AMRFinderPlus (get_amrfinder_organism)
# -------------------------------------------------------------------------

rule get_amrfinder_organism:
    input:
        summary = os.path.join(output_dir, "data", "mash", "mash_taxonomy.tsv")
    output:
        organism_map = os.path.join(output_dir, "data", "amrfinderplus", "amrfinder_organism.tsv")
    conda: "../envs/amrfinder.yaml"
    shell:
        """
        
        python3 workflow/scripts/get_amrfinder_organism.py {input.summary} {output.organism_map}
        
        """
# -----------------------------------------------------
#     Detect AMR + Stress + Virulence  (AMRFinderPlus)
# -----------------------------------------------------

rule amrfinderplus:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        db_init = os.path.join(output_dir, "data", "amrfinderplus", ".afp_db_initialized"),
        organism_map = os.path.join(output_dir, "data", "amrfinderplus", "amrfinder_organism.tsv")
    output:
        afp = os.path.join(output_dir, "data", "amrfinderplus", "{sample}.afp.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.afp.txt")
    conda: "../envs/amrfinder.yaml"
    params:
        organism = lambda wildcards: amrfinder_organism_map.get(wildcards.sample, "Escherichia")
    shell:
        """
        
        # annotate AMR
        amrfinder -n {input.assembly} --plus \
        --name {wildcards.sample} \
        --threads {threads} \
        --organism {params.organism} \
        -o {output.afp} \
        -i 0.8 

       
        """

# ------------------------------------
#     Pull species for Resfinder
# ------------------------------------

rule get_resfinder_species:
    input:
        summary = os.path.join(output_dir, "data", "mash", "mash_taxonomy.tsv")
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
        assembly    = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        species_map = os.path.join(output_dir, "data", "resfinder", "resfinder_species.tsv"),
        res_db      = os.path.join(output_dir, "data", "resfinder", ".res_db_initialized")
    output:
        results_dir = os.path.join(output_dir, "data", "resfinder", "{sample}", "ResFinder_results_tab.txt"),
        pf_result   = os.path.join(output_dir, "data", "resfinder", "{sample}", "PointFinder_results.txt")
    params:
        species = lambda wildcards: resfinder_species_map.get(wildcards.sample, "Other")
    conda: "../envs/resfinder.yaml"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.resfinder.txt")
    shell:
        """
        
        db_res="dbs/resfinder_db"
        db_pt="dbs/pointfinder_db"
        db_dis="dbs/disinfinder_db"
        
        python -m resfinder \
        -ifa {input.assembly} \
        -o $(dirname {output.results_dir}) \
        -s "{params.species}" \
        -t 0.9 \
        -l 0.6 \
        --acquired \
        --point \
        --disinfectant \
        -db_res $db_res \
        -db_point $db_pt \
        -db_disinf $db_dis
        
        
        # Create dummy file if no point-mutations were detected
        if [ ! -f {output.pf_result} ]; then
            echo "# No point mutations detected - writing dummy file" > {output.pf_result}
        fi
        
        """

# ------------------------------------------------------
#     Comprehensively screen MGEs (MobileElementFinder)
# ------------------------------------------------------

rule mef:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta")
    output:
        mef        = os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.csv"),
        temp_fasta = temp(os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.tmp.fasta"))
    resources:
        mef_slots = 1
    conda: "../envs/mobileelementfinder.yaml"
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
