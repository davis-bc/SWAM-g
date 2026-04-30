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
    log:
        os.path.join(output_dir, "logs", "mobsuite", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mobsuite.txt")
    conda: "../envs/mob_suite.yaml"
    params:
        db_dir = "dbs/mob_suite"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] mobsuite"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] mobsuite"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

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

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
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
    log:
        os.path.join(output_dir, "logs", "amrfinderplus", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.afp.txt")
    conda: "../envs/amrfinder.yaml"
    params:
        organism = lambda wildcards: amrfinder_organism_map.get(wildcards.sample, "Escherichia")
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] amrfinderplus"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] amrfinderplus"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        # annotate AMR
        amrfinder -n {input.assembly} --plus \
        --name {wildcards.sample} \
        --threads {threads} \
        --organism {params.organism} \
        -o {output.afp} \
        -i 0.8 

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
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
    log:
        os.path.join(output_dir, "logs", "resfinder", "{sample}.log")
    params:
        species = lambda wildcards: resfinder_species_map.get(wildcards.sample, "Other")
    conda: "../envs/resfinder.yaml"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.resfinder.txt")
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] resfinder"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] resfinder"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

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

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
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
    log:
        os.path.join(output_dir, "logs", "mef", "{sample}.log")
    resources:
        mef_slots = 1
    conda: "../envs/mobileelementfinder.yaml"
    params:
        temp_dir = os.path.join(output_dir, "data", "mobileelementfinder", "{sample}")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mef.txt")
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] mef"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] mef"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        # Preprocess FASTA headers into a temporary file
        sed -E 's/^(>[^ ]+).*/\\1/' {input.assembly} > {output.temp_fasta}

        # Skip mefinder if the assembly contains no sequences (empty FASTA causes
        # BLAST to abort with a "CFastaReader: Near line 1" error).
        # Write a stub CSV that matches mefinder's 5-comment-line + header format
        # so downstream parsers (data_summary.R, skip=5) see a zero-row data frame.
        if [ "$(grep -c '^>' {output.temp_fasta} || true)" -eq 0 ]; then
            printf '# No sequences in assembly - mefinder skipped\n#\n#\n#\n#\nname,type,contig,start,end\n' > {output.mef}
        else
            # Run mefinder on the temporary FASTA file
            mefinder find -c {output.temp_fasta} $(dirname {output.mef})/{wildcards.sample} --temp-dir {params.temp_dir}
        fi

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """
