# ------------------------------------------
#   AMR-oriented annotation and helpers
# ------------------------------------------

# Rules in this module should infer organism context for AMR tools or annotate
# resistance, stress, virulence, and related AMR-focused features.

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

# -----------------------------------------------------
#     Detect AMR + Stress + Virulence  (AMRFinderPlus)
# -----------------------------------------------------

rule amrfinder:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        db_init = os.path.join(output_dir, "data", "amrfinderplus", ".afp_db_initialized"),
        organism_map = os.path.join(output_dir, "data", "amrfinderplus", "amrfinder_organism.tsv")
    output:
        afp = os.path.join(output_dir, "data", "amrfinderplus", "{sample}.afp.tsv")
    log:
        os.path.join(output_dir, "logs", "amrfinder", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.amrfinder.txt")
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
            echo "[swamg-rule] amrfinder"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] amrfinder"
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
