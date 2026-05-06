# ------------------------------------------
#   Taxonomy and general sample typing
# ------------------------------------------

# Rules in this module should classify the isolate at the taxonomy or typing
# level without focusing on downstream AMR or mobility annotation.

# -------------------------------------------------------
#   Screen assembly against EnteroRef MASH sketch (MASH)
# -------------------------------------------------------

rule mash:
    input:
        assembly  = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        msh_ready = enteroref_sketch_ready
    output:
        mash_result = os.path.join(output_dir, "data", "mash", "{sample}.mash_screen.tsv")
    log:
        os.path.join(output_dir, "logs", "mash", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mash.txt")
    conda: "../envs/mash.yaml"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] mash"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] mash"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        mkdir -p $(dirname {output.mash_result})
        mash screen -p {threads} -w dbs/EnteroRef_GTDBSketch_20231003_V2.msh {input.assembly} \
            | sort -grk1,1 > {output.mash_result}

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """

# -------------------------------------------------------
#   Combine per-sample MASH results into taxonomy table
# -------------------------------------------------------

rule mash_taxonomy:
    input:
        mash_results = expand(os.path.join(output_dir, "data", "mash", "{sample}.mash_screen.tsv"), sample=samples)
    output:
        summary = os.path.join(output_dir, "data", "mash", "mash_taxonomy.tsv")
    shell:
        """
        python workflow/scripts/mash_to_taxonomy.py {output.summary} {input.mash_results}
        """


# ---------------------------------
#   Predict sequence type (MLST)
# ---------------------------------

rule mlst:
    input:
        symlink = expand(os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"), sample=samples)
    output:
        mlst = os.path.join(output_dir, "data", "mlst", "mlst.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "mlst.txt")
    conda: "../envs/mlst.yaml"
    shell:
        """
        
        mlst $(dirname {input.symlink[0]})/*.fasta > {output.mlst}
        
        """
