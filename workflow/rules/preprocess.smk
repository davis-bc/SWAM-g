# ------------------------------------------
#   Read preprocessing and input cleanup
# ------------------------------------------

# Rules in this module should prepare raw inputs for downstream assembly or
# analysis without changing the biological interpretation of the sample.

# ------------------------------------------
#        Clean input data (fastp)
# ------------------------------------------

rule fastp:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")
    log:
        os.path.join(output_dir, "logs", "fastp", "{sample}.log")
    params:
        html = "/dev/null/",
        json = "/dev/null/"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.fastp.txt")
    conda: "../envs/assemble.yaml"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] fastp"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] fastp"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        fastp -i {input.r1} \
              -I {input.r2} \
              -o {output.r1_clean} \
              -O {output.r2_clean} \
              --html {params.html} \
              --json {params.json}

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """
