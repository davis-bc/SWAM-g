# ------------------------------------------
#   Assembly retry resource scaling
# ------------------------------------------

# Benchmark-guided retry tiers: the tagged benchmark set peaked just below
# 30 GB RSS for Unicycler, so start near 40 GB and scale upward on retries.
UNICYCLER_MEM_MB_BY_ATTEMPT = (40000, 60000, 80000)
UNICYCLER_RUNTIME_BY_ATTEMPT = ("4h", "6h", "8h")


def unicycler_mem_mb(_, attempt):
    return UNICYCLER_MEM_MB_BY_ATTEMPT[min(attempt, len(UNICYCLER_MEM_MB_BY_ATTEMPT)) - 1]


def unicycler_runtime(_, attempt):
    return UNICYCLER_RUNTIME_BY_ATTEMPT[min(attempt, len(UNICYCLER_RUNTIME_BY_ATTEMPT)) - 1]


# ------------------------------------------
#   Assembly generation and derived inputs
# ------------------------------------------

# Rules in this module should create assemblies or assembly-derived shared
# artifacts that feed multiple downstream analyses.

# ------------------------------------------
#         Assemble (Unicycler)
# ------------------------------------------

rule unicycler:
    input:
        r1_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")
    output:
        unicycler = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta")
    log:
        os.path.join(output_dir, "logs", "unicycler", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.unicycler.txt")
    resources:
        mem_mb = unicycler_mem_mb,
        runtime = unicycler_runtime
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
            echo "[swamg-rule] unicycler"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] unicycler"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        # Unicycler 0.5.x renamed the SPAdes passthrough flag from
        # --spades_opts to --spades_options.
        spades_mem_gb=$(( ({resources.mem_mb} + 1023) / 1024 ))
        echo "[swamg-threads] {threads}"
        echo "[swamg-mem-mb] {resources.mem_mb}"
        echo "[swamg-runtime-min] {resources.runtime}"
        echo "[swamg-spades-mem-gb] $spades_mem_gb"

        unicycler -1 {input.r1_clean} \
                  -2 {input.r2_clean} \
                  -o $(dirname {output.unicycler}) \
                  -t {threads} \
                  --spades_options "-m $spades_mem_gb" \
                  --keep 0

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """

# ------------------------------------------------
#        Create symlink database for assemblies
# ------------------------------------------------

rule symlink:
    input:
        assembly = expand(os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"), sample=samples)
    output:
        symlink = expand(os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"), sample=samples)
    shell:
        """
        
        # Create symbolic links for all assemblies in the batch directory
        for i in $(seq 0 $(expr $(echo {input.assembly} | tr -s ' ' '\n' | wc -l) - 1)); do
            input_assembly=$(echo {input.assembly} | cut -d' ' -f$(expr $i + 1))
            output_symlink=$(echo {output.symlink} | cut -d' ' -f$(expr $i + 1))
            ln -sf "$(readlink -f "$input_assembly")" "$output_symlink"
        done
        
        """
