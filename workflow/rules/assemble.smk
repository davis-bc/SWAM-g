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

# ------------------------------------------
#         Assemble (Unicylcer)
# ------------------------------------------

rule unicylcer:
    input:
        r1_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")
    output:
        unicycler = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta")
    log:
        os.path.join(output_dir, "logs", "unicylcer", "{sample}.log")
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
            echo "[swamg-rule] unicylcer"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] unicylcer"
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

# -----------------------------------------
#    Calculate genomic coverage (minimap2)
# -----------------------------------------

rule coverage:
    input:
        fasta = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        r1 = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2 = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")
    output:
        bam = temp(os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}.bam")),
        coverage = os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}_coverage.tsv")
    log:
        os.path.join(output_dir, "logs", "coverage", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.coverage.txt")
    conda: "../envs/assemble.yaml"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        sort_tmp=$(mktemp -d)
        on_error() {{
            rc=$?
            rm -rf "$sort_tmp"
            echo "[swamg-rule] coverage"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] coverage"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        # Map reads using minimap2
        # Use a unique temp directory for samtools sort to avoid collisions on retries
        minimap2 -ax sr {input.fasta} {input.r1} {input.r2} | \
        samtools view -@ {threads} -bS - | \
        samtools sort -@ {threads} -T "$sort_tmp/sort" -o {output.bam}
        rm -rf "$sort_tmp"

        # Index BAM for coverage calculation
        samtools index {output.bam}

        # Calculate genome coverage with samtools depth
        cov=$(samtools depth {output.bam} | awk '{{sum+=$3; cnt++}} END {{print (cnt > 0 ? sum/cnt : 0)}}')

        # Write result as TSV: sample | coverage
        echo -e "{wildcards.sample}\t$cov" > {output.coverage}

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """
    
