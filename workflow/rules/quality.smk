# ------------------------------------------
#   Quality assessment and QC metrics
# ------------------------------------------

# Rules in this module should quantify assembly or read quality and emit QC
# metrics that inform interpretation, filtering, or reporting.

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

# ---------------------------------------------------------
#        Check completeness and contamination (CheckM2)
# ---------------------------------------------------------

rule checkm2:
    input:
        symlink = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        checkm =  os.path.join(output_dir, "data", "checkm2", ".checkm_initialized")
    output:
        checkm = os.path.join(output_dir, "data", "checkm2", "{sample}", "quality_report.tsv")
    log:
        os.path.join(output_dir, "logs", "checkm2", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.checkm2.txt")
    conda: "../envs/checkm2.yaml"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] checkm2"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] checkm2"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        DB="dbs/CheckM2_database/uniref100.KO.1.dmnd"
        
        # Run CheckM2 on assemblies
        checkm2 predict \
        -i {input.symlink} \
        -o $(dirname {output.checkm}) \
        -x fasta \
        --threads {threads} \
        --database_path "$DB" \
        --force 

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """
