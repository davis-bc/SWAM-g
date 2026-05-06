# ------------------------------------------
#   Serotype- and pathotype-specific typing
# ------------------------------------------

# Rules in this module should run specialized serotype or pathotype tools for
# biologically scoped taxa, such as Salmonella or E. coli.

# ---------------------------------
#   Serotype Salmonella (SeqSero2)
# ---------------------------------

rule seqsero2:
    input:
        r1_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz"),
        mash     = mash_taxonomy_file
    output:
        sq2 = os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "SeqSero_result.tsv")
    log:
        os.path.join(output_dir, "logs", "seqsero2", "{sample}.log")
    threads: 4
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.seqsero2.txt")
    conda: "../envs/seqsero.yaml"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] seqsero2"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] seqsero2"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        outdir=$(dirname {output.sq2})
        mkdir -p "$outdir"

        is_salmonella=$(python3 - "{input.mash}" "{wildcards.sample}" <<'PY'
import csv
import sys

mash_tsv, sample = sys.argv[1], sys.argv[2]
sample_key = sample + ".chromosome"

with open(mash_tsv) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["user_genome"] != sample_key:
            continue
        print("1" if "g__Salmonella" in row["classification"] else "0")
        break
    else:
        raise SystemExit("Sample '%s' missing from %s" % (sample, mash_tsv))
PY
)

        if [ "$is_salmonella" != "1" ]; then
            printf 'Sample name\tOutput directory\tPredicted identification\tPredicted antigenic profile\tPredicted serotype\tWorkflow status\n' > {output.sq2}
            printf '{wildcards.sample}\t%s\tMASH non-Salmonella skip\t\t\tSKIPPED_NOT_SALMONELLA\n' "$outdir" >> {output.sq2}
            trap - ERR
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] SKIPPED"
            echo "[swamg-skip-reason] NOT_SALMONELLA"
            exit 0
        fi

        SeqSero2_package.py -m a \
                            -t 2 \
                            -p {threads} \
                            -i {input.r1_clean} {input.r2_clean} \
                            -d "$outdir" \
                            -n {wildcards.sample}

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """


# -------------------------------
#   Serotype Salmonella (SISTR)
# -------------------------------

rule sistr:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        mash     = mash_taxonomy_file
    output:
        sistr = os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "sistr.tsv")
    log:
        os.path.join(output_dir, "logs", "sistr", "{sample}.log")
    threads: 4
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.sistr.txt")
    conda: "../envs/sistr.yaml"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] sistr"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] sistr"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        outdir=$(dirname {output.sistr})
        mkdir -p "$outdir"

        is_salmonella=$(python3 - "{input.mash}" "{wildcards.sample}" <<'PY'
import csv
import sys

mash_tsv, sample = sys.argv[1], sys.argv[2]
sample_key = sample + ".chromosome"

with open(mash_tsv) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["user_genome"] != sample_key:
            continue
        print("1" if "g__Salmonella" in row["classification"] else "0")
        break
    else:
        raise SystemExit("Sample '%s' missing from %s" % (sample, mash_tsv))
PY
)

        if [ "$is_salmonella" != "1" ]; then
            printf 'genome\tserovar\tserovar_antigen\tserovar_cgmlst\tcgmlst_ST\tserogroup\to_antigen\th1\th2\tqc_status\tqc_messages\n' > {output.sistr}
            printf '{wildcards.sample}\t\t\t\t\t\t\t\t\tSKIPPED_NOT_SALMONELLA\tMASH classified sample as non-Salmonella\n' >> {output.sistr}
            trap - ERR
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] SKIPPED"
            echo "[swamg-skip-reason] NOT_SALMONELLA"
            exit 0
        fi

        sistr_prefix="$outdir/sistr"
        sistr -i {input.assembly} {wildcards.sample} \
              -f tab \
              -o "$sistr_prefix" \
              --qc \
              -t {threads}

        mv "${{sistr_prefix}}.tab" {output.sistr}

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """

# --------------------------------------------
#   Serotype and pathotype E. coli (ECTyper)
# --------------------------------------------

rule ectyper:
    input:
        symlink = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        ecsetup = enteroref_sketch_ready,
        mash    = mash_taxonomy_file
    output:
        ect = os.path.join(output_dir, "data", "serotype", "E.coli", "{sample}", "output.tsv")
    log:
        os.path.join(output_dir, "logs", "ectyper", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.ectyper.txt")
    conda: "../envs/ectyper.yaml"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] ectyper"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] ectyper"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        refseq="dbs/EnteroRef_GTDBSketch_20231003_V2.msh"
        outdir=$(dirname {output.ect})

        mkdir -p "$outdir"

        is_ecoli=$(python3 - "{input.mash}" "{wildcards.sample}" <<'PY'
import csv
import sys

mash_tsv, sample = sys.argv[1], sys.argv[2]
sample_key = sample + ".chromosome"

with open(mash_tsv) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        if row["user_genome"] != sample_key:
            continue
        classification = row["classification"]
        print("1" if "s__Escherichia coli" in classification else "0")
        break
    else:
        raise SystemExit("Sample '%s' missing from %s" % (sample, mash_tsv))
PY
)

        if [ "$is_ecoli" != "1" ]; then
            printf 'Name\tSpecies\tQC\tWarnings\n' > {output.ect}
            printf '{wildcards.sample}\tMASH non-E. coli skip\tSKIPPED_NOT_ECOLI\tPipeline MASH taxonomy did not classify this sample as Escherichia coli\n' >> {output.ect}
            trap - ERR
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] SKIPPED"
            echo "[swamg-skip-reason] NOT_ECOLI"
            exit 0
        fi

        # Copy sketch to local node scratch to avoid I/O contention on the shared
        # filesystem when many ECTyper jobs run simultaneously in a group job.
        LOCAL_MSH="${{TMPDIR:-/tmp}}/EnteroRef_GTDBSketch_20231003_V2.msh"
        cp "$refseq" "$LOCAL_MSH"

        ectyper -i {input.symlink} \
                -o "$outdir" \
                --pathotype \
                -r "$LOCAL_MSH" \
                --cores {threads} \
                --debug

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """        
