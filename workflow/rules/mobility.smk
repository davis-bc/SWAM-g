# ------------------------------------------
#   Mobility and mobile feature detection
# ------------------------------------------

# Rules in this module should reconstruct plasmids or detect mobile, secreted,
# or mobility-associated genomic features.

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

# ------------------------------------------------------
#     Comprehensively screen MGEs (MobileElementFinder)
# ------------------------------------------------------

rule mobileelementfinder:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta")
    output:
        mef        = os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.csv"),
        temp_fasta = temp(os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.tmp.fasta"))
    log:
        os.path.join(output_dir, "logs", "mobileelementfinder", "{sample}.log")
    resources:
        mef_slots = 1
    conda: "../envs/mobileelementfinder.yaml"
    params:
        temp_dir = os.path.join(output_dir, "data", "mobileelementfinder", "{sample}")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mobileelementfinder.txt")
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        tmp_root=$(mktemp -d "{params.temp_dir}/tmp.XXXXXX")
        on_error() {{
            rc=$?
            rm -rf "$tmp_root"
            echo "[swamg-rule] mobileelementfinder"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] mobileelementfinder"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        # MobileElementFinder uses Python's tempdir to place its transient BLAST
        # database under $TMPDIR/mge_finder/database. Give each sample its own
        # temp root to avoid parallel makeblastdb collisions on shared /tmp.
        export TMPDIR="$tmp_root"

        # Rebuild a clean FASTA for MobileElementFinder. This strips a UTF-8 BOM,
        # normalizes line endings, trims header descriptions to the first token,
        # and preserves only sequence data BLAST will accept.
        python - <<'PY' "{input.assembly}" "{output.temp_fasta}"
import sys

input_path, output_path = sys.argv[1:3]
allowed_bases = set("ACGTRYSWKMBDHVN.-*")


def write_record(handle, header, sequence_lines):
    if header is None:
        return

    sequence = "".join(sequence_lines)
    if not sequence:
        raise SystemExit("Empty FASTA sequence for %s in %s" % (header, input_path))

    handle.write(header + "\n")
    for i in range(0, len(sequence), 80):
        handle.write(sequence[i:i + 80] + "\n")


with open(input_path, encoding="utf-8-sig") as src, open(output_path, "w", encoding="ascii", newline="\n") as dst:
    header = None
    sequence_lines = []
    saw_header = False

    for line_number, raw_line in enumerate(src, start=1):
        line = raw_line.rstrip("\r\n")
        if not line:
            continue
        if line.startswith(";"):
            continue

        if line.startswith(">"):
            write_record(dst, header, sequence_lines)
            identifier = line[1:].strip().split(None, 1)[0] if line[1:].strip() else ""
            if not identifier:
                raise SystemExit("Missing FASTA identifier on line %d in %s" % (line_number, input_path))
            header = ">" + identifier
            sequence_lines = []
            saw_header = True
            continue

        if not saw_header:
            raise SystemExit(
                "Encountered FASTA sequence data before the first header on line %d in %s"
                % (line_number, input_path)
            )

        sequence = "".join(line.split()).upper()
        invalid_bases = sorted(set(sequence) - allowed_bases)
        if invalid_bases:
            raise SystemExit(
                "Invalid FASTA character(s) %s on line %d in %s"
                % ("".join(invalid_bases), line_number, input_path)
            )
        sequence_lines.append(sequence)

    if saw_header:
        write_record(dst, header, sequence_lines)
PY

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

        rm -rf "$tmp_root"
        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """

# -------------------------------------------------
#   Scan for bacterial secretion systems (TXSScan)
# -------------------------------------------------

rule txsscan:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        models =   os.path.join(output_dir, "data", "txsscan", ".txsscan.setup.done")
    output:
        prot = os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}.prot.faa"),
        sys  = os.path.join(output_dir, "data", "txsscan", "{sample}", "all_systems.tsv"),
        mods = os.path.join(output_dir, "data", "txsscan", "{sample}", "all_systems.txt")
    log:
        os.path.join(output_dir, "logs", "txsscan", "{sample}.log")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.txsscan.txt")
    conda: "../envs/macsyfinder.yaml"
    shell:
        r"""
        log_dir=$(dirname "{log}")
        mkdir -p "$log_dir"
        exec > "{log}" 2>&1
        set -euo pipefail

        started_at=$(date -Is)
        on_error() {{
            rc=$?
            echo "[swamg-rule] txsscan"
            echo "[swamg-sample] {wildcards.sample}"
            echo "[swamg-host] $(hostname)"
            echo "[swamg-started-at] $started_at"
            echo "[swamg-finished-at] $(date -Is)"
            echo "[swamg-status] FAILED"
            echo "[swamg-exit-code] $rc"
            exit $rc
        }}
        trap 'on_error' ERR

        echo "[swamg-rule] txsscan"
        echo "[swamg-sample] {wildcards.sample}"
        echo "[swamg-host] $(hostname)"
        echo "[swamg-started-at] $started_at"

        models_dir="dbs/macsyfinder/models"

        if [ ! -d "$models_dir/TXSScan" ]; then
            echo "TXSScan models are missing from $models_dir. Run txsscan_init on an internet-connected login/head node." >&2
            exit 1
        fi
        
        # find all proteins in assembly
        prodigal -i {input.assembly} -a {output.prot} -q 
        
        # run TXSScan on proteins using the pre-staged offline model bundle.
        # Upstream TXSScan examples use ordered_replicon mode for protein FASTAs.
        macsyfinder --db-type ordered_replicon \
                    --sequence-db {output.prot} \
                    --models TXSScan all \
                    --models-dir "$models_dir" \
                    -o $(dirname {output.sys}) \
                    --force 

        trap - ERR
        echo "[swamg-finished-at] $(date -Is)"
        echo "[swamg-status] SUCCESS"
        """
