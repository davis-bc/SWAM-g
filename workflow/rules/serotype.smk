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

# ---------------------------------
#   Serotype Salmonella (SeqSero2)
# ---------------------------------

rule seqsero:
    input:
        r1_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz"),
        mash     = mash_taxonomy_file
    output:
        sq2 = os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "SeqSero_result.tsv")
    threads: 4
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.seqsero2.txt")
    conda: "../envs/seqsero.yaml"
    shell:
        r"""
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
            exit 0
        fi

        SeqSero2_package.py -m a \
                            -t 2 \
                            -p {threads} \
                            -i {input.r1_clean} {input.r2_clean} \
                            -d "$outdir" \
                            -n {wildcards.sample}
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
    threads: 4
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.sistr.txt")
    conda: "../envs/sistr.yaml"
    shell:
        r"""
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
            exit 0
        fi

        sistr_prefix="$outdir/sistr"
        sistr -i {input.assembly} {wildcards.sample} \
              -f tab \
              -o "$sistr_prefix" \
              --qc \
              -t {threads}

        mv "${{sistr_prefix}}.tab" {output.sistr}
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
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.txsscan.txt")
    conda: "../envs/macsyfinder.yaml"
    shell:
        """
        
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
        
        """


# --------------------------------------------
#   Serotype and pathotype E. coli (ECTyper)
# --------------------------------------------

rule ectyper:
    input:
        symlink = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        ecsetup = os.path.join(output_dir, "data", "serotype", "E.coli", ".ectyper_setup_done.txt")
    output:
        ect = os.path.join(output_dir, "data", "serotype", "E.coli", "{sample}", "output.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.ectyper.txt")
    conda: "../envs/ectyper.yaml"
    shell:
        """
        refseq="dbs/EnteroRef_GTDBSketch_20231003_V2.msh"

        # Copy sketch to local node scratch to avoid I/O contention on the shared
        # filesystem when many ECTyper jobs run simultaneously in a group job.
        LOCAL_MSH="${{TMPDIR:-/tmp}}/EnteroRef_GTDBSketch_20231003_V2.msh"
        cp "$refseq" "$LOCAL_MSH"

        ectyper -i {input.symlink} \
                -o $(dirname {output.ect}) \
                --pathotype \
                -r "$LOCAL_MSH" \
                --cores {threads} \
                --verify \
                --debug
        """        
