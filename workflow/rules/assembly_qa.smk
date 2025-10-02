# ------------------------------------
#    QA complete assemblies (Quast)
# ------------------------------------

rule assembly_qa:
    input:
        assemblies = expand(os.path.join(output_dir, "samples", "{sample}", "spades", "contigs.fasta"), sample=samples)
    output:
        quast = os.path.join(output_dir, "data", "quast", "report.tsv")
    conda: "../envs/quast.yaml"
    shell:
        r"""
        # Create temp directory for symlinks
        tmp_quast_dir=$(mktemp -d $(dirname {output.quast})/quast_input.XXXXXXXX)
        for infile in {input.assemblies}; do
            sample=$(basename $(dirname $(dirname $infile)))
            ln -sf $(realpath $infile) $tmp_quast_dir/${{sample}}.fasta
        done

        mkdir -p $(dirname {output.quast})
        quast.py $tmp_quast_dir/*.fasta -o $(dirname {output.quast})

        # Clean up temp symlink dir
        rm -rf $tmp_quast_dir
        """
