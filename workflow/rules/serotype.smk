# ---------------------------------
#   Predict sequence type (MLST)
# ---------------------------------

rule mlst:
    input:
        assembly = expand(os.path.join(output_dir, "assemblies", "chromosomes", "{sample}.chromosome.fasta"), sample=samples)
    output:
        mlst = os.path.join(output_dir, "data", "mlst", "mlst.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "mlst.txt")
    conda: "../envs/mlst.yaml"
    shell:
        """
        # MLST
        mkdir -p $(dirname {output.mlst})
        mlst $(dirname {input.assembly[0]})/*.fasta > {output.mlst}
        """

# ---------------------------------
#   Serotype Salmonella (SeqSero2)
# ---------------------------------

rule seqsero:
    input:
        assembly = os.path.join(output_dir, "samples", "{sample}", "unicycler", "assembly.fasta")
    output:
        sq2 = os.path.join(output_dir, "data", "SeqSero2", "{sample}", "SeqSero_result.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.seqsero2.txt")
    conda: "../envs/seqsero.yaml"
    shell:
        """
        mkdir -p $(dirname {output.sq2})
        SeqSero2_package.py -m k -t 4 -i {input.assembly} -d $(dirname {output.sq2})
        """


# ----------------------------------------
#   Serotype/Pathotype E. coli (ECTyper)
# ----------------------------------------

#rule ectyper:
#    input:
#        assembly = os.path.join(output_dir, "assemblies", "chromosomes", "{sample}.chromosome.fasta")
#    output:
#        ect = os.path.join(output_dir, "data", "serotyping", "Ecoli", "{sample}", "output.tsv")
#    resources:
#        mem_mb = 20000,
#        time = "0-10:00:00",
#        threads = 8
#    benchmark:
#        os.path.join(output_dir, "data", "benchmarks", "{sample}.ectyper.txt")
#    conda: "envs/ectyper.yaml"
#    group: "group4"
#    shell:
#        """
#        # ECTyper
#        mkdir -p $(dirname {output.ect})
#        ectyper -i {input.assembly} -o $(dirname {output.ect}) --pathotype
#        """
