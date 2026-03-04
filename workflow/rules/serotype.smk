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
        assembly = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta")
    output:
        sq2 = os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "SeqSero_result.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.seqsero2.txt")
    conda: "../envs/seqsero.yaml"
    shell:
        """
       
        SeqSero2_package.py -m k \
                            -t 4 \
                            -i {input.assembly} \
                            -d $(dirname {output.sq2}) 
        
        # -m 'k'(raw reads and genome assembly k-mer)
        # -t '4' for genome assembly
        
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
        
        models_dir="dbs/macsyfinder/models/"
        
        # find all proteins in assembly
        prodigal -i {input.assembly} -a {output.prot} -q 
        
        # run TXSScan on protiens
        macsyfinder --db-type unordered \
                    --sequence-db {output.prot} \
                    --models TXSScan all \
                    --models-dir $models_dir \
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
    threads: 1
    resources:
        mem_mb = 5000,
        time = "1d"
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

