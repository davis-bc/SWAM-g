# ---------------------------------
#   Predict sequence type (MLST)
# ---------------------------------

rule mlst:
    input:
        assembly = expand(os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"), sample=samples)
    output:
        mlst = os.path.join(output_dir, "data", "mlst", "mlst.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "mlst.txt")
    conda: "../envs/mlst.yaml"
    shell:
        """
        
        # MLST
        mkdir -p $(dirname {output.mlst})
        mlst $(dirname {input.assembly[0]})/*.fasta --quiet > {output.mlst}
        
        """

# ---------------------------------
#   Serotype Salmonella (SeqSero2)
# ---------------------------------

rule seqsero:
    input:
        assembly = os.path.join(output_dir, "data", "unicylcer", "{sample}", "assembly.fasta")
    output:
        sq2 = os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "SeqSero_result.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.seqsero2.txt")
    conda: "../envs/seqsero.yaml"
    shell:
        """
        
        mkdir -p $(dirname {output.sq2})
        SeqSero2_package.py -m k -t 4 -i {input.assembly} -d $(dirname {output.sq2})
        
        """


# ----------------------------------------
#     ECTyper install and database setup
# ----------------------------------------

rule setup_ectyper:
    output:
        ecsetup = os.path.join(output_dir, "bin", ".ectyper_setup_done.txt")
    conda:
        "../envs/ectyper.yaml"
    shell:
        """
        mkdir -p $(dirname {output.ecsetup}) && cd $(dirname {output.ecsetup})
        
        # Clone Ectyper if not already present
        if [ ! -d ecoli_serotyping ]; then
            git clone https://github.com/phac-nml/ecoli_serotyping.git
            cd ecoli_serotyping
            git checkout v2.0.0
        else
            cd ecoli_serotyping
        fi

        # Run pip install and capture output
        INSTALL_LOG=$(pip install . 2>&1)

        # Check for success message and create placeholder file accordingly
        echo "$INSTALL_LOG" | grep "Successfully built ectyper" && \
            echo "Ectyper setup completed successfully" > {output.ecsetup} && \
            echo "Ectyper setup was successful. Environment and software are ready to use." || \
            (echo "$INSTALL_LOG"; echo "Ectyper setup failed!"; exit 1)
        """

# --------------------------------------------
#   Serotype and pathotype E. coli (ECTyper)
# --------------------------------------------

rule ectyper:
    input:
        assembly = os.path.join(output_dir, "data", "unicylcer", "{sample}", "assembly.fasta"),
        ecsetup = os.path.join(output_dir, "bin", ".ectyper_setup_done.txt")
    output:
        ect = os.path.join(output_dir, "data", "serotype", "E.coli", "{sample}", "output.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.ectyper.txt")
    conda: "../envs/ectyper.yaml"
    shell:
        """
        
        # ECTyper
        mkdir -p $(dirname {output.ect})
        ectyper -i {input.assembly} -o $(dirname {output.ect}) --pathotype
        
        """

# --------------------------
#    Initialize MacSyFinder
# --------------------------

rule msf_init:
    output:
        models = os.path.join(output_dir, "bin", "macsyfinder", "models", ".msf.setup.done")
    conda:
        "../envs/macsyfinder.yaml"
    shell:
        """
        
        MSF_DB="$(dirname {output.models})"
        if [ ! -d "$MSF_DB" ]; then
        echo "TXSScan models do not exist, initializing..."
        
        msf_data install TXSScan --target $(dirname {output.models})
        
        fi
        
        touch {output}
        
        """

# -------------------------------------------------
#   Scan for bacterial secretion systems (TXSScan)
# -------------------------------------------------

rule txsscan:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        models =   os.path.join(output_dir, "bin", "macsyfinder", "models", ".msf.setup.done")
    output:
        prot = os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}.prot.faa"),
        sys  = os.path.join(output_dir, "data", "txsscan", "{sample}", "all_systems.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.txsscan.txt")
    conda: "../envs/macsyfinder.yaml"
    shell:
        """
        # find all proteins in assembly
        prodigal -i {input.assembly} -a {output.prot} -q > /dev/null
        
        # run TXSScan on protiens
        macsyfinder --db-type unordered --sequence-db {output.prot} --models TXSScan all --models-dir $(dirname {input.models}) -o $(dirname {output.sys}) > /dev/null
        
        """


