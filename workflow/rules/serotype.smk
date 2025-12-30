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
        
        mlst $(dirname {input.assembly[0]})/*.fasta --quiet > {output.mlst}
        
        """

# ---------------------------------
#   Serotype Salmonella (SeqSero2)
# ---------------------------------

rule seqsero:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta")
    output:
        sq2 = os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "SeqSero_result.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.seqsero2.txt")
    conda: "../envs/seqsero.yaml"
    shell:
        """
        
        SeqSero2_package.py -m k -t 4 -i {input.assembly} -d $(dirname {output.sq2}) > /dev/null 2>&1
        
        """


# --------------------------------------------
#   Serotype and pathotype E. coli (ECTyper)
# --------------------------------------------

rule ectyper:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        ecsetup = os.path.join(output_dir, "data", "serotype", "E.coli", ".ectyper_setup_done.txt")
    output:
        ect = os.path.join(output_dir, "data", "serotype", "E.coli", "{sample}", "output.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.ectyper.txt")
    conda: "../envs/ectyper.yaml"
    shell:
        """

        ectyper -i {input.assembly} -o $(dirname {output.ect}) --pathotype -r $(dirname {input.ecsetup})/refseq.genomes.k21s1000.msh > /dev/null 2>&1
        
        """

# -------------------------------------------------
#   Scan for bacterial secretion systems (TXSScan)
# -------------------------------------------------

rule txsscan:
    input:
        assembly = os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"),
        models =   os.path.join(output_dir, "data", "txsscan", ".txsscan.setup.done")
    output:
        prot = os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}.prot.faa"),
        sys  = os.path.join(output_dir, "data", "txsscan", "{sample}", "all_systems.tsv")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.txsscan.txt")
    conda: "../envs/macsyfinder.yaml"
    shell:
        """
        models_dir="dbs/macsyfinder/models/"
        
        # find all proteins in assembly
        prodigal -i {input.assembly} -a {output.prot} -q > /dev/null 2>&1
        
        # run TXSScan on protiens
        macsyfinder --db-type unordered --sequence-db {output.prot} --models TXSScan all --models-dir $models_dir -o $(dirname {output.sys}) --force > /dev/null 2>&1
        
        """


