# -------------------------------------------------------
#   Screen assembly against EnteroRef MASH sketch (MASH)
# -------------------------------------------------------

rule mash_classify:
    input:
        assembly  = os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"),
        msh_ready = os.path.join(output_dir, "data", "serotype", "E.coli", ".ectyper_setup_done.txt")
    output:
        mash_result = os.path.join(output_dir, "data", "mash", "{sample}.mash_screen.tsv")
    resources:
        mem_mb = 8000,
        time = "0-01:00:00",
        threads = 4
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.mash_classify.txt")
    conda: "../envs/mash.yaml"
    shell:
        """
        mkdir -p $(dirname {output.mash_result})
        mash screen -p {resources.threads} -w dbs/EnteroRef_GTDBSketch_20231003_V2.msh {input.assembly} \
            | sort -grk1,1 > {output.mash_result}
        """

# -------------------------------------------------------
#   Combine per-sample MASH results into taxonomy table
# -------------------------------------------------------

rule mash_taxonomy:
    input:
        mash_results = expand(os.path.join(output_dir, "data", "mash", "{sample}.mash_screen.tsv"), sample=samples)
    output:
        summary = os.path.join(output_dir, "data", "mash", "mash_taxonomy.tsv")
    shell:
        """
        python workflow/scripts/mash_to_taxonomy.py {output.summary} {input.mash_results}
        """

# ------------------------------------------------
#        Create symlink database for assemblies
# ------------------------------------------------

rule symlink:
    input:
        assembly = expand(os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"), sample=samples)
    output:
        symlink = expand(os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"), sample=samples)
    shell:
        """
        
        # Create symbolic links for all assemblies in the batch directory
        for i in $(seq 0 $(expr $(echo {input.assembly} | tr -s ' ' '\n' | wc -l) - 1)); do
            input_assembly=$(echo {input.assembly} | cut -d' ' -f$(expr $i + 1))
            output_symlink=$(echo {output.symlink} | cut -d' ' -f$(expr $i + 1))
            ln -sf "$(readlink -f "$input_assembly")" "$output_symlink"
        done
        
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
    threads: 1
    resources:
        mem_mb = 2000,
        time = "6d"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.checkm2.txt")
    conda: "../envs/checkm2.yaml"
    shell:
        """
        DB="dbs/CheckM2_database/uniref100.KO.1.dmnd"
        
        #tmp="$(dirname {output.checkm})/tmp"
        #mkdir -p $tmp
        
        # Delay before multiprocessing.Manager attempts to create socket files
        #sleep 5
        
        # Run CheckM2 on assemblies
        checkm2 predict \
        -i {input.symlink} \
        -o $(dirname {output.checkm}) \
        -x fasta \
        --threads {threads} \
        --database_path "$DB" \
        --force 
        
        """


