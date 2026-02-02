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
            ln -sf "$input_assembly" "$output_symlink"
        done
        
        """

# ---------------------------------------------------------
#        Check completeness and contamination (CheckM2)
# ---------------------------------------------------------

rule checkm2:
    input:
        symlink = expand(os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"), sample=samples),
        checkm =  os.path.join(output_dir, "data", "checkm2", ".checkm_initialized")
    output:
        checkm = os.path.join(output_dir, "data", "checkm2", "quality_report.tsv")
    resources:
        threads = 32,
        mem_mb = 100000,
        time = "1d"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "checkm2.txt")
    conda: "../envs/classify.yaml"
    shell:
        """
        DB="dbs/CheckM2_database/uniref100.KO.1.dmnd"
        
        # Run CheckM2 on assemblies
        checkm2 predict \
        -i $(dirname {input.symlink[0]}) \
        -o $(dirname {output.checkm}) \
        -x fasta \
        --threads {resources.threads} \
        --database_path "$DB" \
        --force
        
        """


# ------------------------------------------------
#            Classify taxonomy (GTDB-tk)
# ------------------------------------------------

rule gtdbtk:
    input:
        symlink = expand(os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"), sample=samples),
        gtdb = os.path.join(output_dir, "data", "gtdb-tk", ".gtdb_initialized")
    output:
        gtdbtk = protected(os.path.join(output_dir, "data", "gtdb-tk", "gtdbtk.bac120.summary.tsv"))
    resources:
        threads = 32, 
        mem_mb = 150000,
        time = "1d"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "gtdbtk.txt")
    conda: "../envs/classify.yaml"
    shell:
        """
        # Dynamically find the GTDB release directory
        GTDB_PATH=$(find dbs -type d -name 'release*' -maxdepth 1)
        
        # Ensure the GTDB release directory exists
        if [ ! -d "$GTDB_PATH" ]; then
            echo "Error: GTDB release directory not found in dbs/."
            exit 1
        fi

        # Export GTDB data path for GTDB-tk
        export GTDBTK_DATA_PATH=$GTDB_PATH

        echo "Using GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"

        # Run GTDB-tk classify workflow
        gtdbtk classify_wf \
        --genome_dir $(dirname {input.symlink[0]}) \
        --out_dir $(dirname {output.gtdbtk}) \
        -x fasta \
        --cpus {resources.threads} \
        --skip_ani_screen
        
        """


"""


# ------------------------------------------------
#    Compute pairwise ANI distance matrix (FastANI)
# ------------------------------------------------

rule fastani:
    input:
        assemblies = expand(os.path.join(output_dir, "data", "unicycler", "{sample}", "assembly.fasta"), sample=samples)
    output:
        pairwise_matrix = os.path.join(output_dir, "pairwise_ani_matrix.tsv")
    resources:
        mem_mb = 80000,
        time = "1-00:00:00",
        threads = 32
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "fastani.txt")
    conda: "../envs/gtdb.yaml"
    shell:
        
        # Generate a list of input assemblies
        assembly_list=$(mktemp)
        printf "%s\\n" {input.assemblies} > $assembly_list

        # Run FastANI with the --matrix flag to compute the pairwise ANI distance matrix
        fastANI --ql $assembly_list --rl $assembly_list --matrix -o {output.pairwise_matrix} -t {resources.threads}

        # Clean up temporary file
        rm -f $assembly_list
        
        

"""