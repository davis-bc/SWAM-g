# ---------------------------------------------------------
#        Check completeness and contamination (CheckM)
# ---------------------------------------------------------

rule checkm:
    input:
        assembly = expand(os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"), sample=samples),
        checkm = os.path.join(output_dir, "data", "checkm", ".checkm_initialized")
    output:
        checkm = protected(os.path.join(output_dir, "data", "checkm", "genome.stats.tsv")),
        checkm_stats = os.path.join(output_dir, "data", "checkm", "storage", "bin_stats.analyze.tsv")
    resources:
        mem_mb = 100000,
        time = "1-00:00:00",
        threads = 32
    threads: 32
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "checkm.txt")
    conda: "../envs/checkm.yaml"
    shell:
        """
        # Check for the checkm_data directory in dbs/
        CHECKM_DB_PATH=$(find dbs -type d -name 'checkm_data' -maxdepth 1)
        
        # Ensure the checkm_data directory exists
        if [ ! -d "$CHECKM_DB_PATH" ]; then
            echo "Error: CHECKM_DATA_PATH (dbs/checkm_data) not found."
            exit 1
        fi

        # Set CheckM database path environment variable
        export CHECKM_DATA_PATH=$CHECKM_DB_PATH

        echo "Using CHECKM_DATA_PATH=$CHECKM_DATA_PATH"

        # Create output directories
        mkdir -p $(dirname {output.checkm})
        mkdir -p $(dirname {output.checkm_stats})

        # Run CheckM lineage workflow
        checkm lineage_wf $(dirname {input.assembly[0]}) $(dirname {output.checkm}) -x fasta -t {resources.threads} --tmpdir $(dirname {output.checkm})
        
        # Generate QA report
        checkm qa $(dirname {output.checkm})/lineage.ms $(dirname {output.checkm}) -f {output.checkm}
        
        """


# ------------------------------------------------
#            Classify taxonomy (GTDB-tk)
# ------------------------------------------------

rule gtdbtk:
    input:
        assembly = expand(os.path.join(output_dir, "data", "assemblies", "chromosomes", "{sample}.chromosome.fasta"), sample=samples),
        gtdb = os.path.join(output_dir, "data", "gtdb-tk", ".gtdb_initialized")
    output:
        gtdbtk = protected(os.path.join(output_dir, "data", "gtdb-tk", "gtdbtk.bac120.summary.tsv"))
    resources:
        mem_mb = 150000,
        time = "1-00:00:00",
        threads = 32
    threads: 32
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "gtdbtk.txt")
    conda: "../envs/gtdb.yaml"
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

        # Create output directory
        mkdir -p $(dirname {output.gtdbtk})

        # Run GTDB-tk classify workflow
        gtdbtk classify_wf --genome_dir $(dirname {input.assembly[0]}) --out_dir $(dirname {output.gtdbtk}) -x fasta --cpus {resources.threads} --skip_ani_screen
        
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
        
        # Create output directory
        mkdir -p $(dirname {output.pairwise_matrix})

        # Generate a list of input assemblies
        assembly_list=$(mktemp)
        printf "%s\\n" {input.assemblies} > $assembly_list

        # Run FastANI with the --matrix flag to compute the pairwise ANI distance matrix
        fastANI --ql $assembly_list --rl $assembly_list --matrix -o {output.pairwise_matrix} -t {resources.threads}

        # Clean up temporary file
        rm -f $assembly_list
        
        

"""