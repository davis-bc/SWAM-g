# ----------------------------------------
#               Setup GTDB
# ---------------------------------------

rule gtdb_init:
    output:
        touch(os.path.join(output_dir, "data", "gtdb-tk", ".gtdb_initialized"))
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_gtdb.txt")
    shell:
        """
        DB_DIR="dbs"
        GTDB_URL="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz"

        cd $DB_DIR

        # Check if any 'release' directory exists; if not, download and extract the database
        if [ ! -d "$(find . -type d -name 'release*' -maxdepth 1)" ]; then
            echo "GTDB release directory not found. Downloading database..."
            wget -q $GTDB_URL -O gtdbtk_data.tar.gz > /dev/null 2>&1
            echo "Extracting database..."
            tar -xvzf gtdbtk_data.tar.gz > /dev/null 2>&1
            rm gtdbtk_data.tar.gz  # Clean up after extraction
        fi

        # Mark the initialization by touching the output file
        touch {output}
        """

# ----------------------------------------
#               Setup checkM
# ---------------------------------------

rule checkm_init:
    output:
        touch(os.path.join(output_dir, "data", "checkm", ".checkm_initialized"))
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_checkm.txt")
    shell:
        """
        DB_DIR="dbs"
        checkm_URL="https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"

        cd $DB_DIR
        
        if [ ! -d checkm_data ]; then
            echo "checkm_data not found. Downloading database..."
            mkdir checkm_data && cd checkm_data
            wget -q $checkm_URL 
            tar -xvf checkm_data_2015_01_16.tar.gz > /dev/null 2>&1
            rm checkm_data_2015_01_16.tar.gz
        fi

        touch {output}
        
        """

# -----------------------------------
#    Initialize MOB-suite database
# -----------------------------------

rule mob_init:
    output:
        touch(os.path.join(output_dir, "data", "mob-suite", ".mob_suite_db_initialized"))
    conda:
        "../envs/mob_suite.yaml"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_mobsuite.txt")
    shell:
        """
        MOB_DB=$(find .snakemake/conda -type d -path "*/mob_suite/databases" | head -n 1)
        if [ -z "$MOB_DB" ]; then
        echo "mob_suite database does not exist in conda environment, initializing..."
            mob_init > /dev/null 2>&1
        fi
        
        touch {output}
        
        """
        
# --------------------------------------
#   Initialize AMRFinderPlus database
# --------------------------------------

rule afp_init:
    output:
        touch(os.path.join(output_dir, "data", "amrfinderplus", ".afp_db_initialized"))
    conda:
        "../envs/amrfinder.yaml"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_afp.txt")
    shell:
        """
        amrfinder -u 
        
        touch {output}
        
        """
        
# ----------------------------------------
#    Setup ResFinder databases
# ---------------------------------------

rule res_init:
    output:
        touch(os.path.join(output_dir, "data", "resfinder", ".res_db_initialized"))
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_resdb.txt")
    shell:
        """
        DB_DIR="dbs"
        
        cd $DB_DIR

        # Check and clone only if the database directories are missing
        [ ! -d "resfinder_db" ] && git clone https://bitbucket.org/genomicepidemiology/resfinder_db/ > /dev/null 2>&1
        [ ! -d "pointfinder_db" ] && git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/ > /dev/null 2>&1
        [ ! -d "disinfinder_db" ] && git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/ > /dev/null 2>&1

        # Mark the initialization by touching the output file
        
        touch {output}
        
        """    

    

# ----------------------------------------
#     ECTyper database setup
# ----------------------------------------

rule ectyper_init:
    output:
        touch(os.path.join(output_dir, "data", "serotype", "E.coli", ".ectyper_setup_done.txt"))
    conda:
        "../envs/ectyper.yaml"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_ectyper.txt")
    shell:
        """
        DB_DIR="dbs"
        
        cd $DB_DIR

        # Check if the MASH sketch exists, download if missing
        if [ ! -f refseq.genomes.k21s1000.msh ]; then
            echo "RefSeq MASH sketch does not exist, initializing..."
            wget --no-verbose https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
        fi

        touch {output}
        
        """

        
# --------------------------
#    Initialize MacSyFinder
# --------------------------

rule txsscan_init:
    output:
        touch(os.path.join(output_dir, "data", "txsscan", ".txsscan.setup.done"))
    conda:
        "../envs/macsyfinder.yaml"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_txsscan.txt")
    shell:
        """
        MSF_DB="dbs/macysyfinder/models"
        
        if [ ! -d "$MSF_DB" ]; then
            echo "TXSScan models do not exist, initializing..."
            msf_data install TXSScan --target $MSF_DB > /dev/null 2>&1
        fi
        
        touch {output}
        
        """
        
        