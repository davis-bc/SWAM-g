# ----------------------------------------
# ----------------------------------------
#               Setup checkM2
# ---------------------------------------

rule checkm_init:
    output:
        touch(os.path.join(output_dir, "data", "checkm2", ".checkm_initialized"))
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_checkm.txt")
    conda: 
        "../envs/checkm2.yaml"
    shell:
        """
        DB_DIR="dbs"
        DB="dbs/CheckM2_database/uniref100.KO.1.dmnd"

        mkdir -p "$DB_DIR"
        
        if [ ! -f "$DB" ]; then
            echo "CheckM2_database not found. Downloading database..."
            checkm2 database --download --path "$DB_DIR" > /dev/null 2>&1
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
    params:
        db_dir = "dbs/mob_suite"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_mobsuite.txt")
    shell:
        """
        MOB_DB="{params.db_dir}"

        mkdir -p "$MOB_DB"

        # Pre-stage the MOB-suite database in a persistent repo-level directory so
        # compute-node jobs never fall back to the package default and attempt an
        # automatic download without internet access.
        if [ -z "$(find "$MOB_DB" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
            echo "MOB-suite databases not found in $MOB_DB, initializing..."
            mob_init --database_directory "$MOB_DB"
        fi

        if [ -z "$(find "$MOB_DB" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
            echo "MOB-suite initialization did not populate $MOB_DB" >&2
            exit 1
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

        mkdir -p "$DB_DIR"

        # Check and clone only if the database directories are missing
        [ ! -d "$DB_DIR/resfinder_db" ] && git clone https://bitbucket.org/genomicepidemiology/resfinder_db/ "$DB_DIR/resfinder_db" > /dev/null 2>&1
        [ ! -d "$DB_DIR/pointfinder_db" ] && git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/ "$DB_DIR/pointfinder_db" > /dev/null 2>&1
        [ ! -d "$DB_DIR/disinfinder_db" ] && git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/ "$DB_DIR/disinfinder_db" > /dev/null 2>&1

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
        MASH_SKETCH="$DB_DIR/EnteroRef_GTDBSketch_20231003_V2.msh"
        URL="https://zenodo.org/records/13969103/files/EnteroRef_GTDBSketch_20231003_V2.msh"

        mkdir -p "$DB_DIR"

        # Check if the MASH sketch exists, download if missing
        if [ ! -f "$MASH_SKETCH" ]; then
            echo "EnteroRef_GTDB MASH sketch does not exist, initializing..."
            wget -O "$MASH_SKETCH" "$URL"
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
        MSF_ROOT="dbs/macsyfinder"
        MSF_DB="$MSF_ROOT/models"
        TXSSCAN_DIR="$MSF_DB/TXSScan"
        TXSSCAN_META="$TXSSCAN_DIR/metadata.yml"

        mkdir -p "$MSF_DB"

        if command -v msf_data >/dev/null 2>&1; then
            MSF_DATA_BIN=msf_data
        elif command -v macsydata >/dev/null 2>&1; then
            MSF_DATA_BIN=macsydata
        else
            echo "Could not find msf_data or macsydata in the MacSyFinder environment" >&2
            exit 1
        fi

        needs_install=0
        if [ ! -d "$TXSSCAN_DIR" ]; then
            needs_install=1
        elif [ ! -f "$TXSSCAN_META" ]; then
            needs_install=1
        elif ! find "$TXSSCAN_DIR" -type f \\( -name "*.xml" -o -name "*.hmm" \\) -print -quit | grep -q .; then
            needs_install=1
        fi

        if [ "$needs_install" -eq 1 ]; then
            echo "Installing or refreshing TXSScan models in $MSF_DB..."
            rm -rf "$TXSSCAN_DIR"
            "$MSF_DATA_BIN" install --target "$MSF_DB" TXSScan
        fi

        if [ ! -f "$TXSSCAN_META" ]; then
            echo "TXSScan metadata file missing after install: $TXSSCAN_META" >&2
            exit 1
        fi

        if ! find "$TXSSCAN_DIR" -type f \\( -name "*.xml" -o -name "*.hmm" \\) -print -quit | grep -q .; then
            echo "TXSScan install at $TXSSCAN_DIR is incomplete" >&2
            exit 1
        fi
        
        touch {output}
        
        """
        
        
