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
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "init_mobsuite.txt")
    shell:
        """
        # Resolve the databases path within the currently active conda env.
        # The old find-based check matched databases from *any* env (including stale
        # ones), so if the mob_suite env was recreated with a new hash the databases
        # would appear "found" while the new env had none — causing mob_recon to
        # attempt a download on compute nodes (which have no internet access).
        MOB_DB=$(python -c "import mob_suite; import os; print(os.path.join(os.path.dirname(mob_suite.__file__), 'databases'))")
        if [ ! -d "$MOB_DB" ] || [ -z "$(ls -A "$MOB_DB" 2>/dev/null)" ]; then
            echo "MOB-suite databases not found in current conda env, initializing..."
            mob_init
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

        mkdir -p "$MSF_ROOT"
        
        if [ ! -d "$MSF_DB/TXSScan" ]; then
            echo "TXSScan models not found, installing..."
            macsydata install --target "$MSF_DB" TXSScan
        else
            # MacSyFinder 2.x requires model XML version '2.0'. TXSScan-1.x models
            # (installed by older runs) are incompatible. Detect this by checking the
            # metadata.yml that MacSyFinder 2.x model packages include at their root.
            meta="$MSF_DB/TXSScan/metadata.yml"
            if ! grep -q "^vers: [2-9][.]" "$meta" 2>/dev/null; then
                echo "TXSScan v1 models detected (incompatible with MacSyFinder 2.x), reinstalling..."
                rm -rf "$MSF_DB"
                mkdir -p "$MSF_ROOT"
                macsydata install --target "$MSF_DB" TXSScan
            fi
        fi
        
        touch {output}
        
        """
        
        
