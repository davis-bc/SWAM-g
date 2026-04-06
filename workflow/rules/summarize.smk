# ---------------------------
#     Generate Mashtree
# ---------------------------

rule mashtree:
    input:
        symlink = expand(os.path.join(output_dir, "data", "unicycler", "batch", "{sample}.fasta"), sample=samples)
    output:
        tree = os.path.join(output_dir, "mashtree.nwk")
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "mashtree.txt")
    threads: 8
    conda:
        "../envs/mashtree.yaml"
    shell:
        r"""
        inputs=({input.symlink})

        if [ "${{#inputs[@]}}" -eq 1 ]; then
            sample_name=$(basename "${{inputs[0]}}" .fasta)
            printf '%s;\n' "$sample_name" > {output.tree}
        else
            mashtree --numcpus {threads} "${{inputs[@]}}" > {output.tree}
        fi
        """

# ---------------------------
#    Generate Summary File
# ---------------------------

rule pd_isolate_metadata:
    input:
        mash = os.path.join(output_dir, "data", "mash", "mash_taxonomy.tsv")
    output:
        metadata    = os.path.join(output_dir, "data", "pd", "pd_isolate_metadata.tsv"),
        comparators = os.path.join(output_dir, "data", "pd", "pd_cluster_comparators.tsv")
    params:
        samples             = samples,
        enabled             = pd_lookup_enabled,
        backend             = pd_lookup_backend,
        comparator_limit    = pd_comparator_limit,
        sample_metadata_tsv = pd_sample_metadata_tsv,
        pd_isolates_tsv     = pd_isolates_tsv,
        pd_exceptions_tsv   = pd_exceptions_tsv
    conda:
        "../envs/pd_lookup.yaml"
    script:
        "../scripts/pd_isolate_lookup.py"


rule summarize_results:
    input:
        checkm               = expand(os.path.join(output_dir, "data", "checkm2", "{sample}", "quality_report.tsv"), sample=samples),
        mash                 = os.path.join(output_dir, "data", "mash", "mash_taxonomy.tsv"),
        mlst                 = os.path.join(output_dir, "data", "mlst", "mlst.tsv"),
        pd_metadata          = os.path.join(output_dir, "data", "pd", "pd_isolate_metadata.tsv"),
        pd_comparators       = os.path.join(output_dir, "data", "pd", "pd_cluster_comparators.tsv"),
        ectyper_files        = expand(os.path.join(output_dir, "data", "serotype", "E.coli", "{sample}", "output.tsv"), sample=samples),
        resfinder_files      = expand(os.path.join(output_dir, "data", "resfinder", "{sample}", "ResFinder_results_tab.txt"), sample=samples),
        pf_files             = expand(os.path.join(output_dir, "data", "resfinder", "{sample}", "PointFinder_results.txt"), sample=samples),
        afp_files            = expand(os.path.join(output_dir, "data", "amrfinderplus", "{sample}.afp.tsv"), sample=samples),
        seqsero_files        = expand(os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "SeqSero_result.tsv"), sample=samples),
        sistr_files          = expand(os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "sistr.tsv"), sample=samples),
        mobtyper_files       = expand(os.path.join(output_dir, "data", "mob-suite", "{sample}", "mobtyper_results.txt"), sample=samples),
        coverage_files       = expand(os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}_coverage.tsv"), sample=samples),
        mef_files            = expand(os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.csv"), sample=samples),
        contig_files         = expand(os.path.join(output_dir, "data", "mob-suite", "{sample}", "contig_report.txt"), sample=samples),
        mobtyper_blast_files = expand(os.path.join(output_dir, "data", "mob-suite", "{sample}", "biomarkers.blast.txt"), sample=samples),
        txsscan_files        = expand(os.path.join(output_dir, "data", "txsscan", "{sample}", "all_systems.tsv"), sample=samples),
        txsscan_files2       = expand(os.path.join(output_dir, "data", "txsscan", "{sample}", "all_systems.txt"), sample=samples),
        prodigal_files       = expand(os.path.join(output_dir, "data", "unicycler", "{sample}", "{sample}.prot.faa"), sample=samples)
    output:
        xlsx = os.path.join(output_dir, "SWAM-g_results.xlsx"),
        csv  = os.path.join(output_dir, "contig_map.csv")
    params:
        debug = debug_mode
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/data_summary.R"
