# ---------------------------
#    Generate Summary File
# ---------------------------

rule summarize_results:
    input:
        checkm               = os.path.join(output_dir, "data", "checkm", "genome.stats.tsv"),
        checkm_stats         = os.path.join(output_dir, "data", "checkm", "storage", "bin_stats.analyze.tsv"),
        gtdbtk               = os.path.join(output_dir, "data", "gtdb-tk", "gtdbtk.bac120.summary.tsv"),
        mlst                 = os.path.join(output_dir, "data", "mlst", "mlst.tsv"),
        resfinder_files      = expand(os.path.join(output_dir, "data", "resfinder", "{sample}", "ResFinder_results_tab.txt"), sample=samples),
        pf_files             = expand(os.path.join(output_dir, "data", "resfinder", "{sample}", "PointFinder_results.txt"), sample=samples),
        amr_cr_files         = expand(os.path.join(output_dir, "data", "amrfinderplus", "chromosomes", "{sample}.afp.tsv"), sample=samples),
        amr_plas_done        = expand(os.path.join(output_dir, "data", "amrfinderplus", "plasmids", ".{sample}.plasmids.afp.done"), sample=samples),
        seqsero_files        = expand(os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "SeqSero_result.tsv"), sample=samples),
        mobtyper_files       = expand(os.path.join(output_dir, "data", "mob-suite", "{sample}", "mobtyper_results.txt"), sample=samples),
        coverage_files       = expand(os.path.join(output_dir, "data", "bams", "{sample}_coverage.tsv"), sample=samples),
        ectyper_files        = expand(os.path.join(output_dir, "data", "serotype", "E.coli", "{sample}", "output.tsv"), sample=samples),
        mef_files            = expand(os.path.join(output_dir, "data", "mobileelementfinder", "{sample}", "{sample}.csv"), sample=samples),
        contig_files         = expand(os.path.join(output_dir, "data", "mob-suite", "{sample}", "contig_report.txt"), sample=samples),
        mobtyper_blast_files = expand(os.path.join(output_dir, "data", "mob-suite", "{sample}", "biomarkers.blast.txt"), sample=samples),
        txsscan_files        = expand(os.path.join(output_dir, "data", "txsscan", "{sample}", "all_systems.tsv"), sample=samples),
        txsscan_files2       = expand(os.path.join(output_dir, "data", "txsscan", "{sample}", "all_systems.txt"), sample=samples)
    output:
        xlsx = os.path.join(output_dir, "SWAM-g_results.xlsx"),
        csv  = os.path.join(output_dir, "contig_map.csv")
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/data_summary.R"
