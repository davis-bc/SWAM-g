# ---------------------------
#    Generate Summary File
# ---------------------------

rule summarize_results:
    input:
        checkm = os.path.join(output_dir, "data", "checkm", "genome.stats.tsv"),
        checkm_stats = os.path.join(output_dir, "data", "checkm", "storage", "bin_stats.analyze.tsv"),
        gtdbtk = os.path.join(output_dir, "data", "gtdb-tk", "gtdbtk.bac120.summary.tsv"),
        mlst   = os.path.join(output_dir, "data", "mlst", "mlst.tsv"),
        resfinder_files = expand(os.path.join(output_dir, "data", "resfinder", "{sample}", "ResFinder_results_tab.txt"), sample=samples),
        amr_cr_files = expand(os.path.join(output_dir, "data", "amrfinderplus", "chromosomes", "{sample}.afp.tsv"), sample=samples),
        amr_plas_done = expand(os.path.join(output_dir, "data", "amrfinderplus", "plasmids", ".{sample}.plasmids.afp.done"), sample=samples),
        seqsero_files = expand(os.path.join(output_dir, "data", "serotype", "Salmonella", "{sample}", "SeqSero_result.tsv"), sample=samples),
        mobtyper_files = expand(os.path.join(output_dir, "data", "mob-suite", "{sample}", "mobtyper_results.txt"), sample=samples),
        coverage_files = expand(os.path.join(output_dir, "data", "bams", "{sample}_coverage.tsv"), sample=samples),
        ectyper_files = expand(os.path.join(output_dir, "data", "serotype", "E.coli", "{sample}", "output.tsv"), sample=samples)
    output:
        summary = os.path.join(output_dir, "results_summary.csv"),
        seqs2 = os.path.join(output_dir, "salmonella_serotype_all.csv"),
        plasmid_summary = os.path.join(output_dir, "plasmid_data_long.csv"),
        assembly = os.path.join(output_dir, "QA_assembly_stats.csv"),
        ect = os.path.join(output_dir, "ecoli_serotype_all.csv")
    conda:
        "../envs/Renv.yaml"
    script:
        "../scripts/data_summary.R"
