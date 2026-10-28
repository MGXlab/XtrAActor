// =============================================================
// 9. plot_zvalue
// =============================================================
process PLOT_ZVALUE {
    tag "${cluster_name}"
    cpus 1
    memory '10 GB'
    time '1h'

    publishDir "${params.out_dir}/test_phyloglm", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    tuple val(cluster_name), path(csv)

    output:
    path "${cluster_name}_zvalue_plot.png", emit: zvalue_plot
    path "${cluster_name}_zvalue.log", emit: log

    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate phyloglm

    Rscript ${params.shared_bin}/plot_zvalue_featimps.R \\
        ${csv} ${cluster_name}_zvalue_plot.png > ${cluster_name}_zvalue.log 2>&1
    """
}

// =============================================================
// 10. Plots (global, across all clusters)
// =============================================================
process PLOTS {
    cpus 4
    memory { 30.GB * task.attempt }
    time   { 2.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/plots", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    path(rf_files,    stageAs: 'RF_results/raw/*')   // staged from task outputs, not from publishDir
    path(fasta_files, stageAs: 'cluster_fastas/*')
    val biome_1
    val biome_2
    val colour1
    val colour2

    output:
    path "*_seqLogo.png",                         emit: seq_logo_plot,                 optional: true
    path "*_biome_feat_enrich_bars.png",          emit: feature_enrichment_plot,       optional: true
    path "*_biome_feat_enrich_bars_fullprot.png", emit: feature_full_protein_plot,     optional: true
    path "*_AA_radars.png",                       emit: AA_radars_plot,                optional: true
    path "*_AA_enrich_bars.png",                  emit: protein_analysis_plot,         optional: true
    path "*_protein_analysis_subplots.png",       emit: protein_analysis_subplots_png, optional: true
    path "*_protein_analysis_subplots.svg",       emit: protein_analysis_subplots_svg, optional: true
    path "*_protein_scale_flexibility.png",       emit: protein_scale_flexibility,     optional: true
    path "*_protein_scale_GRAVY.png",             emit: protein_scale_GRAVY,           optional: true
    path "*_protein_scale_chargeatpH7.png",       emit: protein_scale_chargeatpH7,     optional: true
    path "Rf_model_report_plots.png",             emit: report_plot_png
    path "Rf_model_report_plots.svg",             emit: report_plot_svg
    path "Rf_model_precision_v_recall.png",       emit: precision_recall
    path "plots.log",                             emit: log

    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate working-env

    results_path=\$(realpath RF_results/raw)
    fasta_path=\$(realpath cluster_fastas)

    python3 ${params.shared_bin}/plots/plotting_script.py \\
        "\$results_path" "\$fasta_path" "${biome_1}" "${biome_2}" "${colour1}" "${colour2}" 50 \$PWD \\
        > plots.log 2>&1
    """
}