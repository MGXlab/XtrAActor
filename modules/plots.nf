// =============================================================
// 10a. RF report plots across all clusters (waits for all RF models)
// =============================================================
process PLOT_RF_REPORT {
    cpus 1
    memory { 16.GB * task.attempt }
    time   { 1.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/plots/overall_cluster_plots", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    path(reports, stageAs: 'RF_results/raw/*')
    val biome_1
    val biome_2
    val colour1
    val colour2

    output:
    path "Rf_model_report_plots.png",       emit: report_plot_png
    path "Rf_model_report_plots.svg",       emit: report_plot_svg
    path "Rf_model_precision_v_recall.png", emit: precision_recall
    path "plot_rf_report.log",              emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    results_path=\$(realpath RF_results/raw)

    python3 ${params.shared_bin}/plots/plotting_script.py \\
        "\$results_path" NONE "${biome_1}" "${biome_2}" "${colour1}" "${colour2}" 50 \$PWD report \\
        > plot_rf_report.log 2>&1
    """
}

// =============================================================
// 10b. Plots per cluster (each starts as soon as its own RF model is done)
// =============================================================
process PLOT_CLUSTER {
    tag "${cluster_name}"
    cpus 1
    memory { 16.GB * task.attempt }
    time   { 1.h * task.attempt }
    maxForks 20

    errorStrategy 'retry'
    maxRetries 2

    // plots/individual_cluster_plots/<cluster>/ (plots and the log of that cluster)
    publishDir "${params.out_dir}/plots/individual_cluster_plots/${cluster_name}", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    tuple val(cluster_name), path(rf_files, stageAs: 'RF_results/raw/*'), path(fasta, stageAs: 'cluster_fastas/*')
    val biome_1
    val biome_2
    val colour1
    val colour2

    output:
    path "${cluster_name}_*.png",       emit: png, optional: true
    path "${cluster_name}_*.svg",       emit: svg, optional: true
    path "${cluster_name}_plots.log",   emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    results_path=\$(realpath RF_results/raw)
    fasta_path=\$(realpath cluster_fastas)

    python3 ${params.shared_bin}/plots/plotting_script.py \\
        "\$results_path" "\$fasta_path" "${biome_1}" "${biome_2}" "${colour1}" "${colour2}" 50 \$PWD ${cluster_name} \\
        > ${cluster_name}_plots.log 2>&1
    """
}

// =============================================================
// 11. Protein properties compared across all clusters (one volcano plot per property)
// =============================================================
process CLUSTER_PROPERTY_PLOTS {
    cpus 1
    memory { 16.GB * task.attempt }
    time   { 2.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/plots/overall_cluster_plots", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    path(data_dfs,    stageAs: 'RF_results/raw/*')   // <cluster>_data_df.tsv
    path(fasta_files, stageAs: 'cluster_fastas/*')
    val biome_1
    val biome_2
    val colour1
    val colour2

    output:
    path "volcano_*.png",                        emit: volcano_plots
    path "protein_properties_per_protein.tsv",   emit: per_protein
    path "protein_properties_cluster_stats.tsv", emit: cluster_stats
    path "cluster_property_plots.log",           emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/plots/plots/protein_property_volcano.py RF_results/raw cluster_fastas \\
        --biome1 "${biome_1}" --biome2 "${biome_2}" \\
        --colour1 "${colour1}" --colour2 "${colour2}" \\
        > cluster_property_plots.log 2>&1
    """
}
