// =============================================================
// 4. RF_MODEL
// =============================================================
process RF_MODEL {
    tag "${cluster_name}"
    cpus 30
    memory { 100.GB * task.attempt }
    time   { 6.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    // flatten RF_results/ so files land directly in RF_results/raw/
    publishDir "${params.out_dir}/RF_results/raw", mode: 'copy',
        saveAs: { fn -> file(fn).name }
    queue 'standard,short,long,fat,gpu-veo'

    input:
    tuple val(cluster_name), path(msa_file), path(metadata_file)

    output:
    tuple val(cluster_name), path("RF_results/*_Feature_importance.tsv"),  emit: features
    tuple val(cluster_name), path("RF_results/*_RandomForest_Report.tsv"), emit: reports
    tuple val(cluster_name), path("RF_results/*_data_df.tsv"),             emit: data_df
    path "${cluster_name}_rf_model.log", emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/rf_model_v3.py \\
        ${msa_file} ${metadata_file} ${cluster_name} > ${cluster_name}_rf_model.log 2>&1
    """
}