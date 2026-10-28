// =============================================================
// 5. Create phyloglm dataset
// =============================================================
process CREATE_PHYLOGLM_DATASET {
    tag "${cluster_name}"
    cpus 2
    memory { 30.GB * task.attempt }
    time   { 1.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/test_phyloglm", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    tuple val(cluster_name), path(feature), path(data)

    output:
    tuple val(cluster_name), path("${cluster_name}_phyloglm_dataset.tsv"), emit: dataset
    // NOTE: assumes create_dataset_phyloglm.py also writes <cluster>_gini.tsv
    tuple val(cluster_name), path("${cluster_name}_gini.tsv"),             emit: gini
    path "${cluster_name}_create_dataset.log", emit: log

    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate working-env

        python3 ${params.shared_bin}/create_dataset_phyloglm.py \\
        ${feature} ${data} -o ${cluster_name}_phyloglm_dataset.tsv \\
        --target_biome "${params.biome2}" \\
        > ${cluster_name}_create_dataset.log 2>&1
    """
}

// =============================================================
// 6. Split gini into chunks (never produces empty chunks)
// =============================================================
process SPLIT_GINI {
    tag "${cluster_name}"
    cpus 1
    memory '2 GB'
    time '30m'

    queue 'standard,short,long,fat'

    input:
    tuple val(cluster_name), path(gini_file)

    output:
    tuple val(cluster_name), path("${cluster_name}_chunk_*.tsv"), emit: chunks, optional: true

    script:
    """
    header=\$(head -n 1 ${gini_file})
    tail -n +2 ${gini_file} | grep -v '^[[:space:]]*\$' > body.tsv || true
    n_lines=\$(wc -l < body.tsv)

    if [ "\$n_lines" -eq 0 ]; then
        echo "No features in ${gini_file}, skipping ${cluster_name}" >&2
        exit 0
    fi

    n_chunks=\$(( n_lines < ${params.n_chunks} ? n_lines : ${params.n_chunks} ))
    per_chunk=\$(( (n_lines + n_chunks - 1) / n_chunks ))
    split -d -a 3 -l \$per_chunk body.tsv tmp_chunk_

    for f in tmp_chunk_*; do
        idx=\${f#tmp_chunk_}
        outname="${cluster_name}_chunk_\${idx}.tsv"
        echo "\$header" > "\$outname"
        cat "\$f" >> "\$outname"
    done
    """
}

// =============================================================
// 7. phyloglm_fit_chunk
// =============================================================
process PHYLOGLM_FIT_CHUNK {
    tag "${gini_chunk.baseName}"
    cpus 20
    memory { 80.GB * task.attempt }
    time   { 4.h * task.attempt }
    maxForks 20

    errorStrategy 'retry'
    maxRetries 2

    queue 'standard,short,long,fat'

    input:
    tuple val(cluster_name), path(dataset), path(tree), path(gini_chunk)

    output:
    tuple val(cluster_name), path("${gini_chunk.baseName}_partial.csv"), emit: partial
    path "${gini_chunk.baseName}_fit.log", emit: log

    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate phyloglm

    Rscript ${params.shared_bin}/phyloglm_fit.R \\
        ${dataset} ${tree} ${gini_chunk} ${gini_chunk.baseName}_partial.csv \\
        > ${gini_chunk.baseName}_fit.log 2>&1
    """
}

// =============================================================
// 8. merge_phyloglm
// =============================================================
process MERGE_PHYLOGLM {
    tag "${cluster_name}"
    cpus 1
    memory { 20.GB * task.attempt }
    time   { 1.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/test_phyloglm", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    tuple val(cluster_name), path(partials)

    output:
    tuple val(cluster_name), path("${cluster_name}_phyloglm.csv"), emit: csv
    path "${cluster_name}_phyloglm_plot.png", emit: plot
    path "${cluster_name}_merge.log", emit: log

    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate phyloglm

    Rscript ${params.shared_bin}/phyloglm_merge.R \\
        ${partials} ${cluster_name}_phyloglm.csv ${cluster_name}_phyloglm_plot.png \\
        > ${cluster_name}_merge.log 2>&1
    """
}

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
