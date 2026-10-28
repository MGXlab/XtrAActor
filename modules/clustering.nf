// =============================================================
// 0b. MMseqs2 clustering (all genomes pooled)
// =============================================================
process MMSEQ2 {
    cpus 10
    memory { 50.GB * task.attempt }
    time   { 6.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/mmseq2", mode: 'copy', pattern: "mmseq2*"
    queue 'standard,short,long,fat'

    input:
    path(prediction_file)
    val(min_id)
    val(min_cov)

    output:
    path "mmseq2_all_seqs.fasta", emit: all_fasta
    path "mmseq2_cluster.tsv",    emit: tsv
    path "mmseq2_rep_seq.fasta",  emit: rep_fasta
    path "mmseq2.log",            emit: log

    script:
    """
    source ${params.conda_sh_veo}
    conda activate mmseq2_v14.7e284
    mmseqs easy-cluster ${prediction_file} mmseq2 tmp \\
        --min-seq-id ${min_id} -c ${min_cov} --threads ${task.cpus} > mmseq2.log 2>&1
    """
}

// =============================================================
// 0c. Orthogroups: per-cluster fasta + per-cluster metadata
// =============================================================
process ORTHOGROUPS {
    cpus 2
    memory { 20.GB * task.attempt }
    time   { 2.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    // orthogroups.py writes into the task root -> sort into sub-folders on publish
    publishDir "${params.out_dir}", mode: 'copy',
        saveAs: { fn ->
            if (fn.endsWith('.fasta'))      "cluster_fastas/${fn}"
            else if (fn.endsWith('.tsv'))   "cluster_biomes/${fn}"
            else                            "logs/${fn}"
        }
    queue 'standard,short,long,fat'

    input:
    path cluster_tsv
    path proteins_fasta
    path metadata_file

    output:
    path "cluster_fastas/*.fasta",       emit: cluster_fastas
    path "cluster_biomes/*_biome.tsv",   emit: cluster_biomes
    path "cluster_size_histogram.png",   emit: size_hist
    path "orthogroups.log",              emit: log

    script:
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate working-env

    python3 ${projectDir}/bin/orthogroups.py \\
        ${cluster_tsv} ${proteins_fasta} ${metadata_file} \\
        --fasta_out_dir cluster_fastas --biome_out_dir cluster_biomes 2> orthogroups.log
    """
}