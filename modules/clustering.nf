// =============================================================
// 0b. Dereplication: remove near-identical proteins WITHIN each biome
// =============================================================
process MMSEQ2_DEREP {
    tag "${biome_set}"
    publishDir "${params.out_dir}/mmseq2_derep", mode: 'copy'

    input:
    tuple val(biome_set), path(proteins)
    val(min_id)
    val(min_cov)

    output:
    tuple val(biome_set), path("${biome_set}_derep_rep_seq.fasta"), emit: rep_fasta
    path "${biome_set}_derep_cluster.tsv",                           emit: tsv
    path "${biome_set}_derep.log",                                   emit: log

    script:
    """
    source ${params.conda_sh_veo}
    conda activate mmseq2_v14.7e284
    mmseqs easy-cluster ${proteins} ${biome_set}_derep tmp \\
        --min-seq-id ${min_id} -c ${min_cov} --threads ${task.cpus} > ${biome_set}_derep.log 2>&1
    """
}

// =============================================================
// 0c. MMseqs2 clustering (dereplicated proteins, all biomes pooled)
// =============================================================
process MMSEQ2 {
    publishDir "${params.out_dir}/mmseq2", mode: 'copy', pattern: "mmseq2*"

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
// 0d. Orthogroups: per-cluster fasta + per-cluster biome table
// =============================================================
process ORTHOGROUPS {
    // outputs are already in sub-folders -> keep paths; logs go to logs/
    publishDir "${params.out_dir}", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }

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
    source ${params.conda_sh_ri}
    conda activate working-env
    # random comment line

    python3 ${projectDir}/bin/orthogroups_sel_by_biome.py \\
        ${cluster_tsv} ${proteins_fasta} ${metadata_file} \\
        --fasta_out_dir cluster_fastas --biome_out_dir cluster_biomes 2> orthogroups.log
    """
}