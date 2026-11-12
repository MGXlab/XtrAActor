// =============================================================
// Single-biome cluster analysis
//   all-v-all DIAMOND blastp -> Dice distances + genome tree (dice_tree.py)
//   -> phyloglm of every single-biome cluster against the genome tree
//   -> report of the clusters with padj < params.sb_padj
// =============================================================

// -------------------------------------------------------------
// S1. DIAMOND database of all predicted proteins
// -------------------------------------------------------------
process DIAMOND_MAKEDB {
    cpus 8
    memory { 20.GB * task.attempt }
    time   { 2.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2
    queue 'standard,short,long,fat'

    input:
    path(all_genes)

    output:
    path "all_genes_diamond_db.dmnd", emit: db
    path "diamond_makedb.log",        emit: log

    script:
    def activate = params.diamond_env ? "source ${params.conda_sh_veo}\n    conda activate ${params.diamond_env}" : ''
    """
    ${activate}

    ${params.diamond_bin} makedb --in ${all_genes} -d ./all_genes_diamond_db \\
        --threads ${task.cpus} > diamond_makedb.log 2>&1
    """
}

// -------------------------------------------------------------
// S2. all-v-all blastp, one task per query chunk
// -------------------------------------------------------------
process DIAMOND_BLASTP {
    tag "${query_chunk.baseName}"
    cpus 20
    memory { 60.GB * task.attempt }
    time   { 12.h * task.attempt }
    maxForks 20

    errorStrategy 'retry'
    maxRetries 2
    queue 'standard,short,long,fat'

    input:
    path(query_chunk)
    path(db)

    output:
    path "${query_chunk.baseName}_diamond.out", emit: hits
    path "${query_chunk.baseName}_diamond.log", emit: log

    script:
    def activate = params.diamond_env ? "source ${params.conda_sh_veo}\n    conda activate ${params.diamond_env}" : ''
    """
    ${activate}

    ${params.diamond_bin} blastp -d ${db} -q ${query_chunk} -o ${query_chunk.baseName}_diamond.out \\
        -k 0 --evalue ${params.diamond_evalue} --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \\
        ${params.diamond_args} > ${query_chunk.baseName}_diamond.log 2>&1
    """
}

// -------------------------------------------------------------
// S3. Dice distances + neighbour-joining genome tree
// -------------------------------------------------------------
process DICE_TREE {
    cpus 2
    memory { 100.GB * task.attempt }
    time   { 24.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/single_biome/tree", mode: 'copy'
    queue 'standard,long,fat'

    input:
    path(hits, stageAs: 'hits/*')

    output:
    path "dice_distances.tsv", emit: distances
    path "dice_tree.nwk",      emit: tree
    path "dice_tree.log",      emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    cat hits/*_diamond.out > all_v_all_diamond.out

    python3 ${projectDir}/bin/dice_tree.py -i all_v_all_diamond.out \\
        -out_csv dice_distances.tsv -out_tree dice_tree.nwk \\
        -e ${params.dice_evalue} -p ${params.dice_pident} ${params.dice_args} \\
        > dice_tree.log 2>&1
    """
}

// -------------------------------------------------------------
// S4. Clusters of the 50/50 clustering with members from one biome only
// -------------------------------------------------------------
process SINGLE_BIOME_CLUSTERS {
    cpus 1
    memory { 30.GB * task.attempt }
    time   { 2.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/single_biome", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }
    queue 'standard,short,long,fat'

    input:
    path(cluster_tsv)
    path(derep_tsvs, stageAs: 'derep/*')    // empty list if --dereplicate false
    path(all_genes)
    path(metadata_file)

    output:
    path "single_biome_clusters.tsv", emit: clusters
    path "single_biome_clusters.log", emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/single_biome_clusters.py \\
        ${cluster_tsv} ${metadata_file} --derep_dir derep --proteins_fasta ${all_genes} \\
        --min_genomes ${params.sb_min_genomes} \\
        -o single_biome_clusters.tsv 2> single_biome_clusters.log
    """
}

// -------------------------------------------------------------
// S5. phyloglm per single-biome cluster, one task per chunk of clusters
// -------------------------------------------------------------
process PHYLOGLM_SINGLE_BIOME {
    tag "${cluster_chunk.baseName}"
    cpus 10
    memory { 40.GB * task.attempt }
    time   { 4.h * task.attempt }
    maxForks 20

    errorStrategy 'retry'
    maxRetries 2
    queue 'standard,short,long,fat'

    input:
    path(cluster_chunk)
    path(tree)
    path(metadata_file)

    output:
    path "${cluster_chunk.baseName}_partial.tsv", emit: partial
    path "${cluster_chunk.baseName}_fit.log",     emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate phyloglm

    Rscript ${projectDir}/bin/phyloglm_single_biome.R \\
        ${cluster_chunk} ${tree} ${metadata_file} "${params.biome2}" \\
        ${cluster_chunk.baseName}_partial.tsv > ${cluster_chunk.baseName}_fit.log 2>&1
    """
}

// -------------------------------------------------------------
// S6. merge chunks, FDR once over all clusters, report padj < threshold
// -------------------------------------------------------------
process SINGLE_BIOME_REPORT {
    cpus 1
    memory { 10.GB * task.attempt }
    time   { 1.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/single_biome", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }
    queue 'standard,short,long,fat'

    input:
    path(partials, stageAs: 'partials/*')

    output:
    path "single_biome_phyloglm_all.tsv",    emit: all
    path "single_biome_cluster_report.tsv",  emit: report
    path "single_biome_report.log",          emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate phyloglm

    Rscript ${projectDir}/bin/phyloglm_single_biome_report.R \\
        partials ${params.sb_padj} \\
        single_biome_phyloglm_all.tsv single_biome_cluster_report.tsv \\
        > single_biome_report.log 2>&1
    """
}
