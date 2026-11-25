// =============================================================
// 1. MSA (mafft)
// =============================================================
process MAFFT {
    tag "${cluster_name}"
    cpus 10
    memory { 50.GB * task.attempt }
    time   { 6.h * task.attempt }
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/mafft", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    tuple val(cluster_name), path(fasta_file)

    output:
    tuple val(cluster_name), path("${cluster_name}_mafft_aa.fa"), emit: aligned
    path "${cluster_name}_mafft.log", emit: log

    script:
    """
    # Split fasta into ref.fa (longest sequence) and remain.fa (everything else)
    awk -v reff="ref.fa" -v restf="remain.fa" '
    /^>/ {
        if (seq != "") { n++; headers[n] = header; seqs[n] = seq }
        header = \$0
        sub(/\\r\$/, "", header)
        seq = ""
        next
    }
    {
        line = \$0
        gsub(/\\r/, "", line)
        gsub(/[ \\t]+\$/, "", line)
        seq = seq line
    }
    END {
        if (seq != "") { n++; headers[n] = header; seqs[n] = seq }
        maxlen = 0
        for (i = 1; i <= n; i++) {
            if (length(seqs[i]) > maxlen) { maxlen = length(seqs[i]); best = i }
        }
        print headers[best] "\\n" seqs[best] > reff
        for (i = 1; i <= n; i++) {
            if (i != best) print headers[i] "\\n" seqs[i] > restf
        }
    }
    ' ${fasta_file}

    ${params.mafft_bin} --thread ${task.cpus} --keeplength --add remain.fa ref.fa \\
        > ${cluster_name}_mafft_aa.fa 2> ${cluster_name}_mafft.log
    """
}

// =============================================================
// 2. Optional trimAl
// =============================================================
process TRIMAL {
    tag "${cluster_name}"
    cpus 1
    memory { 10.GB * task.attempt }
    time   { 1.h * task.attempt }

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/mafft", mode: 'copy'
    queue 'standard,short,long,fat,gpu-veo'

    input:
    tuple val(cluster_name), path(msa_file)

    output:
    tuple val(cluster_name), path("${cluster_name}_mafft_aa_trimmed.fa"), emit: faa
    path "${cluster_name}_trimal.log", emit: log

    script:
    """
    ${params.trimal_bin} -in ${msa_file} \\
        -out ${cluster_name}_mafft_aa_trimmed.fa \\
        -gappyout \\
        -colnumbering > ${cluster_name}_trimal.log 2>&1
    """
}

// =============================================================
// 3. FastTree
// =============================================================
process FASTTREE {
    tag "${cluster_name}"
    cpus 10
    memory { 50.GB * task.attempt }
    time   { 6.h * task.attempt }
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/fasttree", mode: 'copy'
    queue 'standard,short,long,fat'

    input:
    tuple val(cluster_name), path(msa)

    output:
    tuple val(cluster_name), path("${cluster_name}_combined_fastree_lg_pseudo.tree"), emit: tree
    path "${cluster_name}_fasttree.log", emit: log

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    ${params.fasttree_bin} -lg -pseudo < ${msa} \\
        > ${cluster_name}_combined_fastree_lg_pseudo.tree 2> ${cluster_name}_fasttree.log
    """
}
