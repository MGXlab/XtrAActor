// Both processes emit proteins renamed to <genome_id>_<n>, so the downstream
// protein IDs look the same whichever tool produced them.
// TODO Remember to adjust pharokka_args to default to --fast for large number of sequences

process PHAROKKA {
    tag "${genome_id}"
    publishDir "${params.out_dir}/annotation/pharokka", mode: 'copy'

    input:
    tuple val(genome_id), path(genome)

    output:
    tuple val(genome_id), path("${genome_id}_pharokka"),      emit: annotation
    tuple val(genome_id), path("${genome_id}_proteins.faa"),  emit: proteins
    path("${genome_id}_pharokka.log"),                        emit: log

    script:
    def db = params.pharokka_db ? "-d ${params.pharokka_db}" : ''
    """
    source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
    conda activate ${params.pharokka_env}

    pharokka run -i ${genome} -o ${genome_id}_pharokka ${db} \\
        -t ${task.cpus} -g ${params.gene_predictor} -p ${genome_id} \\
        ${params.pharokka_args} > ${genome_id}_pharokka.log 2>&1

    # CDS proteins are written as <predictor>.faa (possibly with prefix)
    faa=\$(ls ${genome_id}_pharokka/*${params.gene_predictor}.faa 2>/dev/null | head -n 1)
    if [ -z "\$faa" ]; then
        echo "No ${params.gene_predictor}.faa found in pharokka output" >&2
        exit 1
    fi
    cp "\$faa" ${genome_id}_proteins.faa
    """
}

process BAKTA {
    tag "${genome_id}"
    publishDir "${params.out_dir}/annotation/bakta", mode: 'copy'

    input:
    tuple val(genome_id), path(genome)

    output:
    tuple val(genome_id), path("${genome_id}_bakta"),         emit: annotation
    tuple val(genome_id), path("${genome_id}_proteins.faa"),  emit: proteins
    path("${genome_id}_bakta.log"),                           emit: log

    script:
    def db = params.bakta_db ? "--db ${params.bakta_db}" : ''
    """
    source ${params.conda_sh_ri}
    conda activate ${params.bakta_env}

    bakta ${db} --output ${genome_id}_bakta --prefix ${genome_id} \\
        --threads ${task.cpus} --force ${params.bakta_args} ${genome} \\
        > ${genome_id}_bakta.log 2>&1

    # locus tag -> contig from the .tsv, then rename proteins to <contig>_CDS_<n>
    awk -F'\\t' '
        NR == FNR { if (\$0 !~ /^#/) contig[\$6] = \$1; next }
        /^>/ {
            split(substr(\$0, 2), h, " "); id = h[1]
            c = (id in contig) ? contig[id] : "unknown_contig"
            n[c]++
            printf(">%s_CDS_%04d %s\\n", c, n[c], substr(\$0, 2)); next
        }
        { print }
    ' ${genome_id}_bakta/${genome_id}.tsv ${genome_id}_bakta/${genome_id}.faa > ${genome_id}_proteins.faa
    """
}