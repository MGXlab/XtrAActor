// Protein IDs are <genome>_CDS_<n> for both tools, so downstream steps can take the genome
// from the part before "_CDS_":
//   PHAROKKA: one phage contig = one genome, the genome is the contig name
//   BAKTA:    one input file  = one genome (any number of contigs), the genome is the file name
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
    path("${genome_id}_protein2contig.tsv"),                  emit: protein2contig
    path("${genome_id}_bakta.log"),                           emit: log

    script:
    def db = params.bakta_db ? "--db ${params.bakta_db}" : ''
    """
    source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh
    conda activate ${params.bakta_env}

    bakta ${db} --output ${genome_id}_bakta --prefix ${genome_id} \\
        --threads ${task.cpus} --force ${params.bakta_args} ${genome} \\
        > ${genome_id}_bakta.log 2>&1

    # rename proteins to <genome_id>_CDS_<n> (one input file = one genome, numbered in file order)
    # header text = product only; protein -> contig + locus tag goes to <genome_id>_protein2contig.tsv
    awk -F'\\t' -v g="${genome_id}" -v map="${genome_id}_protein2contig.tsv" '
        BEGIN { print "protein\\tcontig\\tlocus_tag" > map }
        NR == FNR { if (\$0 !~ /^#/) contig[\$6] = \$1; next }
        /^>/ {
            split(substr(\$0, 2), h, " "); id = h[1]
            product = substr(\$0, length(id) + 3)
            n++
            pid = sprintf("%s_CDS_%05d", g, n)
            printf(">%s %s\\n", pid, product)
            print pid "\\t" ((id in contig) ? contig[id] : "unknown_contig") "\\t" id > map
            next
        }
        { print }
    ' ${genome_id}_bakta/${genome_id}.tsv ${genome_id}_bakta/${genome_id}.faa > ${genome_id}_proteins.faa
    """
}