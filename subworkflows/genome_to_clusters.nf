include { PHAROKKA; BAKTA }     from '../modules/annotation'
include { MMSEQ2; ORTHOGROUPS } from '../modules/clustering'

workflow GENOME_TO_CLUSTERS {
    take:
    genomes_ch      // path: one genome fasta per item
    metadata_ch     // value channel: metadata csv

    main:
    named_ch = genomes_ch.map { f -> tuple(f.baseName, f) }

    def annot
    if (params.genome_type == 'phage') {
        annot = PHAROKKA(named_ch)
    } else {
        annot = BAKTA(named_ch)
    }

    // pool all genomes so MMseqs2 clusters across genomes
    all_proteins = annot.proteins.collectFile(name: 'all_proteins.faa', sort: { f -> f.name })

    MMSEQ2(all_proteins, params.id, params.cov)
    ORTHOGROUPS(MMSEQ2.out.tsv, all_proteins, metadata_ch)

    emit:
    fastas      = ORTHOGROUPS.out.cluster_fastas
                    .flatten()
                    .map { f -> tuple(f.baseName, f) }
    biomes      = ORTHOGROUPS.out.cluster_biomes
                    .flatten()
                    .map { f -> tuple(f.baseName.replaceFirst(/_biome$/, ''), f) }
    annotations = annot.annotation
}