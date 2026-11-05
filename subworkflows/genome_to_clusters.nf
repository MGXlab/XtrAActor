include { PHAROKKA; BAKTA }                     from '../modules/annotation'
include { MMSEQ2_DEREP; MMSEQ2; ORTHOGROUPS }   from '../modules/clustering'

workflow GENOME_TO_CLUSTERS {
    take:
    genomes_ch      // path: one multi-fasta per NEW biome (annotated + dereplicated in this run)
    metadata_ch     // value channel: combined metadata file (all biomes)

    main:
    biome_sets_ch = genomes_ch.map { f -> tuple(f.baseName, f) }   // e.g. ('alkaline', alkaline.fna)

    // 1. gene calling / annotation for the new biome set(s)
    def annot
    if (params.genome_type == 'phage') {
        annot = PHAROKKA(biome_sets_ch)
    } else {
        annot = BAKTA(biome_sets_ch)
    }

    // 2. dereplicate within each new biome set
    def new_proteins_ch
    if (params.dereplicate) {
        MMSEQ2_DEREP(annot.proteins, params.derep_id, params.derep_cov)
        new_proteins_ch = MMSEQ2_DEREP.out.rep_fasta.map { name, f -> f }
    } else {
        new_proteins_ch = annot.proteins.map { name, f -> f }
    }

    // 3. already dereplicated proteins of other biomes (e.g. freshwater), used as they are
    precomputed_ch = params.precomputed_proteins
        ? Channel.fromPath(params.precomputed_proteins, checkIfExists: true)
        : Channel.empty()

    // 4. pool everything for the cross-biome clustering
    proteins_for_clustering = new_proteins_ch
        .mix(precomputed_ch)
        .collectFile(name: 'proteins_for_clustering.faa', sort: { f -> f.name })

    MMSEQ2(proteins_for_clustering, params.id, params.cov)
    ORTHOGROUPS(MMSEQ2.out.tsv, proteins_for_clustering, metadata_ch)

    emit:
    fastas      = ORTHOGROUPS.out.cluster_fastas
                    .flatten()
                    .map { f -> tuple(f.baseName, f) }
    biomes      = ORTHOGROUPS.out.cluster_biomes
                    .flatten()
                    .map { f -> tuple(f.baseName.replaceFirst(/_biome$/, ''), f) }
    annotations = annot.annotation
}