include {
    ENRICH_INPUTS
    TREE_SIMILARITY
    PREP_GENOMES
    VIPTREEGEN
    MASH_DISTANCES
    DIAMOND_MAKEDB
    DIAMOND_BLASTP
    SUM_BITSCORES
    DICE_DISTANCES
    PYSEER
    ENRICHMENT_REPORT
} from '../modules/enrichment'

// Protein enrichment analysis: which protein clusters are enriched in one biome, after correcting
// for population structure. Runs next to the multi-biome analysis (MAFFT / RF / phyloglm); only the
// final presence/absence step waits for the list of clusters that went into the RF analysis.
//
//   clusters ─ ENRICH_INPUTS (presence/absence, phenotypes) ─────────────────────────┐
//   population structure, one of:                                                    │
//     --genome_tree              ─ TREE_SIMILARITY ─────────────── (similarity, LMM) ┤
//     --enrich_distances viptreegen ─ PREP_GENOMES ─ VIPTREEGEN ── (distances)       ├─ PYSEER ─ ENRICHMENT_REPORT
//     --enrich_distances mash       ─ (PREP_GENOMES) ─ MASH_DISTANCES                 │
//     --enrich_distances dice       ─ DIAMOND ─ SUM_BITSCORES ─ DICE_DISTANCES ──────┘

workflow PROTEIN_ENRICHMENT {
    take:
    genomes_ch            // path: genome fastas of this run (phage: per biome; non-phage: per genome)
    proteins_ch           // tuple(biome_set, proteins.faa) of ALL biomes, before dereplication (Dice only)
    clustered_proteins    // proteins_for_clustering.faa (annotation column)
    cluster_tsv_ch        // mmseq2_cluster.tsv
    derep_tsv_ch          // value channel: list of *_derep_cluster.tsv ([] without dereplication)
    rf_clusters_ch        // value channel: list of cluster names in the multi-biome (RF) analysis
    metadata_ch           // value channel: metadata file

    main:
    ENRICH_INPUTS(cluster_tsv_ch, derep_tsv_ch, clustered_proteins, metadata_ch, rf_clusters_ch)

    // genomes for ViPTreeGen / mash: this run's genomes + --precomputed_genomes
    all_genomes_ch = genomes_ch.mix(
        params.precomputed_genomes ? Channel.fromPath(params.precomputed_genomes, checkIfExists: true)
                                   : Channel.empty())

    def structure_ch
    def model
    if (params.genome_tree) {
        // user tree: shared branch length from the root -> linear mixed model
        TREE_SIMILARITY(file(params.genome_tree, checkIfExists: true), ENRICH_INPUTS.out.pheno)
        structure_ch = TREE_SIMILARITY.out.similarity
        model = 'lmm'
    } else {
        def method = params.enrich_distances ?: (params.genome_type == 'phage' ? 'viptreegen' : 'mash')
        model = 'fixed'
        if (method == 'viptreegen') {
            PREP_GENOMES(all_genomes_ch.collect())
            VIPTREEGEN(PREP_GENOMES.out.fasta, PREP_GENOMES.out.names)
            structure_ch = VIPTREEGEN.out.distances
        } else if (method == 'mash') {
            if (params.genome_type == 'phage') {          // one contig = one genome: sketch each sequence
                PREP_GENOMES(all_genomes_ch.collect())
                MASH_DISTANCES(PREP_GENOMES.out.fasta, PREP_GENOMES.out.names)
            } else {                                      // one file = one genome
                no_names = Channel.of('').collectFile(name: 'NO_NAMES')
                MASH_DISTANCES(all_genomes_ch.collect(), no_names)
            }
            structure_ch = MASH_DISTANCES.out.distances
        } else if (method == 'dice') {
            all_genes_ch = proteins_ch
                .map { name, f -> f }
                .collectFile(name: 'all_genes.faa', sort: { f -> f.name })
                .first()
            DIAMOND_MAKEDB(all_genes_ch)
            DIAMOND_BLASTP(all_genes_ch.splitFasta(by: params.diamond_chunk_size, file: true),
                           DIAMOND_MAKEDB.out.db)
            SUM_BITSCORES(DIAMOND_BLASTP.out.hits)
            DICE_DISTANCES(SUM_BITSCORES.out.sums.collect())
            structure_ch = DICE_DISTANCES.out.distances
        } else {
            error "--enrich_distances must be 'viptreegen', 'mash' or 'dice' (got: ${method})"
        }
    }

    PYSEER(ENRICH_INPUTS.out.pres, ENRICH_INPUTS.out.pheno, structure_ch, model)
    ENRICHMENT_REPORT(PYSEER.out.assoc, PYSEER.out.patterns, ENRICH_INPUTS.out.info)

    emit:
    all     = ENRICHMENT_REPORT.out.all
    report  = ENRICHMENT_REPORT.out.report
    summary = ENRICHMENT_REPORT.out.summary
}
