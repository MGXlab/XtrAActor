include {
    DIAMOND_MAKEDB
    DIAMOND_BLASTP
    DICE_TREE
    SINGLE_BIOME_CLUSTERS
    PHYLOGLM_SINGLE_BIOME
    SINGLE_BIOME_REPORT
} from '../modules/single_biome'

// Runs next to CLUSTER_ANALYSIS / PHYLOGLM_ANALYSIS (the multi-biome analysis); nothing in
// here waits for MAFFT / FastTree / RF, and the tree branch does not wait for the clustering.
//
//   predicted proteins ─ DIAMOND_MAKEDB ─ DIAMOND_BLASTP (chunks) ─ DICE_TREE ───────────┐
//   mmseq2_cluster.tsv ─ SINGLE_BIOME_CLUSTERS ──────────────────────────────────────────┴ PHYLOGLM_SINGLE_BIOME (chunks) ─ SINGLE_BIOME_REPORT

workflow SINGLE_BIOME_ANALYSIS {
    take:
    proteins_ch       // tuple(biome_set, proteins.faa): predicted proteins BEFORE dereplication
    cluster_tsv_ch    // mmseq2_cluster.tsv of the 50/50 clustering
    derep_tsv_ch      // value channel: list of <biome_set>_derep_cluster.tsv ([] without dereplication)
    metadata_ch       // value channel: combined metadata file (all biomes)

    main:
    // proteins of biomes that were not annotated in this run; without them the genomes
    // of that biome would be missing from the tree
    precomputed_ch = (params.precomputed_proteins && params.sb_include_precomputed)
        ? Channel.fromPath(params.precomputed_proteins, checkIfExists: true)
        : Channel.empty()

    all_genes_ch = proteins_ch
        .map { name, f -> f }
        .mix(precomputed_ch)
        .collectFile(name: 'all_genes.faa', sort: { f -> f.name })
        .first()

    // 1. all-v-all blastp
    DIAMOND_MAKEDB(all_genes_ch)
    query_chunks_ch = all_genes_ch.splitFasta(by: params.diamond_chunk_size, file: true)
    DIAMOND_BLASTP(query_chunks_ch, DIAMOND_MAKEDB.out.db)

    // 2. Dice distances + genome tree
    DICE_TREE(DIAMOND_BLASTP.out.hits.collect())

    // 3. single-biome clusters (independent of the blastp / tree)
    SINGLE_BIOME_CLUSTERS(cluster_tsv_ch, derep_tsv_ch, all_genes_ch, metadata_ch)

    // 4. phyloglm of each single-biome cluster against the genome tree
    cluster_chunks_ch = SINGLE_BIOME_CLUSTERS.out.clusters
        .splitText(by: params.sb_chunk_size, keepHeader: true, file: true)
    PHYLOGLM_SINGLE_BIOME(cluster_chunks_ch, DICE_TREE.out.tree, metadata_ch)

    // 5. FDR over all clusters, keep padj < params.sb_padj
    SINGLE_BIOME_REPORT(PHYLOGLM_SINGLE_BIOME.out.partial.collect())

    emit:
    tree      = DICE_TREE.out.tree
    distances = DICE_TREE.out.distances
    clusters  = SINGLE_BIOME_CLUSTERS.out.clusters
    all       = SINGLE_BIOME_REPORT.out.all
    report    = SINGLE_BIOME_REPORT.out.report
}
