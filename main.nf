#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================
// Parameters
// =============================================================
params.out_dir         = "output_dir_1"
params.entry_point     = null          // "genome" or "clusters"

// --- entry point: genome ---
params.genome_dir      = null          // dir with one genome fasta per file (*.fna/*.fa/*.fasta)
params.metadata_file   = null
params.precomputed_proteins = null    // mmseq2_derep folder of an earlier run (already dereplicated proteins of other biomes)
params.genome_type     = null          // "phage" (pharokka) or "non-phage" (bakta)
params.gene_predictor  = "prodigal-gv" // phage only: "phanotate", "prodigal" or "prodigal-gv"
params.pharokka_db     = "/vast/ri65fin/pharokka_db"          // optional if pharokka's DB is in its default location
params.bakta_db        = "/veodata/03/databases/bakta/v20240516/db"          // optional if BAKTA_DB is set in the environment
params.pharokka_args   = "--meta --skip_extra_annotations --skip_mash"
params.bakta_args      = "--skip-plot --gram -"
params.id              = 0.5           // MMseqs2 --min-seq-id (fraction)
params.cov             = 0.5           // MMseqs2 -c (fraction)
params.dereplicate     = true          // remove near-identical proteins within each biome first
params.derep_id        = 0.99          // MMseqs2 --min-seq-id for dereplication (fraction)
params.derep_cov       = 0.9           // MMseqs2 -c for dereplication (fraction)

// --- entry point: clusters ---
params.fasta_dir       = null
params.metadata_dir    = "metadata"

// --- downstream ---
params.min_seqs        = 4
params.do_trimal       = false
params.n_chunks        = 15
params.biome1          = null
params.biome2          = null
params.colour1         = null
params.colour2         = null
params.help            = false

// --- protein enrichment analysis with pyseer (entry point: genome) ---
params.enrichment             = true      // run the protein enrichment analysis
params.enrich_min_proteins    = 10        // test clusters with at least this many proteins (members in mmseq2_cluster.tsv)
params.genome_tree            = null      // optional Newick tree of all genomes (tip labels = genome IDs) -> pyseer LMM
params.enrich_distances       = null      // without a tree: 'viptreegen', 'mash' or 'dice' (default: phage viptreegen, non-phage mash)
params.precomputed_genomes    = null      // genome fastas of the --precomputed_proteins biome(s), for viptreegen / mash; glob allowed
params.pyseer_env             = "pyseer"  // conda env with pyseer and mash, in /vast/ri65fin/miniconda3 (conda_sh_ri)
params.pyseer_min_af          = 0.01      // pyseer --min-af: min. fraction of genomes carrying a cluster
params.pyseer_max_af          = 0.99      // pyseer --max-af
params.pyseer_max_dimensions  = 10        // MDS dimensions of the distance matrix (fixed-effects model)
params.pyseer_args            = ""        // extra pyseer arguments
params.enrich_alpha           = 0.05      // significance level (Bonferroni over distinct presence/absence patterns)
params.viptreegen_bin         = "/home/groups/VEO/tools/ViPTreeGen/v1.1.3/ViPTreeGen/ViPTreeGen"
params.viptreegen_conda_sh    = "/vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh"
params.viptreegen_env         = "ViPTreeGen_v1.1.3"                       // conda env for ViPTreeGen
params.viptreegen_blast_dir   = "/home/groups/VEO/tools/blast/v2.14.0/bin"  // BLAST (tblastx) added to the PATH
params.viptreegen_args        = ""        // extra ViPTreeGen arguments, e.g. "--mode mmseqs-tblastx" (newer versions)
params.mash_sketch_size       = 10000     // mash sketch -s
params.r_env                  = "phyloglm"    // conda env with R + ape (phangorn for midpoint rooting), for --genome_tree
params.diamond_bin            = "/home/groups/VEO/tools/diamond/v2.1.8/diamond"
params.diamond_env            = null      // conda env with diamond (activated via conda_sh_veo); null = use diamond_bin as is
params.diamond_evalue         = 1e-5      // same as dice_evalue: weaker hits are not used for the distances
params.diamond_pident         = 40        // diamond --id; same as dice_pident
params.diamond_args           = ""        // extra diamond blastp arguments
params.diamond_chunk_size     = 100000    // query proteins per blastp task
params.dice_evalue            = 1e-5      // dice_tree_fast.py -e
params.dice_pident            = 40        // dice_tree_fast.py -p
params.dice_args              = ""        // extra filters for the Dice step, e.g. "-qc 70 -tc 70"

// --- tools / environments ---
params.mafft_bin       = "/home/groups/VEO/tools/mafft/v7.505/bin/mafft"
params.trimal_bin      = "/home/groups/VEO/tools/trimal/v1.5.0/trimal/source/trimal"
params.fasttree_bin    = "/home/groups/VEO/tools/fastTreeMP/v2.1.11/FastTreeMP"
params.conda_sh_veo    = "/vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh"
params.conda_sh_ri     = "/vast/ri65fin/miniconda3/etc/profile.d/conda.sh"
params.pharokka_env    = "pharokka"
params.bakta_env       = "bakta_v1.9.3"
params.shared_bin      = "/vast/ri65fin/nextflow_AAcidx_v2/bin"

// =============================================================
// Includes
// =============================================================
include { GENOME_TO_CLUSTERS } from './subworkflows/genome_to_clusters'
include { CLUSTER_ANALYSIS }   from './subworkflows/cluster_analysis'
include { PHYLOGLM_ANALYSIS }  from './subworkflows/phyloglm_analysis'
include { PROTEIN_ENRICHMENT } from './subworkflows/protein_enrichment'
include { PLOT_RF_REPORT; PLOT_CLUSTER; CLUSTER_PROPERTY_PLOTS } from './modules/plots'
include { RANK_CLUSTERS }      from './modules/ranking'

// =============================================================
// Help & validation
// =============================================================
if (params.help) {
    log.info """
    Usage:
        # phage genomes
        nextflow run main.nf --entry_point genome --genome_dir genomes \\
            --genome_type phage --gene_predictor phanotate \\
            --biome1 "..." --biome2 "..." --colour1 "#0bd600" --colour2 "#8e00fa"

        # non-phage genomes
        nextflow run main.nf --entry_point genome --genome_dir genomes \\
            --genome_type non-phage --bakta_db /path/to/bakta/db \\
            --biome1 "..." --biome2 "..." --colour1 "#0bd600" --colour2 "#8e00fa"

        # pre-built clusters
        nextflow run main.nf --entry_point clusters --fasta_dir cluster_fastas \\
            --metadata_dir metadata --biome1 "..." --biome2 "..." --colour1 "#0bd600" --colour2 "#8e00fa"

    Options:
        --entry_point          'genome' or 'clusters' (required)
        --out_dir              Output directory (default: output_dir_1)

      genome entry point:
        --genome_dir           Directory of genome fasta files (*.fna, *.fa, *.fasta)
        --metadata_file        Metadata CSV passed to orthogroups.py
        --precomputed_proteins mmseq2_derep folder of an earlier run: the already dereplicated proteins
                               of other biomes (<biome>_derep_rep_seq.fasta + <biome>_derep_cluster.tsv),
                               added directly to the clustering
        --genome_type          'phage' (Pharokka) or 'non-phage' (Bakta) (required)
        --gene_predictor       Phage only: 'phanotate'  or 'prodigal-gv' (default = meta mode)
        --pharokka_db          Pharokka database dir (optional)
        --bakta_db             Bakta database dir (optional if BAKTA_DB is set)
        --pharokka_args        Extra Pharokka arguments (default: none)
        --bakta_args           Extra Bakta arguments (default: --skip-plot)
        --id                   MMseqs2 min. sequence identity, fraction (default 0.5)
        --cov                  MMseqs2 coverage, fraction (default 0.7)
        --dereplicate          Dereplicate proteins within each biome first (default true)
        --derep_id             Dereplication min. sequence identity, fraction (default 0.99)
        --derep_cov            Dereplication coverage, fraction (default 0.9)

      protein enrichment analysis with pyseer (genome entry point only):
        --enrichment           Run the protein enrichment analysis (default true)
        --enrich_min_proteins  Test clusters with at least this many proteins (default 10)
        --genome_tree          Newick tree of the genomes (tips = genome IDs); used instead of
                               distances, with pyseer's linear mixed model
        --enrich_distances     Without a tree: 'viptreegen', 'mash' or 'dice'
                               (default: viptreegen for phages, mash for non-phage)
        --precomputed_genomes  Genome fastas of the --precomputed_proteins biome(s) (viptreegen / mash)
        --pyseer_min_af / --pyseer_max_af   pyseer frequency filter (default 0.01 / 0.99)
        --pyseer_max_dimensions  MDS dimensions for the distance matrix (default 10)
        --enrich_alpha         Significance level (default 0.05)
        --pyseer_env / --r_env                      conda environments (in --conda_sh_ri)
        --viptreegen_bin / --viptreegen_args        ViPTreeGen executable / extra arguments
        --viptreegen_conda_sh / --viptreegen_env    conda installation and env for ViPTreeGen
        --viptreegen_blast_dir                      BLAST bin folder added to the PATH for ViPTreeGen
        --diamond_bin / --diamond_env               diamond (dice only)
        --diamond_evalue / --diamond_pident         all-v-all blastp e-value / min. identity (default 1e-5 / 40)
        --diamond_chunk_size   Query proteins per blastp task (default 100000)
        --dice_evalue / --dice_pident / --dice_args  Dice thresholds (default 1e-5 / 40 / none)

      clusters entry point:
        --fasta_dir            Per-cluster fasta files (<cluster>.faa or .fasta)
        --metadata_dir         Per-cluster tables (<cluster>_biome.tsv or <cluster>_metadata.tsv)

      downstream:
        --min_seqs             Skip clusters with fewer sequences (default 4)
        --do_trimal            Run trimAl on the MAFFT MSAs (default false)
        --n_chunks             Max phyloglm chunks per cluster (default 15)
        --biome1 / --biome2    Biome labels (required)
        --colour1 / --colour2  Plot colours (required)
    """
    exit 0
}

def validateParams() {
    def errors = []
    if (!(params.entry_point in ['genome', 'clusters']))
        errors << "--entry_point must be 'genome' or 'clusters' (got: ${params.entry_point})"

    if (params.entry_point == 'genome') {
        if (!params.genome_dir)    errors << "--genome_dir is required for --entry_point genome"
        if (!params.metadata_file) errors << "--metadata_file is required for --entry_point genome"
        if (!(params.genome_type in ['phage', 'non-phage']))
            errors << "--genome_type must be 'phage' or 'non-phage' (got: ${params.genome_type})"
        if (params.genome_type == 'phage' && !(params.gene_predictor in ['phanotate', 'prodigal', 'prodigal-gv']))
            errors << "--gene_predictor must be 'phanotate', 'prodigal' or 'prodigal-gv' for phage genomes"
        if (params.genome_dir) {
            def n_files = files("${params.genome_dir}/*.{fna,fa,fasta}").size()
            def n_pre   = params.precomputed_proteins ? files("${params.precomputed_proteins}/*_rep_seq.fasta").size() : 0
            if (n_files < 1)
                errors << "--genome_dir contains no fasta files"
            if (params.precomputed_proteins && !file(params.precomputed_proteins).isDirectory())
                errors << "--precomputed_proteins must be a folder (mmseq2_derep of an earlier run): ${params.precomputed_proteins}"
            else if (params.precomputed_proteins && n_pre == 0)
                errors << "--precomputed_proteins contains no *_rep_seq.fasta: ${params.precomputed_proteins}"
            if (n_files + n_pre < 2)
                errors << "need at least two biome sets: fastas in --genome_dir plus biomes in --precomputed_proteins"
        }
        if (params.enrichment) {
            def method = params.enrich_distances ?: (params.genome_type == 'phage' ? 'viptreegen' : 'mash')
            if (!params.genome_tree && !(method in ['viptreegen', 'mash', 'dice']))
                errors << "--enrich_distances must be 'viptreegen', 'mash' or 'dice' (got: ${method})"
            if (params.genome_tree && !file(params.genome_tree).exists())
                errors << "--genome_tree not found: ${params.genome_tree}"
            if (!params.genome_tree && method in ['viptreegen', 'mash'] && params.precomputed_proteins && !params.precomputed_genomes)
                errors << "--enrich_distances ${method} needs the genomes of the precomputed biome(s): " +
                          "give --precomputed_genomes, a --genome_tree, or use --enrich_distances dice"
        }
    }
    if (params.entry_point == 'clusters') {
        if (!params.fasta_dir)    errors << "--fasta_dir is required for --entry_point clusters"
        if (!params.metadata_dir) errors << "--metadata_dir is required for --entry_point clusters"
    }
    ['biome1', 'biome2', 'colour1', 'colour2'].each { p ->
        if (!params[p]) errors << "--${p} is required"
    }
    if (errors) error("Parameter errors:\n  - " + errors.join("\n  - "))
}

// =============================================================
// Main workflow
// =============================================================
workflow {
    validateParams()

    def fastas_ch
    def biomes_ch

    if (params.entry_point == 'genome') {
        genomes_ch  = Channel.fromPath("${params.genome_dir}/*.{fna,fa,fasta}", checkIfExists: true)
        metadata_ch = Channel.value(file(params.metadata_file, checkIfExists: true))

        GENOME_TO_CLUSTERS(genomes_ch, metadata_ch)
        fastas_ch = GENOME_TO_CLUSTERS.out.fastas
        biomes_ch = GENOME_TO_CLUSTERS.out.biomes

        // protein enrichment (pyseer); runs in parallel with the multi-biome analysis below
        if (params.enrichment) {
            rf_clusters_ch = GENOME_TO_CLUSTERS.out.fastas
                .map { name, f -> name }
                .collect()
                .ifEmpty([])
            PROTEIN_ENRICHMENT(
                genomes_ch,
                GENOME_TO_CLUSTERS.out.proteins,
                GENOME_TO_CLUSTERS.out.clustering_fasta,
                GENOME_TO_CLUSTERS.out.cluster_tsv,
                GENOME_TO_CLUSTERS.out.derep_tsv,
                rf_clusters_ch,
                metadata_ch
            )
        }
    } else {
        fastas_ch = Channel.fromPath("${params.fasta_dir}/*.{faa,fasta}", checkIfExists: true)
            .map { f -> tuple(f.baseName, f) }
        biomes_ch = Channel.fromPath("${params.metadata_dir}/*_{biome,metadata}.tsv", checkIfExists: true)
            .map { f -> tuple(f.baseName.replaceFirst(/_(biome|metadata)$/, ''), f) }
    }

    CLUSTER_ANALYSIS(fastas_ch, biomes_ch)

    PHYLOGLM_ANALYSIS(
        CLUSTER_ANALYSIS.out.features,
        CLUSTER_ANALYSIS.out.data_df,
        CLUSTER_ANALYSIS.out.trees
    )

    // plots per cluster: each starts as soon as that cluster's RF model is done
    // (cluster, [Feature_importance, RandomForest_Report, data_df], fasta)
    plot_cluster_ch = CLUSTER_ANALYSIS.out.features
        .join(CLUSTER_ANALYSIS.out.reports)
        .join(CLUSTER_ANALYSIS.out.data_df)
        .join(CLUSTER_ANALYSIS.out.fastas)
        .map { name, feat, report, data, fasta -> tuple(name, [feat, report, data].flatten(), fasta) }

    PLOT_CLUSTER(
        plot_cluster_ch,
        params.biome1,
        params.biome2,
        params.colour1,
        params.colour2
    )

    // RF report plots across all clusters (after the last RF model)
    PLOT_RF_REPORT(
        CLUSTER_ANALYSIS.out.reports.map { name, f -> f }.flatten().collect(),
        params.biome1,
        params.biome2,
        params.colour1,
        params.colour2
    )

    plot_fastas_ch = CLUSTER_ANALYSIS.out.fastas
        .map { name, f -> f }
        .collect()

    // protein properties compared across all clusters: one volcano plot per property
    data_df_ch = CLUSTER_ANALYSIS.out.data_df
        .map { name, f -> f }
        .collect()

    CLUSTER_PROPERTY_PLOTS(
        data_df_ch,
        plot_fastas_ch,
        params.biome1,
        params.biome2,
        params.colour1,
        params.colour2
    )

    // rank clusters by RF accuracy, with annotations from the cluster fasta headers
    reports_ch = CLUSTER_ANALYSIS.out.reports
        .map { name, f -> f }
        .collect()

    RANK_CLUSTERS(reports_ch, plot_fastas_ch)
}