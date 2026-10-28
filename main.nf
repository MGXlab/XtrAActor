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
params.genome_type     = null          // "phage" (pharokka) or "non-phage" (bakta)
params.gene_predictor  = "prodigal-gv"   // phage only: "phanotate" or "prodigal"
params.pharokka_db     = "/vast/ri65fin/pharokka_db"          // optional if pharokka's DB is in its default location
params.bakta_db        = "/veodata/03/databases/bakta/v20240516/db"          // optional if BAKTA_DB is set in the environment
params.pharokka_args   = "--meta --skip_extra_annotations --skip_mash"
params.bakta_args      = "--skip-plot"
params.id              = 0.5           // MMseqs2 --min-seq-id (fraction)
params.cov             = 0.5           // MMseqs2 -c (fraction)

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

// --- tools / environments ---
params.mafft_bin       = "/home/groups/VEO/tools/mafft/v7.505/bin/mafft"
params.trimal_bin      = "/home/groups/VEO/tools/trimal/v1.5.0/trimal/source/trimal"
params.fasttree_bin    = "/home/groups/VEO/tools/fastTreeMP/v2.1.11/FastTreeMP"
params.conda_sh_veo    = "/vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh"
params.conda_sh_lydia  = "/vast/groups/VEO/shared_data/Lydia_Varada/miniconda3/etc/profile.d/conda.sh"
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
include { PLOTS }              from './modules/plots'

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
        --genome_type          'phage' (Pharokka) or 'non-phage' (Bakta) (required)
        --gene_predictor       Phage only: 'phanotate'  or 'prodigal-gv' (default = meta mode)
        --pharokka_db          Pharokka database dir (optional)
        --bakta_db             Bakta database dir (optional if BAKTA_DB is set)
        --pharokka_args        Extra Pharokka arguments (default: none)
        --bakta_args           Extra Bakta arguments (default: --skip-plot)
        --id                   MMseqs2 min. sequence identity, fraction (default 0.5)
        --cov                  MMseqs2 coverage, fraction (default 0.7)

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
        if (params.genome_type == 'phage' && !(params.gene_predictor in ['phanotate', 'prodigal','prodigal-gv']))
            errors << "--gene_predictor must be 'phanotate' or 'prodigal' for phage genomes"
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

    // global plots across all clusters
    rf_files_ch = CLUSTER_ANALYSIS.out.features
        .mix(CLUSTER_ANALYSIS.out.reports, CLUSTER_ANALYSIS.out.data_df)
        .map { name, f -> f }
        .flatten()
        .collect()

    plot_fastas_ch = CLUSTER_ANALYSIS.out.fastas
        .map { name, f -> f }
        .collect()

    PLOTS(
        rf_files_ch,
        plot_fastas_ch,
        params.biome1.tokenize(';').last(),
        params.biome2.tokenize(';').last(),
        params.colour1,
        params.colour2
    )
}