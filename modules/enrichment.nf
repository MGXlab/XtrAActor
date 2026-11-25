// =============================================================
// Protein enrichment analysis (pyseer)
//   protein clusters (presence/absence per genome) ~ biome,
//   corrected for population structure from one of:
//     --genome_tree           user tree -> similarity matrix -> pyseer linear mixed model
//     ViPTreeGen  (phage)     1 - SG score            -> pyseer fixed effects (MDS)
//     mash        (non-phage) mash distances          -> pyseer fixed effects (MDS)
//     dice                    DIAMOND + Dice distance -> pyseer fixed effects (MDS)
// =============================================================

def activateEnv(env) {
    return "source ${params.conda_sh_ri}\n    conda activate ${env}"
}

// -------------------------------------------------------------
// E1. presence/absence matrix, phenotypes, cluster info
// -------------------------------------------------------------
process ENRICH_INPUTS {
    cpus 1
    memory { 30.GB * task.attempt }
    time   { 2.h * task.attempt }

    publishDir "${params.out_dir}/protein_enrichment/inputs", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.Rtab') ? null : fn }       // the matrix can be very large
    queue 'standard,short,long,fat'

    input:
    path(cluster_tsv)
    path(derep_tsvs, stageAs: 'derep/*')    // empty list without dereplication
    path(clustered_proteins)
    path(metadata_file)
    val(rf_clusters)                         // cluster names of the multi-biome (RF) analysis

    output:
    path "clusters.Rtab",            emit: pres
    path "phenotypes.tsv",           emit: pheno
    path "cluster_info.tsv",         emit: info
    path "enrich_inputs.log",        emit: log

    script:
    def rf_list = rf_clusters ? rf_clusters.join('\n') : ''
    """
    source ${params.conda_sh_ri}
    conda activate working-env

cat <<'RF_CLUSTERS' > rf_clusters.txt
${rf_list}
RF_CLUSTERS

    python3 ${projectDir}/bin/enrichment_inputs.py presence \\
        --clusters ${cluster_tsv} --metadata ${metadata_file} \\
        --biome1 "${params.biome1}" --biome2 "${params.biome2}" \\
        --derep_dir derep --proteins_fasta ${clustered_proteins} \\
        --rf_clusters rf_clusters.txt --min_proteins ${params.enrich_min_proteins} \\
        2> enrich_inputs.log
    """
}

// -------------------------------------------------------------
// E2a. user tree -> similarity matrix (for --lmm)
// -------------------------------------------------------------
process TREE_SIMILARITY {
    cpus 1
    memory { 30.GB * task.attempt }
    time   { 2.h * task.attempt }

    publishDir "${params.out_dir}/protein_enrichment/structure", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? fn : null }
    queue 'standard,short,long,fat'

    input:
    path(tree)
    path(phenotypes)

    output:
    path "similarity.tsv",       emit: similarity
    path "tree_similarity.log",  emit: log

    script:
    """
    ${activateEnv(params.r_env)}

    tail -n +2 ${phenotypes} | cut -f1 > genome_ids.txt
    Rscript ${projectDir}/bin/tree_to_similarity.R ${tree} genome_ids.txt similarity.tsv \\
        > tree_similarity.log 2>&1
    """
}

// -------------------------------------------------------------
// E2b. genomes -> one record per genome (ViPTreeGen, mash for phages)
// -------------------------------------------------------------
process PREP_GENOMES {
    cpus 1
    memory { 8.GB * task.attempt }
    time   { 2.h * task.attempt }
    queue 'standard,short,long,fat'

    input:
    path(genomes, stageAs: 'genomes/*')

    output:
    path "genomes.fna",       emit: fasta
    path "genome_names.tsv",  emit: names

    script:
    // phage: one contig = one genome; non-phage: one file = one genome
    def per_file = params.genome_type == 'phage' ? '' : '--one_per_file'
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/enrichment_inputs.py genomes -i genomes/* ${per_file} \\
        -o genomes.fna --names genome_names.tsv
    """
}

// -------------------------------------------------------------
// E2c. ViPTreeGen distance matrix (1 - SG score)
// -------------------------------------------------------------
process VIPTREEGEN {
    cpus 30
    memory { 60.GB * task.attempt }
    time   { 24.h * task.attempt }
    queue 'standard,long,fat'

    publishDir "${params.out_dir}/protein_enrichment/structure", mode: 'copy',
        saveAs: { fn -> fn == 'viptreegen.log' ? fn : null }

    input:
    path(genomes_fasta)
    path(names)

    output:
    path "distances.tsv",   emit: distances
    path "viptreegen.log",  emit: log

    script:
    """
    # ViPTreeGen: its own conda installation + BLAST (tblastx) on the PATH
    source ${params.viptreegen_conda_sh}
    conda activate ${params.viptreegen_env}
    export PATH=\$PATH:${params.viptreegen_blast_dir}

    ${params.viptreegen_bin} --notree --ncpus ${task.cpus} ${params.viptreegen_args} \\
        ${genomes_fasta} viptree_out > viptreegen.log 2>&1

    # back to the pyseer env for renaming the matrix to genome IDs
    conda deactivate

    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/enrichment_inputs.py viptree \\
        -i viptree_out/result/all.dist.matrix --names ${names} -o distances.tsv
    """
}

// -------------------------------------------------------------
// E2d. mash distances
// -------------------------------------------------------------
process MASH_DISTANCES {
    cpus 8
    memory { 30.GB * task.attempt }
    time   { 6.h * task.attempt }
    queue 'standard,short,long,fat'

    publishDir "${params.out_dir}/protein_enrichment/structure", mode: 'copy',
        saveAs: { fn -> fn == 'mash.log' ? fn : null }

    input:
    path(genomes, stageAs: 'genomes/*')   // non-phage: one file per genome; phage: genomes.fna
    path(names)                           // phage: genome_names.tsv; non-phage: NO_NAMES placeholder

    output:
    path "distances.tsv", emit: distances
    path "mash.log",      emit: log

    script:
    def phage = params.genome_type == 'phage'
    def sketch_in = phage ? "-i genomes/*" : "genomes/*"
    def name_arg  = phage ? "--names ${names}" : ''
    """
    ${activateEnv(params.pyseer_env)}

    mash sketch -p ${task.cpus} -s ${params.mash_sketch_size} -o sketch ${sketch_in} > mash.log 2>&1
    mash dist -p ${task.cpus} sketch.msh sketch.msh \\
        | python3 ${projectDir}/bin/enrichment_inputs.py mash ${name_arg} -o distances.tsv 2>> mash.log
    """
}

// -------------------------------------------------------------
// E2e. DIAMOND all-v-all + Dice distances
// -------------------------------------------------------------
process DIAMOND_MAKEDB {
    cpus 8
    memory { 20.GB * task.attempt }
    time   { 2.h * task.attempt }
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

process DIAMOND_BLASTP {
    tag "${query_chunk.baseName}"
    cpus 20
    memory { 60.GB * task.attempt }
    time   { 12.h * task.attempt }
    maxForks 20
    queue 'standard,short,long,fat'

    input:
    path(query_chunk)
    path(db)

    output:
    path "${query_chunk.baseName}_diamond.out", emit: hits
    path "${query_chunk.baseName}_diamond.log", emit: log

    script:
    def activate = params.diamond_env ? "source ${params.conda_sh_veo}\n    conda activate ${params.diamond_env}" : ''
    // hits below this identity are thrown away by the Dice step anyway (dice_pident)
    def min_id   = params.diamond_pident ? "--id ${params.diamond_pident}" : ''
    """
    ${activate}

    ${params.diamond_bin} blastp -d ${db} -q ${query_chunk} -o ${query_chunk.baseName}_diamond.out \\
        -k 0 --evalue ${params.diamond_evalue} ${min_id} --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \\
        ${params.diamond_args} > ${query_chunk.baseName}_diamond.log 2>&1
    """
}

process SUM_BITSCORES {
    tag "${hits.baseName}"
    cpus 1
    memory { 16.GB * task.attempt }
    time   { 2.h * task.attempt }
    maxForks 20
    queue 'standard,short,long,fat'

    input:
    path(hits)

    output:
    path "${hits.baseName}_pairsums.tsv", emit: sums
    path "${hits.baseName}_pairsums.log", emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/dice_tree_fast.py sum -i ${hits} \\
        -o ${hits.baseName}_pairsums.tsv \\
        -e ${params.dice_evalue} -p ${params.dice_pident} ${params.dice_args} \\
        > ${hits.baseName}_pairsums.log 2>&1
    """
}

process DICE_DISTANCES {
    cpus 2
    memory { 60.GB * task.attempt }
    time   { 6.h * task.attempt }
    queue 'standard,short,long,fat'

    publishDir "${params.out_dir}/protein_enrichment/structure", mode: 'copy',
        saveAs: { fn -> fn == 'distances.tsv' ? null : fn }

    input:
    path(sums, stageAs: 'sums/*')

    output:
    path "distances.tsv",       emit: distances
    path "dice_distances.tsv",  emit: pairs        // long table, pairs with dice < 1 only
    path "dice_distances.log",  emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/dice_tree_fast.py tree -i sums/*_pairsums.tsv \\
        -out_csv dice_distances.tsv --out_matrix distances.tsv --no_tree \\
        > dice_distances.log 2>&1
    """
}

// -------------------------------------------------------------
// E3. pyseer
// -------------------------------------------------------------
process PYSEER {
    cpus 8
    memory { 60.GB * task.attempt }
    time   { 12.h * task.attempt }
    queue 'standard,short,long,fat'

    publishDir "${params.out_dir}/protein_enrichment", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? "logs/${fn}" : fn }

    input:
    path(pres)
    path(phenotypes)
    path(structure)       // distances.tsv (fixed effects) or similarity.tsv (--lmm)
    val(model)            // 'fixed' or 'lmm'

    output:
    path "pyseer.assoc",   emit: assoc
    path "patterns.txt",   emit: patterns
    path "pyseer.log",     emit: log

    script:
    def structure_args = model == 'lmm'
        ? "--lmm --similarity ${structure}"
        : "--distances ${structure} --max-dimensions ${params.pyseer_max_dimensions}"
    """
    ${activateEnv(params.pyseer_env)}

    pyseer --phenotypes ${phenotypes} --pres ${pres} ${structure_args} \\
        --min-af ${params.pyseer_min_af} --max-af ${params.pyseer_max_af} \\
        --cpu ${task.cpus} --output-patterns patterns.txt ${params.pyseer_args} \\
        > pyseer.assoc 2> pyseer.log
    """
}

// -------------------------------------------------------------
// E4. report: pattern-based Bonferroni threshold, BH for reference
// -------------------------------------------------------------
process ENRICHMENT_REPORT {
    cpus 1
    memory { 8.GB * task.attempt }
    time   { 1.h * task.attempt }
    queue 'standard,short,long,fat'

    publishDir "${params.out_dir}/protein_enrichment", mode: 'copy'

    input:
    path(assoc)
    path(patterns)
    path(info)

    output:
    path "enrichment_all.tsv",     emit: all
    path "enrichment_report.tsv",  emit: report
    path "enrichment_summary.txt", emit: summary

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/enrichment_report.py --pyseer ${assoc} --patterns ${patterns} \\
        --info ${info} --biome1 "${params.biome1}" --biome2 "${params.biome2}" \\
        --alpha ${params.enrich_alpha}
    """
}
