include { PHAROKKA; BAKTA }                     from '../modules/annotation'
include { MMSEQ2_DEREP; MMSEQ2; ORTHOGROUPS }   from '../modules/clustering'

// genome -> biome from the metadata file (header line; genome <whitespace> biome), as in
// orthogroups_sel_by_biome.py
def readBiomes(metadata_file) {
    def lookup = [:]
    file(metadata_file).readLines().drop(1).each { line ->
        def parts = line.trim().split(/\s+/, 2)
        if (parts.size() == 2 && !lookup.containsKey(parts[0])) lookup[parts[0]] = parts[1].trim()
    }
    return lookup
}

// file-name-safe biome label: the part after the last ';' ("Hot (42-90C)" -> "Hot_42_90C");
// the full label if two biomes share that part
def biomeNames(biomes) {
    def safe  = { String s -> s.replaceAll(/[^A-Za-z0-9]+/, '_').replaceAll(/^_+|_+$/, '') }
    def short_ = biomes.collectEntries { b -> [ b, safe(b.split(';')[-1]) ] }
    def clash = short_.values().countBy { it }.any { k, v -> v > 1 }
    return clash ? biomes.collectEntries { b -> [ b, safe(b) ] } : short_
}

workflow GENOME_TO_CLUSTERS {
    take:
    genomes_ch      // path: phage: one multi-fasta per NEW biome; non-phage: one fasta per genome
    metadata_ch     // value channel: combined metadata file (all biomes)

    main:
    biome_sets_ch = genomes_ch.map { f -> tuple(f.baseName, f) }   // e.g. ('alkaline', alkaline.fna)

    // 1. gene calling / annotation for the new biome set(s)
    //    phage:     one file per biome, PHAROKKA proteins are already one file per biome
    //    non-phage: one file per genome, BAKTA proteins are pooled per biome (metadata file)
    def annot
    def biome_proteins_ch       // tuple(biome_set, proteins.faa), one per biome
    if (params.genome_type == 'phage') {
        annot = PHAROKKA(biome_sets_ch)
        biome_proteins_ch = annot.proteins
    } else {
        annot = BAKTA(biome_sets_ch)

        def genome2biome = readBiomes(params.metadata_file)
        def biome_name   = biomeNames(genome2biome.values().toList().unique())
        biome_proteins_ch = annot.proteins
            .filter { genome_id, f ->
                if (genome2biome.containsKey(genome_id)) return true
                log.warn "Genome ${genome_id} is not in ${params.metadata_file}: left out of dereplication and clustering"
                return false
            }
            .map { genome_id, f -> tuple("${biome_name[genome2biome[genome_id]]}_proteins.faa", f) }
            .collectFile(storeDir: "${params.out_dir}/annotation/by_biome", sort: { f -> f.name })
            .map { f -> tuple(f.name.replaceFirst(/_proteins\.faa$/, ''), f) }
    }

    // 2. dereplicate within each new biome set
    def new_proteins_ch
    if (params.dereplicate) {
        MMSEQ2_DEREP(biome_proteins_ch, params.derep_id, params.derep_cov)
        new_proteins_ch = MMSEQ2_DEREP.out.rep_fasta.map { name, f -> f }
    } else {
        new_proteins_ch = biome_proteins_ch.map { name, f -> f }
    }

    // 3. already dereplicated proteins of other biomes (e.g. freshwater), used as they are
    //    --precomputed_proteins = mmseq2_derep folder of an earlier run:
    //    <biome>_derep_rep_seq.fasta + <biome>_derep_cluster.tsv
    def pre_dir = params.precomputed_proteins ? file(params.precomputed_proteins) : null
    precomputed_ch = pre_dir
        ? Channel.fromPath("${pre_dir}/*_rep_seq.fasta", checkIfExists: true)
        : Channel.empty()

    // 4. pool everything for the cross-biome clustering
    proteins_for_clustering = new_proteins_ch
        .mix(precomputed_ch)
        .collectFile(name: 'proteins_for_clustering.faa', sort: { f -> f.name })

    MMSEQ2(proteins_for_clustering, params.id, params.cov)
    ORTHOGROUPS(MMSEQ2.out.tsv, proteins_for_clustering, metadata_ch)

    // ---- for PROTEIN_ENRICHMENT ----
    // dereplication tables of this run and of the precomputed biome(s); empty list if there are none
    new_derep_tsv_ch = params.dereplicate ? MMSEQ2_DEREP.out.tsv : Channel.empty()
    pre_derep_tsv_ch = pre_dir
        ? Channel.fromPath("${pre_dir}/*_cluster.tsv", checkIfExists: true)
        : Channel.empty()
    derep_tsv_ch = new_derep_tsv_ch.mix(pre_derep_tsv_ch).toList()

    // proteins BEFORE dereplication of the precomputed biome(s): <biome>_proteins.faa as published
    // by PHAROKKA / BAKTA. Looked for in annotation/<pharokka|bakta>/ of
    //   1. this run's --out_dir
    //   2. the folder that contains the mmseq2_derep folder
    // If it is in neither, the dereplicated proteins are used for the genome tree instead.
    def pre_full = []
    if (pre_dir) {
        files("${pre_dir}/*_rep_seq.fasta").each { rep ->
            def biome = rep.name.replaceFirst(/(_derep)?_rep_seq\.fasta$/, '')
            def annotation_dirs = [ file(params.out_dir), pre_dir.parent ]
            def full = annotation_dirs.collectMany { d -> files("${d}/annotation/*/${biome}_proteins.faa") }
            if (full) {
                log.info "Precomputed biome '${biome}': genome tree uses ${full[0]}"
                pre_full << tuple(biome, full[0])
            } else {
                log.warn "Precomputed biome '${biome}': no annotation/*/${biome}_proteins.faa in " +
                         "${annotation_dirs.join(' or ')}, the genome tree uses the dereplicated proteins (${rep.name})"
                pre_full << tuple(biome, rep)
            }
        }
    }
    all_proteins_ch = biome_proteins_ch.mix(Channel.fromList(pre_full))

    emit:
    proteins    = all_proteins_ch       // tuple(biome_set, proteins.faa) of ALL biomes, before dereplication
    cluster_tsv = MMSEQ2.out.tsv        // mmseq2_cluster.tsv of the 50/50 clustering
    clustering_fasta = proteins_for_clustering   // all (dereplicated) proteins that were clustered
    derep_tsv   = derep_tsv_ch
    fastas      = ORTHOGROUPS.out.cluster_fastas
                    .flatten()
                    .map { f -> tuple(f.baseName, f) }
    biomes      = ORTHOGROUPS.out.cluster_biomes
                    .flatten()
                    .map { f -> tuple(f.baseName.replaceFirst(/_biome$/, ''), f) }
    annotations = annot.annotation
}