include { MAFFT; TRIMAL; FASTTREE } from '../modules/alignment'
include { RF_MODEL }                from '../modules/rf_model'

workflow CLUSTER_ANALYSIS {
    take:
    fastas_ch       // tuple(cluster_name, fasta)
    biomes_ch       // tuple(cluster_name, biome_tsv)

    main:
    // skip clusters too small for MSA / tree / RF
    filtered_ch = fastas_ch.filter { name, f ->
        def n = f.countFasta()
        if (n < params.min_seqs) {
            log.warn "Skipping cluster ${name}: ${n} sequence(s) < --min_seqs ${params.min_seqs}"
            return false
        }
        return true
    }

    MAFFT(filtered_ch)
    msa_ch = params.do_trimal ? TRIMAL(MAFFT.out.aligned).faa : MAFFT.out.aligned

    FASTTREE(msa_ch)

    rf_in_ch = msa_ch
        .join(biomes_ch, remainder: true)
        .filter { name, msa, biome ->
            if (msa == null) return false
            if (biome == null) {
                log.warn "No biome table found for cluster ${name}, skipping RF_MODEL"
                return false
            }
            return true
        }
    RF_MODEL(rf_in_ch)

    emit:
    fastas   = filtered_ch
    trees    = FASTTREE.out.tree
    features = RF_MODEL.out.features
    reports  = RF_MODEL.out.reports
    data_df  = RF_MODEL.out.data_df
}