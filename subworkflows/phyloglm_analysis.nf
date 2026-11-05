include {
    CREATE_PHYLOGLM_DATASET
    SPLIT_GINI
    PHYLOGLM_FIT_CHUNK
    MERGE_PHYLOGLM
    PLOT_ZVALUE
} from '../modules/phyloglm'

workflow PHYLOGLM_ANALYSIS {
    take:
    features_ch     // tuple(cluster_name, feature_importance.tsv)
    data_df_ch      // tuple(cluster_name, data_df.tsv)
    trees_ch        // tuple(cluster_name, tree)

    main:
    CREATE_PHYLOGLM_DATASET(features_ch.join(data_df_ch))

    chunks_ch = SPLIT_GINI(CREATE_PHYLOGLM_DATASET.out.gini).chunks
        .flatMap { name, chunks ->
            (chunks instanceof List ? chunks : [chunks]).collect { c -> tuple(name, c) }
        }

    fit_input_ch = chunks_ch
        .combine(CREATE_PHYLOGLM_DATASET.out.dataset, by: 0)
        .combine(trees_ch, by: 0)
        .map { name, chunk, dataset, tree -> tuple(name, dataset, tree, chunk) }

    PHYLOGLM_FIT_CHUNK(fit_input_ch)
    MERGE_PHYLOGLM(PHYLOGLM_FIT_CHUNK.out.partial.groupTuple())

 // only plot clusters where at least one feature was fitted (more than the header line)
    csv_with_results = MERGE_PHYLOGLM.out.csv.filter { name, f ->
        def has_rows = f.readLines().size() > 1
        if (!has_rows) log.warn "No phyloglm results for cluster ${name}, skipping z-value plot"
        has_rows
    }
    PLOT_ZVALUE(csv_with_results.join(features_ch))   // -> (cluster, phyloglm.csv, Feature_importance.tsv)
    // PLOT_ZVALUE(MERGE_PHYLOGLM.out.csv)

    emit:
    csv = MERGE_PHYLOGLM.out.csv
}