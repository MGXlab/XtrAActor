#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================
// Parameters
// =============================================================
params.out_dir       = "output_dir_1"
params.fasta_dir     = null           // dir containing one *.faa file per cluster/group
params.metadata_dir  = "metadata"     // dir containing metadata.csv
params.biome1        = null           // e.g. "some;path;style;label1"
params.biome2        = null
params.colour1       = null
params.colour2       = null
params.do_trimal     = false          // whether to run trimAl on the mafft MSAs before FASTTREE/RF_MODEL
params.help          = false
params.gene_prediction = "prodigal"
params.entry_point     = null   // z.B. "genome" oder "clusters"
params.id              = 0.5     //identity in percent
params.cov             = 0.7     //coverage in percent

if (params.help) {
    log.info """
    Usage:
        nextflow run pipeline_merged.nf --fasta_dir cluster_fastas --metadata_dir metadata \\
            --biome1 "..." --biome2 "..." --colour1 "#1f77b4" --colour2 "#ff7f0e"

    Options:
        --out_dir              Output directory (default: output_dir_1)
        --fasta_dir            Directory of per-cluster/group fasta files (*.faa)
        --metadata_dir         Directory containing metadata.csv
        --biome1 / --biome2    Biome labels used for plotting/enrichment
        --colour1 / --colour2  Colours used for plotting
        --do_trimal            Whether to run trimAl on the mafft MSAs (default: false)
        --gene_prediction      Choose the tool, prodigal or phanotate
        --entry_point          Choose the entry point of the pipeline, genome or cluster
        --id                   Specify the ID for MMSeq2, default 0.5
        --cov                  Specifiy the coverage for MMSeq2, default 0.7
    """
    exit 0
}

// Prodigal or Phanotate, user decides
process PRODIGAL{
    cpus 10
    memory '50 GB'
    time '6h'
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/prodigal", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path(fasta_file)

    output:
    path("genbank_prodigal.gbk"), emit: genbank
    path("proteins_prodigal.faa"), emit: fasta
    path "prodigal.log", emit: log

    script:
    """
    prodigal="/home/groups/VEO/tools/prodigal/v2.6.3/prodigal"
    \$prodigal -i ${fasta_file} -o genbank_prodigal.gbk -a proteins_prodigal.faa 2> prodigal.log
    """
}

process PHANOTATE{
    cpus 10
    memory '50 GB'
    time '6h'
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/phanotate", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path(fasta_file)

    output:
    path("phanotate.faa"), emit: fasta
    path "phanotate.log", emit: log

    script:
    """
    source /home/groups/VEO/tools/phanotate/v1.5.0/myenv/bin/activate
    phanotate.py ${fasta_file} -o phanotate.faa -f fasta 2> phanotate.log
    """
}

process MMSEQ2{
    cpus 10
    memory '50 GB'
    time '6h'
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/mmseq2", mode: 'copy'

    queue 'standard,short,long,fat'

    input:
    path(prediction_file)
    val id
    val cov

    output:
    path "mmseq2_all_seqs.fasta", emit: all_fasta
    path "mmseq2_cluster.tsv", emit: tsv
    path "mmseq2_rep_seq.fasta", emit: rep_fasta
    path "mmseq2.log", emit: log

    script:
    """
    source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
    conda activate mmseq2_v14.7e284
    mmseqs easy-cluster ${prediction_file} mmseq2 tmp --min-seq-id ${id} -c ${cov} --threads ${task.cpus} 2> mmseq2.log
    """
}

process ORTHOGROUPS{
    cpus 10
    memory '50 GB'
    time '6h'
    maxForks 4

    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.out_dir}/cluster_fastas", mode: 'copy', pattern: "*.fasta"
    publishDir "${params.out_dir}/cluster_biomes", mode: 'copy', pattern: "*.tsv"

    queue 'standard,short,long,fat'

    input:
    path cluster_tsv
    path proteins_fasta
    path metadata_file

    output:
    path "*.fasta", emit: cluster_fastas
    path "*.tsv", emit: cluster_biomes
    path "orthogroups.log", emit: log

    script:
    """
    source /vast/groups/VEO/shared_data/Lydia_Varada/miniconda3/etc/profile.d/conda.sh
    conda activate working-env
    
    python3 ${projectDir}/bin/orthogroups.py \\
        ${cluster_tsv} ${proteins_fasta} ${metadata_file} \\
        --fasta_out_dir cluster_fastas --biome_out_dir cluster_biomes 2> orthogroups.log
    """
}

// // =============================================================
// // 1. MSA (mafft)
// // =============================================================
// process MAFFT {
//     cpus 10
//     memory '50 GB'
//     time '6h'
//     maxForks 4

//     errorStrategy 'retry'
//     maxRetries 2

//     publishDir "${params.out_dir}/mafft", mode: 'copy'

//     queue 'standard,short,long,fat'

//     input:
//     tuple val(cluster_name), path(fasta_file)

//     output:
//     tuple val(cluster_name), path("${cluster_name}_mafft_aa.fa"), emit: aligned
//     path "${cluster_name}_mafft.log", emit: log

//     script:
//     """
//     mafft="/home/groups/VEO/tools/mafft/v7.505/bin/mafft"

//     # Fasta split in ref.fa (longest sequence) and remain.fa (what remains)
//     awk -v reff="ref.fa" -v restf="remain.fa" '
//     /^>/ {
//         if (seq != "") {
//             n++; headers[n] = header; seqs[n] = seq
//         }
//         header = \$0
//         seq = ""
//         next
//     }
//     {
//         line = \$0
//         gsub(/\\r/, "", line)
//         gsub(/[ \\t]+\$/, "", line)
//         seq = seq line
//     }
//     END {
//         if (seq != "") { n++; headers[n] = header; seqs[n] = seq }
//         maxlen = 0
//         for (i = 1; i <= n; i++) {
//             if (length(seqs[i]) > maxlen) { maxlen = length(seqs[i]); best = i }
//         }
//         print headers[best]"\\n"seqs[best] > reff
//         for (i = 1; i <= n; i++) {
//             if (i != best) {
//                 print headers[i]"\\n"seqs[i] >> restf
//             }
//         }
//     }
//     ' ${fasta_file}

//     \$mafft --keeplength --add remain.fa ref.fa > ${cluster_name}_mafft_aa.fa 2> ${cluster_name}_mafft.log
//     """
// }

// // =============================================================
// // 2. Optional trimAl
// // =============================================================
// process TRIMAL {
//     cpus 4
//     memory '10 GB'
//     time '1h'

//     publishDir "${params.out_dir}/mafft", mode: 'copy'

//     queue 'standard,short,long,fat,gpu-veo'

//     input:
//     tuple val(cluster_name), path(msa_file)

//     output:
//     tuple val(cluster_name), path("${cluster_name}_mafft_aa_trimmed.fa"), emit: faa
//     path "${cluster_name}_trimal.log", emit: log

//     script:
//     """
//     /home/groups/VEO/tools/trimal/v1.5.0/trimal/source/trimal -in ${msa_file} \\
//            -out ${cluster_name}_mafft_aa_trimmed.fa \\
//            -gappyout \\
//            -colnumbering > ${cluster_name}_trimal.log 2>&1
//     """
// }

// // =============================================================
// // 3. FastTree
// // =============================================================
// process FASTTREE {
//     cpus 10
//     memory '50 GB'
//     time '6h'
//     maxForks 4

//     errorStrategy 'retry'
//     maxRetries 2

//     publishDir "${params.out_dir}/fasttree", mode: 'copy'

//     queue 'standard,short,long,fat'

//     input:
//     tuple val(cluster_name), path(msa)

//     output:
//     tuple val(cluster_name), path("${cluster_name}_combined_fastree_lg_pseudo.tree"), emit: tree
//     path "${cluster_name}_fasttree.log", emit: log

//     script:
//     """
//     fastTreeMP="/home/groups/VEO/tools/fastTreeMP/v2.1.11/FastTreeMP"

//     \$fastTreeMP -lg -pseudo < ${msa} > ${cluster_name}_combined_fastree_lg_pseudo.tree 2> ${cluster_name}_fasttree.log
//     """
// }

// // =============================================================
// // 4. RF_MODEL
// // =============================================================
// process RF_MODEL {
//     cpus 30
//     memory '100 GB'
//     time '6h'

//     errorStrategy 'retry'
//     maxRetries 2

//     publishDir "${params.out_dir}/RF_results/raw", mode: 'copy'

//     queue 'standard,short,long,fat,gpu-veo'

//     input:
//     tuple val(cluster_name), path(msa_file)

//     output:
//     path "RF_results/*_Feature_importance.tsv", emit: features
//     path "RF_results/*_RandomForest_Report.tsv", emit: reports
//     path "RF_results/*_data_df.tsv", emit: data_df
//     path "${cluster_name}_rf_model.log", emit: log

//     script:
//     """
//     source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
//     conda activate working-env

//     python3 ${projectDir}/bin/rf_model_v3.py \\
//         ${msa_file} \\
//         ${params.metadata_dir}/${cluster_name}_metadata.tsv \\
//         ${cluster_name} >> ${cluster_name}_rf_model.log 2>&1
//     """
// }

// // =============================================================
// // 5. Create phyloglm dataset
// // =============================================================
// process CREATE_PHYLOGLM_DATASET {
//     publishDir "${params.out_dir}/test_phyloglm", mode: 'copy'
//     cpus 10
//     memory '30 GB'
//     time '1h'

//     input:
//     path feature
//     path data

//     output:
//     tuple val(cluster_name), path("${cluster_name}_phyloglm_dataset.tsv"), emit: dataset
//     tuple val(cluster_name), path("${cluster_name}_gini.tsv"), emit: gini
//     path "create_dataset.log", emit: log

//     script:
//     cluster_name = data.baseName.replaceAll(/_data_df.*/, '')
//     """
//     source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
//     conda activate working-env

//     python3 /vast/groups/VEO/shared_data/Lydia_Varada/bin/create_dataset_phyloglm.py \\
//         ${feature} ${data} -o ${cluster_name}_phyloglm_dataset.tsv >> create_dataset.log 2>&1
//     """
// }

// // =============================================================
// // 6. Split gini
// // =============================================================
// process SPLIT_GINI {
//     cpus 1
//     memory '5 GB'
//     time '30m'

//     input:
//     tuple val(cluster_name), path(gini_file)

//     output:
//     tuple val(cluster_name), path("${cluster_name}_chunk_*.tsv"), emit: chunks

//     script:
//     """
//     n_chunks=15
//     header=\$(head -1 ${gini_file})
//     tail -n +2 ${gini_file} | split -d -n l/\$n_chunks - tmp_chunk_
//     for f in tmp_chunk_*; do
//         idx=\${f#tmp_chunk_}
//         outname="${cluster_name}_chunk_\${idx}.tsv"
//         echo "\$header" > "\$outname"
//         cat "\$f" >> "\$outname"
//     done
//     """
// }

// // =============================================================
// // 7. phyloglm_fit_chunk
// // =============================================================
// process PHYLOGLM_FIT_CHUNK {
//     cpus 20
//     memory '80 GB'
//     time '4h'
//     maxForks 20

//     errorStrategy 'retry'
//     maxRetries 2

//     queue 'standard,short,long,fat'

//     input:
//     tuple val(cluster_name), path(dataset), path(tree), path(gini_chunk)

//     output:
//     tuple val(cluster_name), path("${gini_chunk.baseName}_partial.csv"), emit: partial
//     path "${gini_chunk.baseName}_fit.log", emit: log

//     script:
//     """
//     source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
//     conda activate phyloglm

//     Rscript /vast/groups/VEO/shared_data/Lydia_Varada/bin/phyloglm_fit.R \\
//         ${dataset} ${tree} ${gini_chunk} ${gini_chunk.baseName}_partial.csv >> ${gini_chunk.baseName}_fit.log 2>&1
//     """
// }

// // =============================================================
// // 8. merge_phyloglm
// // =============================================================
// process MERGE_PHYLOGLM {
//     publishDir "${params.out_dir}/test_phyloglm", mode: 'copy'
//     cpus 4
//     memory '20 GB'
//     time '1h'

//     input:
//     tuple val(cluster_name), path(partials)

//     output:
//     tuple val(cluster_name), path("${cluster_name}_phyloglm.csv"), emit: csv
//     path "${cluster_name}_phyloglm_plot.png", emit: plot
//     path "${cluster_name}_merge.log", emit: log

//     script:
//     """
//     source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
//     conda activate phyloglm

//     Rscript /vast/groups/VEO/shared_data/Lydia_Varada/bin/phyloglm_merge.R \\
//         ${partials} ${cluster_name}_phyloglm.csv ${cluster_name}_phyloglm_plot.png >> ${cluster_name}_merge.log 2>&1
//     """
// }

// // =============================================================
// // 9. plot_zvalue
// // =============================================================
// process PLOT_ZVALUE {
//     publishDir "${params.out_dir}/test_phyloglm", mode: 'copy'
//     cpus 2
//     memory '10 GB'
//     time '1h'

//     input:
//     tuple val(cluster_name), path(csv)

//     output:
//     path "${cluster_name}_zvalue_plot.png", emit: zvalue_plot
//     path "${cluster_name}_zvalue.log", emit: log

//     script:
//     """
//     source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh
//     conda activate phyloglm

//     Rscript /vast/groups/VEO/shared_data/Lydia_Varada/bin/plot_zvalue_featimps.R \\
//         ${csv} ${cluster_name}_zvalue_plot.png >> ${cluster_name}_zvalue.log 2>&1
//     """
// }

// =============================================================
// 10. Plots (global, across all clusters)
// =============================================================
// process PLOTS {
//     publishDir "${params.out_dir}/plots", mode: 'copy'
//     cpus 10
//     memory '30 GB'
//     time '2h'

//     errorStrategy 'retry'
//     maxRetries 2

//     queue 'standard,short,long,fat'

//     input:
//     path results_dir
//     path fasta_dir
//     val biome_1
//     val biome_2
//     val colour1
//     val colour2

//     output:
//     path "*_seqLogo.png", emit: seq_logo_plot
//     path "*_biome_feat_enrich_bars.png", emit: feature_enrichment_plot
//     path "*_biome_feat_enrich_bars_fullprot.png", emit: feature_full_protein_plot
//     path "*_AA_radars.png", emit: AA_radars_plot
//     path "*_AA_enrich_bars.png", emit: protein_analysis_plot
//     path "*_protein_analysis_subplots.png", emit: protein_analysis_subplots_png
//     path "*_protein_analysis_subplots.svg", emit: protein_analysis_subplots_svg
//     path "*_protein_scale_flexibility.png", emit: protein_scale_flexibility
//     path "*_protein_scale_GRAVY.png", emit: protein_scale_GRAVY
//     path "*_protein_scale_chargeatpH7.png", emit: protein_scale_chargeatpH7
//     path "Rf_model_report_plots.png", emit: report_plot_png
//     path "Rf_model_report_plots.svg", emit: report_plot_svg
//     path "Rf_model_precision_v_recall.png", emit: precision_recall
//     path "plots.log", emit: log

//     script:
//     """
//     source /vast/groups/VEO/shared_data/Lydia_Varada/miniconda3/etc/profile.d/conda.sh
//     conda activate working-env

//     results_path=\$(realpath "${results_dir}")
//     fasta_path=\$(realpath "${fasta_dir}")

//     python3 /vast/groups/VEO/shared_data/Lydia_Varada/bin/plots/plotting_script.py \\
//         "\$results_path" "\$fasta_path" "${biome_1}" "${biome_2}" "${colour1}" "${colour2}" 50 \$PWD >> plots.log 2>&1
//     """
// }

// =============================================================
// Workflow
// =============================================================
workflow {
    //Entry point: genomes
    if (params.entry_point == "genome") {
        input_ch = Channel.fromPath("${params.fasta_file}/*.fna")
        if (params.gene_prediction == "phanotate") {
            prediction_ch = PHANOTATE(input_ch)
        } else {
            prediction_ch = PRODIGAL(input_ch)
        }

        //prediction_ch = Channel.fromPath("/vast/groups/VEO/shared_data/Lydia_Varada/first_part_of_pipe/phanotate/phanotate.faa")
        //MMseq2
        mmseq2_ch = MMSEQ2(prediction_ch.fasta, params.id, params.cov)
        
        //mmseq2_ch = Channel.fromPath("/vast/groups/VEO/shared_data/Lydia_Varada/first_part_of_pipe/mmseq2/mmseq2_cluster.tsv")
        //prediction_ch = Channel.fromPath("/vast/groups/VEO/shared_data/Lydia_Varada/first_part_of_pipe/prodigal/proteins_prodigal.faa")
        metadata = Channel.fromPath("/vast/groups/VEO/shared_data/Lydia_Varada/first_part_of_pipe/hic_viral_target.csv")
        
        orthogroups_ch = ORTHOGROUPS(mmseq2_ch.tsv , prediction_ch.fasta , metadata)

    } else if (params.entry_point == "clusters") {
        // Entry point: one fasta file per cluster/group
        // fasta_ch = Channel.fromPath("${params.fasta_dir}/*.faa")
        //     .map { f -> tuple(f.baseName, f) }

    } else {
        error "Please specify --entry_point as either 'genome' or 'clusters'"
    }

    // // 1. MSA
    // mafft_ch = MAFFT(fasta_ch).aligned

    // // 2. Optional trimAl
    // msa_ch = params.do_trimal ? TRIMAL(mafft_ch).faa : mafft_ch

    // // 3. FastTree — msa_ch can be reused directly, no into() needed
    // fasttree_ch = FASTTREE(msa_ch).tree

    // // 4. RF_MODEL — msa_ch reused again
    // rf_ch = RF_MODEL(msa_ch)

    // // 5. Create phyloglm dataset — rf_ch.data_df reused below too
    // dataset_ch = CREATE_PHYLOGLM_DATASET(rf_ch.features, rf_ch.data_df)

    // // 6. Split gini
    // chunks_ch = SPLIT_GINI(dataset_ch.gini)
    //     .transpose()

    // fit_input_ch = chunks_ch
    //     .combine(dataset_ch.dataset, by: 0)
    //     .combine(fasttree_ch, by: 0)
    //     .map { name, chunk, dataset, tree -> tuple(name, dataset, tree, chunk) }

    // // 7. phyloglm_fit_chunk
    // partial_ch = PHYLOGLM_FIT_CHUNK(fit_input_ch)

    // grouped_partials_ch = partial_ch.partial.groupTuple()

    // // 8. merge_phyloglm
    // merged_ch = MERGE_PHYLOGLM(grouped_partials_ch)

    // // 9. plot_zvalue
    // zvalue_ch = PLOT_ZVALUE(merged_ch.csv)

    // // 10. Plots
    // biome1_short = params.biome1.tokenize(';').last()
    // biome2_short = params.biome2.tokenize(';').last()

    // results_dir_ch = rf_ch.data_df.collect().map { file("${params.out_dir}/RF_results/raw") }
    // fasta_dir_ch   = Channel.value(file(params.fasta_dir))

    // plots_ch = PLOTS(results_dir_ch, fasta_dir_ch, biome1_short, biome2_short, params.colour1, params.colour2)
}