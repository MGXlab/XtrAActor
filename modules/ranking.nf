process RANK_CLUSTERS {
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path(reports, stageAs: 'reports/*')
    path(fastas,  stageAs: 'cluster_fastas/*')

    output:
    path "cluster_ranking.tsv", emit: ranking
    path "rank_clusters.log",   emit: log

    script:
    """
    source ${params.conda_sh_ri}
    conda activate working-env

    python3 ${projectDir}/bin/rank_clusters.py reports cluster_fastas \\
        -o cluster_ranking.tsv 2> rank_clusters.log
    """
}