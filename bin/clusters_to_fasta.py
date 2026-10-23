import pandas as pd
import numpy as np
from collections import defaultdict
from Bio import SeqIO
import sys
from pathlib import Path
import argparse

def get_inputs():
    parser = argparse.ArgumentParser(description="This script will take the clusters.tsv file, parse the clusters, add the biome and taxonomy metadata. " \
    "Based on which taxonomy and biomes you want to compare, it will filter the clusters and produce fasta files per cluster to turn into a MSA" \
    "If the cluster has more than 30 members for each specified biome, it creates its own cluster" \
    "If the cluster has between 5 and 30 members for either specified biome, it goes into a mixed cluster where the members are put together" \
    "Clusters with <5 members for either biome are not considered for further analyses" \
    "Each cluster should not have a disproportionate number of members from one biome, if biome1 = 3*biome2 (or vise versa), the cluster is not considered")

    parser.add_argument("-clusters_file","--clusters_file",
        help="Path to clusters.tsv file from foldseek"
    )
    parser.add_argument("-fasta_file","--fasta_file", 
        help="Path to fasta file that contains ALL sequences"
    )
    parser.add_argument("-biomes_file","--biomes_file", 
        help="Path to biomes.csv file, e.g. 'best_biome_per_mcp_protein.csv'"
    )
    parser.add_argument("-taxa_file","--taxa_file", 
        help="Path to taxonomy.tsv file, e.g. 'clean_taxonomy_per_protein.tsv'"
    )
    parser.add_argument(
        "--num_seqs", type=int, default=25,
        help="Number of members you want per cluster as a first filter (default: 25)"
    )
    parser.add_argument(
        "--min_seqs_own_model", type=int, default=30,
        help="minimum number of sequences per biome required for a cluster to be its own model"
    )
    parser.add_argument(
        "-o", "--out_dir", required=True,
        help="Output directory where the cluster FASTAs will go. "
             "Also produces 'list_of_paths_to_cluster_fasta.txt' in the root folder."
    )
    parser.add_argument(
        "-t", "--wanted_taxonomy", required=True, default=None,
        help='Filter by taxonomy (e.g. "Caudoviricetes" or "Microviridae")'
    )
    parser.add_argument(
        "-b1", "--biome1", required=True, default=None,
        help='Primary biome filter (e.g. "Environmental;Aquatic;Freshwater;Lake")'
    )
    parser.add_argument(
        "-b2", "--biome2", required=True, default=None,
        help='Secondary biome filter (e.g. "Environmental;Aquatic;Thermal springs;Hot (42-90C)")'
    )
    return parser.parse_args()


def clean_clusters_file(path_to_clusters,num_seqs=25):
    """ Clean a clusters.tsv file and filter clusters by minimum size """

    df = pd.read_csv(
        path_to_clusters,
        sep='\t',
        header=None,
        names=['cluster_rep', 'member']
    )

    # --- Cleaning functions ---
    '''Change the values with "LAB|" to "LAB_" '''
    def replace_first_two_pipes(val):
        return val.replace('|', '_', 1).replace('|', '_', 1)

    def replace_colons(val):
        return val.replace(':', '_', 1)

     # --- Apply cleaning ---
    for col in ['cluster_rep', 'member']:
        df[col] = df[col].where(~df[col].str.contains('LAB'), df[col].map(replace_first_two_pipes))
        df[col] = df[col].str.split('|', expand=True)[0]
        df[col] = df[col].where(~df[col].str.contains('SULFURCAVE'), df[col].map(replace_colons))

    # --- Count cluster sizes ---
    cluster_sizes = df.groupby('cluster_rep').size().sort_values(ascending=False)
    print(cluster_sizes)

    # --- Filter by minimum size ---
    valid_clusters = cluster_sizes[cluster_sizes >= num_seqs].index

    df_filtered = df[df['cluster_rep'].isin(valid_clusters)]

    return df_filtered

def convert_dataframe_to_dict(df_with_biome_taxa, biome1, biome2):
    
    # most common taxonomy
    df_mode_taxonomy = (
        df_with_biome_taxa.groupby('cluster_rep')['taxonomy']
        .agg(lambda x: pd.Series.mode(x).iloc[0] if not pd.Series.mode(x).empty else np.nan)
        .reset_index(name='mode_CR_taxonomy')
    )
    # rename columns
    df_mode_taxonomy.rename(columns = {'taxonomy' : 'mode_CR_taxonomy'},inplace=True)
    df_mode_taxonomy

    # convert to a dictionary
    mode_taxa_dict = {
        cluster_rep : mode_taxonomy
        for cluster_rep, mode_taxonomy in zip(df_mode_taxonomy['cluster_rep'],df_mode_taxonomy['mode_CR_taxonomy'])
    }
    # get biomes in each cluster and their numbers
    df_biomes_per_cluster = df_with_biome_taxa.groupby(['cluster_rep','biome']).count().reset_index()
    df_biomes_per_cluster.drop(columns=['ecosystem','ecosystem_category','ecosystem_type','ecosystem_subtype','specific_ecosystem'],inplace=True)
    df_biomes_per_cluster.rename(columns={'taxonomy' : 'biome_count'}, inplace=True)
    df_biomes_per_cluster

    biome_dict = defaultdict(dict)
    for _,row in df_biomes_per_cluster.iterrows():
        cluster = row['cluster_rep']
        biome = row['biome']
        count = row['biome_count']
        biome_dict[cluster][biome] = count
    
    # add mode taxonomy into biome_dict for querying
    for cluster_rep in biome_dict:
        if cluster_rep in mode_taxa_dict:
            biome_dict[cluster_rep]['mode_taxonomy'] = mode_taxa_dict[cluster_rep]

    return biome_dict

def query_dictionary(biome_dict, wanted_taxonomy, biome1, biome2, min_seqs_own_model):
    own_models_list = []
    mixed_models_list = []
    print(biome1)
    print(biome2)
    biome1_short = biome1.split(';')[3]
    biome2_short = biome2.split(';')[3]
    for cluster_rep, biomes in biome_dict.items():
        mode_taxonomy = biomes.get('mode_taxonomy') # skips any that are missing a mode_taxonomy
        if wanted_taxonomy in str(mode_taxonomy):
            if biome1 in biome_dict[cluster_rep] and biome2 in biome_dict[cluster_rep]:
                biome1_members = biome_dict[cluster_rep][biome1]
                biome2_members = biome_dict[cluster_rep][biome2]
                if (biome1_members > biome2_members and biome1_members < 100*biome2_members) or (biome2_members > biome1_members and biome2_members < 100*biome1_members):
                    if biome1_members >= min_seqs_own_model and biome2_members >= min_seqs_own_model:
                        own_models_list.append(cluster_rep)
                        print(f"own model {len(own_models_list)} -> {cluster_rep}, {biome1_short} members:{biome1_members}, {biome2_short} members:{biome2_members}")
                    else:
                    # elif 5 <= biome1_members < min_seqs_own_model and 5 <= biome2_members < min_seqs_own_model:
                        print(f"mixed model -> {cluster_rep}, {biome1_short} members:{biome1_members}, {biome2_short} members:{biome2_members}")
                        mixed_models_list.append(cluster_rep)
    return own_models_list, mixed_models_list

def get_model_members(df_with_biome_taxa, own_models_list, mixed_models_list):
    # getting own models members into a dictionary
    member_dict = {}
    i=0
    for rep in own_models_list:
        # print(rep)
        df_list = df_with_biome_taxa.loc[df_with_biome_taxa["cluster_rep"] == str(rep)]
        temp_list = df_list['member'].tolist()
        member_dict[i] = temp_list
        if i == 2:
            print(member_dict[i])
        i+=1

    df_list = df_with_biome_taxa[df_with_biome_taxa['cluster_rep'].isin(mixed_models_list)]
    # print(df_list)
    temp_list = df_list['member'].tolist()
    member_dict['mixed_model'] = temp_list

    return member_dict

def get_sequences(member_dict,fasta_file,out_dir):
    # Reverse lookup: seqid -> cluster
    seq_to_cluster = {
        seqid: cluster for cluster, seqs in member_dict.items() for seqid in seqs}
    
    cluster_records = {cluster: [] for cluster in member_dict}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) >= 400:
            seqid = record.id.split("|")[0]
            if seqid in seq_to_cluster:
                cluster_num = seq_to_cluster[seqid]
                #print(len(seq_to_cluster))
                if cluster_num is not None:
                    cluster_records[cluster_num].append(record)     # add records to the cluster_records dict

    # writing one fasta per cluster in one go
    for cluster_num, records in cluster_records.items():
        if len(records) >= 10:
            out_filename = f"{out_dir}/cluster_{cluster_num}.faa"
            with open(out_filename, "w") as out_file:
                SeqIO.write(records, out_file, "fasta")
        else:
            print(f"SKIPPED cluster {cluster_num}: only {len(records)} records after length filtering (need >=10)")
    return


if __name__ == '__main__':
    args=get_inputs()

    clusters_file = args.clusters_file
    num_seqs = args.num_seqs
    min_seqs_own_model = args.min_seqs_own_model    # minimum number of sequences per biome required for a cluster to be its own model
    biomes_file=args.biomes_file
    taxa_file=args.taxa_file
    out_dir=args.out_dir
    wanted_taxonomy = args.wanted_taxonomy   #"Caudoviricetes"
    biome1 = args.biome1    #"Environmental;Aquatic;Freshwater;Lake"
    biome2 = args.biome2    #"Environmental;Aquatic;Thermal springs;Hot (42-90C)"
    fasta_file = args.fasta_file


    # "./clusters_50_50.tsv"
    # get filtered clusters.tsv file. Only clusters with more members than specified in "num_seqs" are considered
    df_filtered = clean_clusters_file(clusters_file, num_seqs)

    # read biomes in
    # "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/amino_acid_analyis/imgvr_mgnify_combined/best_biome_per_mcp_protein.csv"

    df_biomes = pd.read_csv(biomes_file)
    df_biomes.drop('cluster_rep',axis=1,inplace=True)
    df_biomes['biome'] = df_biomes['biome'].str.replace(':', ";")
    print(df_biomes)

    # read taxonomy in
    # "/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/amino_acid_analyis/imgvr_mgnify_combined/clean_taxonomy_per_protein.tsv"

    df_taxonomy = pd.read_csv(taxa_file, sep='\t',header=0)
    print(df_taxonomy)

    # add biomes to the clusters.tsv file
    df_temp1 = df_filtered.set_index('member').join(df_biomes.set_index('protein'))
    df_with_biome_taxa = df_temp1.join(df_taxonomy.set_index('protein')).reset_index()

    # --- filter the df_with_biome_taxa for only members that are either biome1 or biome2 or that the member is a lab phage --- 
    # --- these will be the only ones included in the final cluster_x.faa ---
    df_with_biome_taxa = df_with_biome_taxa.loc[(df_with_biome_taxa['biome'] == biome1) | 
                                                (df_with_biome_taxa['biome'] == biome2) | 
                                                (df_with_biome_taxa['member'].str.contains("LAB"))]

    biome_dict = convert_dataframe_to_dict(df_with_biome_taxa, biome1, biome2)
    own_models_list, mixed_models_list = query_dictionary(biome_dict, wanted_taxonomy, biome1, biome2, min_seqs_own_model)
    member_dict = get_model_members(df_with_biome_taxa, own_models_list, mixed_models_list)
    
    Path(out_dir).mkdir(parents=True, exist_ok=True)    # make the output directory if it doesn't exist
    get_sequences(member_dict,fasta_file,out_dir)
    print("UPDATE: writing fasta files for each cluster complete")

    with open(f"{out_dir}/list_of_paths_to_cluster_fasta.txt",'a') as fh:
        for path in Path(out_dir).glob('**/*.faa'):
            print(path)
            fh.write(f"{path}\n")
    print("UPDATE: writing txt file of fasta files complete \n" \
        "UPDATE: Cluster selection script complete")
    