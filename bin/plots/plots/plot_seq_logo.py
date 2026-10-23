import pandas as pd
import logomaker
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import seaborn as sns

# def make_matrices(cluster_num, biome1, biome2, x=50):

#     # read in the data_df for the cluster
#     data_df_merged_balanced = pd.read_csv(f"RF_results/raw/{cluster_num}_data_df.tsv", sep='\t')
#     data_df_merged_balanced['biome'] = data_df_merged_balanced['biome'].str.replace(":",";")

#     temp = data_df_merged_balanced.groupby("ecosystem_subtype").sum().reset_index()
#     temp_1 = temp.drop(columns=['protein','biome','ecosystem','ecosystem_category','ecosystem_type','specific_ecosystem'])
#     temp_1_melted = temp_1.melt(id_vars="ecosystem_subtype", var_name="pos_aa", value_name="count")
#     temp_1_melted[["position", "aa"]] = temp_1_melted['pos_aa'].str.split("_",expand=True)
#     temp_1_melted['position'] = temp_1_melted['position'].astype(int)

#     # create separate matrices for each 
#     ecosystem_matrices = {
#         ecosys: subdf.pivot(index="position", columns="aa", values="count").fillna(0).astype(int)
#         for ecosys, subdf in temp_1_melted.groupby("ecosystem_subtype")
#     }

#     # Example: access the 'Lake' matrix
#     biome1_matrix = ecosystem_matrices[biome1]

#     # Example: access the 'Oceanic' matrix
#     biome2_matrix = ecosystem_matrices[biome2]

#     # Read in feature importance scores
#     feat_imp_df = pd.read_csv(f"RF_results/raw/{cluster_num}_Feature_importance.tsv", sep='\t', usecols=[0,2],names=['col_name','feature_importance_vals'],header=0)
#     top_50_ft_df = feat_imp_df.head(n=x)

#     top_50_ft_df[['position','aa']] = top_50_ft_df['col_name'].str.split("_",expand=True)
#     top_50_ft_df['position'] = top_50_ft_df['position'].astype(int)

#     biome1_matrix_filtered = biome1_matrix.loc[biome1_matrix.index.intersection(top_50_ft_df['position'])]
#     print(biome1_matrix_filtered)
#     biome2_matrix_filtered = biome2_matrix.loc[biome2_matrix.index.intersection(top_50_ft_df['position'])]

#     data_df_merged_balanced_1 = data_df_merged_balanced.groupby('ecosystem_subtype').apply(lambda x: (x==1).sum()).T
#     data_df_merged_balanced_1
#     data_df_merged_balanced_1[f"Proportion_{biome1}"] = (
#         data_df_merged_balanced_1[biome1] / data_df_merged_balanced[data_df_merged_balanced['ecosystem_subtype'] == biome1].shape[0] * 100)
#     data_df_merged_balanced_1[f"Proportion_{biome2}"] = (
#         data_df_merged_balanced_1[biome2] / data_df_merged_balanced[data_df_merged_balanced['ecosystem_subtype'] == biome2].shape[0] * 100)

#     # print(data_df_merged_balanced_1)
#     data_df_merged_balanced_1 = data_df_merged_balanced_1.T

#     top_50_ft_df_t = top_50_ft_df.set_index('col_name').T
#     top_features = data_df_merged_balanced_1[top_50_ft_df_t.columns]
#     # print(top_features.T.head(n=20))

#     # only proportions
#     top_features_melted = top_features.reset_index().melt(id_vars="ecosystem_subtype")
#     top_features_melted_1 = top_features_melted.sort_values(['variable','value'],ascending=[True,False]).groupby('variable').head(n=1)
#     top_feat_bars = top_features_melted_1.merge(top_50_ft_df,left_on='variable',right_on='col_name')
#     print(top_feat_bars)
    
#     return biome1_matrix_filtered, biome2_matrix_filtered, top_feat_bars

def make_matrices(path_to_rf_results,cluster_num, biome1, biome2, x=50):

    # read in the data_df for the cluster
    data_file = path_to_rf_results.joinpath(f"{cluster_num}_data_df.tsv")
    data_df_merged_balanced = pd.read_csv(data_file, sep='\t')
    data_df_merged_balanced['biome'] = data_df_merged_balanced['biome'].str.replace(":", ";")

    # count amino acids per ecosystem_subtype
    temp = data_df_merged_balanced.groupby("ecosystem_subtype").sum().reset_index()
    temp_1 = temp.drop(columns=['protein','biome','ecosystem',
                                'ecosystem_category','ecosystem_type','specific_ecosystem'])
    temp_1_melted = temp_1.melt(id_vars="ecosystem_subtype",
                                var_name="pos_aa", value_name="count")
    temp_1_melted[["position", "aa"]] = temp_1_melted['pos_aa'].str.split("_", expand=True)
    temp_1_melted['position'] = temp_1_melted['position'].astype(int)

    # create matrices for each ecosystem
    ecosystem_matrices = {
        ecosys: subdf.pivot(index="position", columns="aa", values="count").fillna(0).astype(int)
        for ecosys, subdf in temp_1_melted.groupby("ecosystem_subtype")
    }

    # just the 20 amino acids without (X, B, Z, *, -, ...)
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    ecosystem_matrices = {
        ecosys: matrix[[c for c in matrix.columns if c in valid_aa]]
        for ecosys, matrix in ecosystem_matrices.items()
    }

    biome1_matrix = ecosystem_matrices[biome1]
    biome2_matrix = ecosystem_matrices[biome2]

    # feature importance scores
    feat_file = path_to_rf_results.joinpath(f"{cluster_num}_Feature_importance.tsv")
    feat_imp_df = pd.read_csv(
        feat_file,
        sep='\t', usecols=[0,2],
        names=['col_name','feature_importance_vals'], header=0
    )
    top_50_ft_df = feat_imp_df.head(n=x)
    top_50_ft_df[['position','aa']] = top_50_ft_df['col_name'].str.split("_", expand=True)
    top_50_ft_df['position'] = top_50_ft_df['position'].astype(int)

    # filter matrices to only top feature positions
    biome1_matrix_filtered = biome1_matrix.loc[
        biome1_matrix.index.intersection(top_50_ft_df['position'])
    ]
    biome2_matrix_filtered = biome2_matrix.loc[
        biome2_matrix.index.intersection(top_50_ft_df['position'])
    ]

    # calculate proportions per biome (normalize counts by number of rows per biome)
    counts_per_biome = data_df_merged_balanced.groupby('ecosystem_subtype').apply(lambda x: (x==1).sum()).T

    proportion_df = pd.DataFrame()
    proportion_df[biome1] = (
        counts_per_biome[biome1] /
        data_df_merged_balanced[data_df_merged_balanced['ecosystem_subtype'] == biome1].shape[0] * 100
    )
    proportion_df[biome2] = (
        counts_per_biome[biome2] /
        data_df_merged_balanced[data_df_merged_balanced['ecosystem_subtype'] == biome2].shape[0] * 100
    )
    proportion_df = proportion_df.T  # biome names stay clean here

    # restrict to top features
    top_50_ft_df_t = top_50_ft_df.set_index('col_name').T
    top_features = proportion_df[top_50_ft_df_t.columns]

    # melt into long form and merge with feature importance
    top_features_melted = (
    top_features
    .reset_index()
    .rename(columns={'index': 'ecosystem_subtype'})
    .melt(id_vars="ecosystem_subtype")  
    )

    top_features_melted_1 = (
        top_features_melted
        .sort_values(['variable','value'], ascending=[True, False])
        .groupby('variable').head(n=1)
    )
    top_feat_bars = top_features_melted_1.merge(
        top_50_ft_df, left_on='variable', right_on='col_name'
    )
    print(top_feat_bars)
    return biome1_matrix_filtered, biome2_matrix_filtered, top_feat_bars


def plot_seq_logo(biome1, biome2, biome1_matrix_filtered, biome2_matrix_filtered, top_feat_bars,colour1, colour2):

    # Create the new figure and subplots
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(30, 5), dpi = 300)

    # Create the marine logo directly on the first subplot
    try:
        AA_logo_biome2 = logomaker.Logo(
            biome2_matrix_filtered,
            ax=axes[1],  # Use the first subplot's axes
            font_name='Arial Rounded MT Bold',
            color_scheme='chemistry',#'dmslogo_charge'
            stack_order='small_on_top',  # Ensures visibility of larger letters
            width=2.5
        )
    except Exception as e:
        print("Fehlertyp:", type(e))
        print("Fehlerinhalt:", repr(e))
        import traceback
        traceback.print_exc()
        raise

    # Create the freshwater logo directly on the second subplot
    AA_logo_biome1 = logomaker.Logo(
        biome1_matrix_filtered*-1,
        ax=axes[2],  # Use the second subplot's axes
        font_name='Arial Rounded MT Bold',
        color_scheme='dmslogo_charge',
        flip_below=False,
        baseline_width = 10,
        stack_order='small_on_top',  # Ensures visibility of larger letters
        width=2.5
    )

    # create marine bars above first seqlogo (axes[0])
    top_feat_bars = top_feat_bars[top_feat_bars["position"].isin(biome2_matrix_filtered.index)]
    data_m=top_feat_bars.loc[top_feat_bars['ecosystem_subtype'] == biome2]

    axes[0].bar(
        data_m["position"], 
        data_m["feature_importance_vals"], 
        color=colour2,
        width=2  # Adjust width to match sequence logo spacing
    )

    # create freshwater bars below last seqlogo (axes[3])
    top_feat_bars = top_feat_bars[top_feat_bars["position"].isin(biome1_matrix_filtered.index)]
    data_f = top_feat_bars.loc[top_feat_bars['ecosystem_subtype'] == biome1]

    axes[3].bar(
        data_f["position"], 
        data_f["feature_importance_vals"], 
        color=colour1,
        width=2
    )
    axes[3].invert_yaxis()

    # for ax in axes:
    #     ax.set_xticks([0, 500])

    # Ensure Logomaker uses the same ticks
    AA_logo_biome2.ax.set_xticks(biome1_matrix_filtered.index)
    AA_logo_biome1.ax.set_xticks(biome1_matrix_filtered.index)

    axes[-1].set_xticklabels(biome1_matrix_filtered.index, rotation=90, fontsize=10)

    axes[0].tick_params(labelsize=10)
    axes[1].tick_params(labelsize=10)
    axes[2].tick_params(labelsize=10)
    axes[3].tick_params(labelsize=10)


    # Adjust layout
    plt.tight_layout()
    # Display the combined figure
    # plt.show()

    return plt

def plot(path_to_rf_results,cluster_num,biome1, biome2, colour1, colour2,x=50):
    """ takes variables from outside script here sends them to internal functions"""
    biome1_matrix_filtered, biome2_matrix_filtered, top_feat_bars = make_matrices(path_to_rf_results,cluster_num, biome1, biome2, x)
    plt = plot_seq_logo(biome1, biome2, biome1_matrix_filtered, biome2_matrix_filtered, top_feat_bars,colour1, colour2)
    print("subscript done")
    return plt