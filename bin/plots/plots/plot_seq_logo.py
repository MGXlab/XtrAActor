import pandas as pd
import logomaker
import numpy as np
import matplotlib.pyplot as plt

# Expected inputs (from rf_model_v3.py, in path_to_rf_results):
#   <cluster>_data_df.tsv            : gene | <pos_aa feature columns, 0/1> | biome (full biome string)
#   <cluster>_Feature_importance.tsv : feature | gini_imp | cluster_num


def make_matrices(path_to_rf_results, cluster_num, biome1, biome2, x=50):

    # read in the data_df for the cluster
    data_file = path_to_rf_results.joinpath(f"{cluster_num}_data_df.tsv")
    data_df_merged_balanced = pd.read_csv(data_file, sep='\t')
    data_df_merged_balanced['biome'] = data_df_merged_balanced['biome'].str.replace(":", ";")

    # feature columns = everything except the gene ID and the label
    feature_cols = [c for c in data_df_merged_balanced.columns if c not in ('gene', 'biome')]

    for b in (biome1, biome2):
        if b not in set(data_df_merged_balanced['biome']):
            raise ValueError(f"Biome '{b}' not found in {data_file.name}. "
                             f"Biomes present: {sorted(data_df_merged_balanced['biome'].unique())}")

    # count amino acids per biome (only the feature columns)
    temp = data_df_merged_balanced.groupby("biome")[feature_cols].sum().reset_index()
    temp_melted = temp.melt(id_vars="biome", var_name="pos_aa", value_name="count")
    temp_melted[["position", "aa"]] = temp_melted['pos_aa'].str.split("_", n=1, expand=True)
    temp_melted['position'] = temp_melted['position'].astype(int)

    # create matrices for each biome: rows = alignment position, columns = amino acid
    ecosystem_matrices = {
        ecosys: subdf.pivot(index="position", columns="aa", values="count").fillna(0).astype(int)
        for ecosys, subdf in temp_melted.groupby("biome")
    }

    # just the 20 amino acids without (X, B, Z, *, -, ...)
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    ecosystem_matrices = {
        ecosys: matrix[[c for c in matrix.columns if c in valid_aa]]
        for ecosys, matrix in ecosystem_matrices.items()
    }

    biome1_matrix = ecosystem_matrices[biome1]
    biome2_matrix = ecosystem_matrices[biome2]

    # feature importance scores: column 0 = feature, column 1 = gini_imp
    feat_file = path_to_rf_results.joinpath(f"{cluster_num}_Feature_importance.tsv")
    feat_imp_df = pd.read_csv(
        feat_file,
        sep='\t', usecols=[0, 1],
        names=['col_name', 'feature_importance_vals'], header=0
    )
    top_50_ft_df = feat_imp_df.head(n=x).copy()
    top_50_ft_df[['position', 'aa']] = top_50_ft_df['col_name'].str.split("_", n=1, expand=True)
    top_50_ft_df['position'] = top_50_ft_df['position'].astype(int)

    # filter matrices to only top feature positions
    biome1_matrix_filtered = biome1_matrix.loc[
        biome1_matrix.index.intersection(top_50_ft_df['position'])
    ]
    biome2_matrix_filtered = biome2_matrix.loc[
        biome2_matrix.index.intersection(top_50_ft_df['position'])
    ]

    # calculate proportions per biome (normalize counts by number of rows per biome)
    counts_per_biome = data_df_merged_balanced.groupby('biome')[feature_cols].apply(lambda d: (d == 1).sum()).T

    proportion_df = pd.DataFrame()
    proportion_df[biome1] = (
        counts_per_biome[biome1] /
        (data_df_merged_balanced['biome'] == biome1).sum() * 100
    )
    proportion_df[biome2] = (
        counts_per_biome[biome2] /
        (data_df_merged_balanced['biome'] == biome2).sum() * 100
    )
    proportion_df = proportion_df.T  # rows = biomes, columns = features

    # restrict to top features (only those present as columns)
    top_cols = [c for c in top_50_ft_df['col_name'] if c in proportion_df.columns]
    top_features = proportion_df[top_cols]

    # melt into long form and merge with feature importance
    top_features_melted = (
        top_features
        .reset_index()
        .rename(columns={'index': 'biome'})
        .melt(id_vars="biome")
    )

    # for each feature: the biome where it is more frequent
    top_features_melted_1 = (
        top_features_melted
        .sort_values(['variable', 'value'], ascending=[True, False])
        .groupby('variable').head(n=1)
    )
    top_feat_bars = top_features_melted_1.merge(
        top_50_ft_df, left_on='variable', right_on='col_name'
    )
    print(top_feat_bars)
    return biome1_matrix_filtered, biome2_matrix_filtered, top_feat_bars


def plot_seq_logo(biome1, biome2, biome1_matrix_filtered, biome2_matrix_filtered, top_feat_bars, colour1, colour2):

    # Create the new figure and subplots
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(30, 5), dpi=300)

    # biome2 logo on the second subplot
    AA_logo_biome2 = logomaker.Logo(
        biome2_matrix_filtered,
        ax=axes[1],
        font_name='Arial Rounded MT Bold',
        color_scheme='chemistry',
        stack_order='small_on_top',  # Ensures visibility of larger letters
        width=2.5
    )

    # biome1 logo (mirrored) on the third subplot
    AA_logo_biome1 = logomaker.Logo(
        biome1_matrix_filtered * -1,
        ax=axes[2],
        font_name='Arial Rounded MT Bold',
        color_scheme='dmslogo_charge',
        flip_below=False,
        baseline_width=10,
        stack_order='small_on_top',
        width=2.5
    )

    # biome2 importance bars above the logos (axes[0])
    bars_b2 = top_feat_bars[top_feat_bars["position"].isin(biome2_matrix_filtered.index)]
    data_m = bars_b2.loc[bars_b2['biome'] == biome2]

    axes[0].bar(
        data_m["position"],
        data_m["feature_importance_vals"],
        color=colour2,
        width=2  # Adjust width to match sequence logo spacing
    )

    # biome1 importance bars below the logos (axes[3])
    bars_b1 = top_feat_bars[top_feat_bars["position"].isin(biome1_matrix_filtered.index)]
    data_f = bars_b1.loc[bars_b1['biome'] == biome1]

    axes[3].bar(
        data_f["position"],
        data_f["feature_importance_vals"],
        color=colour1,
        width=2
    )
    axes[3].invert_yaxis()

    # Ensure Logomaker uses the same ticks
    AA_logo_biome2.ax.set_xticks(biome1_matrix_filtered.index)
    AA_logo_biome1.ax.set_xticks(biome1_matrix_filtered.index)

    axes[-1].set_xticklabels(biome1_matrix_filtered.index, rotation=90, fontsize=10)

    for ax in axes:
        ax.tick_params(labelsize=10)

    plt.tight_layout()
    return plt


def plot(path_to_rf_results, cluster_num, biome1, biome2, colour1, colour2, x=50):
    """ takes variables from outside script here sends them to internal functions"""
    biome1_matrix_filtered, biome2_matrix_filtered, top_feat_bars = make_matrices(
        path_to_rf_results, cluster_num, biome1, biome2, x)
    plt = plot_seq_logo(biome1, biome2, biome1_matrix_filtered, biome2_matrix_filtered,
                        top_feat_bars, colour1, colour2)
    print("subscript done")
    return plt