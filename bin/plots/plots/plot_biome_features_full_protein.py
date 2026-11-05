import pandas as pd
import numpy as np
import matplotlib as mpl
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

    # feature importance scores: column 0 = feature, column 1 = gini_imp
    feat_file = path_to_rf_results.joinpath(f"{cluster_num}_Feature_importance.tsv")
    feat_imp_df = pd.read_csv(feat_file, sep='\t', usecols=[0, 1],
                              names=['col_name', 'feature_importance_vals'],
                              header=0).set_index('col_name')

    # number of sequences carrying each feature, per biome (rows = features, columns = biomes)
    counts = data_df_merged_balanced.groupby('biome')[feature_cols].apply(lambda d: (d == 1).sum()).T

    # proportion (%) of sequences per biome carrying each feature; biome1 negative -> bars point down
    n_biome1 = (data_df_merged_balanced['biome'] == biome1).sum()
    n_biome2 = (data_df_merged_balanced['biome'] == biome2).sum()
    proportions = pd.DataFrame({
        f"Proportion_{biome1}": counts[biome1] / n_biome1 * 100 * -1,
        f"Proportion_{biome2}": counts[biome2] / n_biome2 * 100,
    })

    # add feature importance and reshape to long form
    temp = proportions.join(feat_imp_df, how='left').reset_index(names='pos_AA')
    temp_m = pd.melt(temp, id_vars=['pos_AA', 'feature_importance_vals'],
                     var_name='biome', value_name='proportion').copy()
    temp_m[['position', 'AA']] = temp_m['pos_AA'].str.split('_', n=1, expand=True)

    # per position and biome: keep only the most important feature
    df_final = (temp_m.sort_values("feature_importance_vals", ascending=False)
                      .drop_duplicates(['position', 'biome'])
                      .copy())

    imp_min = df_final['feature_importance_vals'].min()
    imp_max = df_final['feature_importance_vals'].max()
    df_final['feat_imp_normal'] = ((df_final['feature_importance_vals'] - imp_min)
                                   / (imp_max - imp_min + 1E-12)).astype(float)

    df_final = df_final.dropna(subset=["position", "proportion"])
    df_final["proportion"] = df_final["proportion"].astype(float)
    df_final["position"] = df_final["position"].astype(int)

    return df_final


def plot_bars(df_final, cluster_num, biome1, biome2, colour1, colour2, k):
    # colormaps from pale to the biome colour; darker = more important feature
    cmap_b1 = mpl.colors.LinearSegmentedColormap.from_list("biome1_cmap", ["#fafaf0", colour1])
    cmap_b2 = mpl.colors.LinearSegmentedColormap.from_list("biome2_cmap", ["#fafaf0", colour2])

    # normalize feature importance to the colour range
    norm = mpl.colors.Normalize(vmin=df_final["feature_importance_vals"].min(),
                                vmax=df_final["feature_importance_vals"].max())

    # separate data by biome
    b1 = df_final[df_final["biome"] == f"Proportion_{biome1}"].sort_values("position")
    b2 = df_final[df_final["biome"] == f"Proportion_{biome2}"].sort_values("position")

    fig, ax = plt.subplots(figsize=(30, 6))

    # biome1 bars (negative proportions, pointing down)
    ax.bar(
        b1["position"],
        b1["proportion"],
        color=cmap_b1(norm(b1["feature_importance_vals"])),
        width=0.5,
        label=biome1
    )

    # biome2 bars (positive proportions, pointing up)
    ax.bar(
        b2["position"],
        b2["proportion"],
        color=cmap_b2(norm(b2["feature_importance_vals"])),
        width=0.5,
        label=biome2
    )

    # labels & formatting
    ax.set_xticks(np.arange(df_final["position"].min(), df_final["position"].max() + 1, 50))
    ax.set_xlabel("MSA Position", fontsize=14)
    ax.set_ylabel("Proportion of most important feature\nper biome (%)", fontsize=14)
    ax.set_title(cluster_num, fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.xticks(rotation=90)

    # colorbars to show importance gradients
    sm_b1 = mpl.cm.ScalarMappable(cmap=cmap_b1, norm=norm)
    sm_b2 = mpl.cm.ScalarMappable(cmap=cmap_b2, norm=norm)
    cbar1 = plt.colorbar(sm_b1, ax=ax, fraction=0.03, pad=0.01)
    cbar2 = plt.colorbar(sm_b2, ax=ax, fraction=0.03, pad=0.01)
    cbar1.set_label(f"Feature importance\n({biome1})")
    cbar2.set_label(f"Feature importance\n({biome2})")

    plt.tight_layout()
    return plt


def plot(path_to_rf_results, cluster_num, biome1, biome2, colour1, colour2, k):
    df_final = make_matrices(path_to_rf_results, cluster_num, biome1, biome2, k)
    plt = plot_bars(df_final, cluster_num, biome1, biome2, colour1, colour2, k)
    return plt