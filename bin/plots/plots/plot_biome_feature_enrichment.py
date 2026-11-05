import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
                              names=['col_name', 'feature_importance_vals'], header=0)
    top_50_ft_df = feat_imp_df.head(n=x).copy()
    top_50_ft_df[['position', 'aa']] = top_50_ft_df['col_name'].str.split("_", n=1, expand=True)
    top_50_ft_df['position'] = top_50_ft_df['position'].astype(int)

    # number of sequences carrying each feature, per biome (rows = features, columns = biomes)
    counts = data_df_merged_balanced.groupby('biome')[feature_cols].apply(lambda d: (d == 1).sum()).T

    # proportion (%) of sequences per biome carrying each feature
    n_biome1 = (data_df_merged_balanced['biome'] == biome1).sum()
    n_biome2 = (data_df_merged_balanced['biome'] == biome2).sum()
    proportions = pd.DataFrame({
        f"Proportion_{biome1}": counts[biome1] / n_biome1 * 100 * -1,   # negative -> bars to the left
        f"Proportion_{biome2}": counts[biome2] / n_biome2 * 100,
    })

    # keep only the top features, in importance order
    top_cols = [c for c in top_50_ft_df['col_name'] if c in proportions.index]
    top_features = proportions.loc[top_cols].reset_index(names='pos_AA')

    top_features_T_m = top_features.melt(id_vars="pos_AA", var_name="biome", value_name="value")
    return top_features_T_m


def plot_bars(top_features_T_m, cluster_num, biome1, biome2, colour1, colour2, x=50):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, x / 2))

    b = sns.barplot(
        data=top_features_T_m,
        x="value",
        y="pos_AA",
        hue="biome",
        hue_order=[f"Proportion_{biome1}", f"Proportion_{biome2}"],
        palette=[colour1, colour2],
        dodge=False,
        width=0.5,
    )
    b.tick_params(axis='x', labelsize=16)
    b.tick_params(axis='y', labelsize=16)
    b.set_title(f"{cluster_num} Top{x} features", fontsize=16)
    b.set_xlabel("Proportion of feature in sequences from either biome (%)", fontsize=16)
    b.set_ylabel(f"Top {x} Important features", fontsize=16)

    plt.tight_layout()
    return plt


def plot(path_to_rf_results, cluster_num, biome1, biome2, colour1, colour2, k):
    top_features_T_m = make_matrices(path_to_rf_results, cluster_num, biome1, biome2, k)
    plt = plot_bars(top_features_T_m, cluster_num, biome1, biome2, colour1, colour2, k)
    return plt