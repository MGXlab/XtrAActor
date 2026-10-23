import pandas as pd
import logomaker
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import seaborn as sns

def make_matrices(path_to_rf_results,cluster_num, biome1, biome2, x=50):

    # read in the data_df for the cluster
    data_file = path_to_rf_results.joinpath(f"{cluster_num}_data_df.tsv")
    data_df_merged_balanced = pd.read_csv(data_file, sep='\t')
    data_df_merged_balanced['biome'] = data_df_merged_balanced['biome'].str.replace(":",";")
    print(data_df_merged_balanced)
    # print(data_df_merged_balanced.groupby('ecosystem_subtype')['179_G'].sum())

    # Read in feature importance scores
    feat_file = path_to_rf_results.joinpath(f"{cluster_num}_Feature_importance.tsv")
    feat_imp_df = pd.read_csv(feat_file, sep='\t', usecols=[0,2],names=['col_name','feature_importance_vals'],header=0).set_index('col_name')
    print(feat_imp_df)

    # top_50_ft_df = feat_imp_df
    # # top_50_ft_df = feat_imp_df.head(n=x)

    # # top_50_ft_df[['position','aa']] = top_50_ft_df['col_name'].str.split("_",expand=True)
    # # top_50_ft_df['position'] = top_50_ft_df['position'].astype(int)
 

    data_df_merged_balanced_1 = data_df_merged_balanced.groupby('ecosystem_subtype').apply(lambda x: (x==1).sum()).T
    data_df_merged_balanced_1[f"Proportion_{biome1}"] = (
        data_df_merged_balanced_1[biome1] / data_df_merged_balanced[data_df_merged_balanced['ecosystem_subtype'] == biome1].shape[0] * 100)*-1
    data_df_merged_balanced_1[f"Proportion_{biome2}"] = (
        data_df_merged_balanced_1[biome2] / data_df_merged_balanced[data_df_merged_balanced['ecosystem_subtype'] == biome2].shape[0] * 100)
    data_df_merged_balanced_1.drop(['protein','biome','ecosystem','ecosystem_category','ecosystem_type','ecosystem_subtype','specific_ecosystem'], axis=0, inplace=True)
    print(data_df_merged_balanced_1)

    temp = data_df_merged_balanced_1.join(feat_imp_df,how='left').reset_index(names='pos_AA')
    temp_m = pd.melt(temp,id_vars=['pos_AA','feature_importance_vals'],var_name='biome',value_name='proportion')
    print(temp_m)
    temp_m = temp_m.loc[temp_m['biome'].str.contains("Proportion")]
    temp_m[['position','AA']] = temp_m['pos_AA'].str.split('_',expand=True)
    df_final = temp_m.sort_values("feature_importance_vals", ascending=False).drop_duplicates(['position','biome'])
    imp_min = df_final['feature_importance_vals'].min()
    imp_max = df_final['feature_importance_vals'].max()
    df_final['feat_imp_normal'] = (df_final['feature_importance_vals'] - imp_min) / (imp_max - imp_min + 1E-12)
    df_final['feat_imp_normal'] = df_final['feat_imp_normal'].astype(float)
    df_final = df_final.dropna(subset=["position", "proportion"])
    df_final["proportion"] = df_final["proportion"].astype(float)
    df_final["position"] = df_final["position"].astype(int)
    # print(df_final["position"].head(10))
    # print(df_final["position"].dtype)
    # print(df_final["position"].apply(type).value_counts())
    # print(df_final[["position", "proportion", "feat_imp_normal", "biome"]].isna().sum())
    # print(df_final[["position", "proportion", "feat_imp_normal"]].describe())

    print(df_final.dtypes)

    return df_final

# def plot_bars(df_final, cluster_num,colour1, colour2, x=50):
#     sns.set_theme(style="whitegrid")
#     my_palette = sns.diverging_palette(h_neg=100, h_pos=200,s=100,l=60, n=2)
#     min,max = 2, 50
#     plt.figure(figsize=(30,5))

#     # b = sns.barplot(
#     #     data=top_features_T_m,
#     #     x="position",
#     #     y="value",
#     #     hue="ecosystem_subtype",
#     #     palette=[colour1,colour2],
#     #     dodge=False,
#     #     width=3  
#     # )

#     b = sns.scatterplot(
#         data=df_final,
#         x="position",
#         y="proportion",
#         hue="feat_imp_normal",
#         palette="flare",
#         size="feat_imp_normal",
#         sizes = (20,200),
#         alpha = 0.75)

# #     alphas = (df_final["feature_importance_vals"] - df_final["feature_importance_vals"].min()) / \
# #             (df_final["feature_importance_vals"].max() - df_final["feature_importance_vals"].min())
# #     alphas = 0.2 + alphas * 0.8

# #     plt.scatter(
# #         x=df_final["position"],
# #         y=df_final["proportion"],
# #         c=df_final["biome"].map({"Biome1": colour1, "Biome2": colour2}),  # adjust to your mapping
# #         s=df_final["feat_imp_normal"] * 200,  # scale size if desired
# # )

#     # show the graph
#     # b.set_title(f"{cluster_num} Top{x} features",fontsize=14)
#     # b.set_xlabel("Proportion of feature in sequences from either biome (%)",fontsize=14)
#     # b.set_ylabel(f"Top {x} Important features",fontsize=14)

#     # plt.show()
#     return plt

import matplotlib as mpl

def plot_scatter(df_final, cluster_num,colour1, colour2,k):
    # Define the two colormaps
    cmap_lake = mpl.colors.LinearSegmentedColormap.from_list("lake_cmap", ["#fafaf0", colour1])
    cmap_ocean = mpl.colors.LinearSegmentedColormap.from_list("ocean_cmap", ["#fafaf0", colour2])

    def biome_color(row):
        if row["biome"] == "Proportion_Lake":
            return cmap_lake(row["feat_imp_normal"])
        elif row["biome"] == "Proportion_Oceanic":
            return cmap_ocean(row["feat_imp_normal"])
        else:
            return (0.8, 0.8, 0.8, 1.0)  # fallback gray

    df_final["color"] = df_final.apply(biome_color, axis=1)

    fig, ax = plt.subplots(figsize=(30,5))

    sc1 = ax.scatter(
        df_final.loc[df_final["biome"]=="Proportion_Lake", "position"],
        df_final.loc[df_final["biome"]=="Proportion_Lake", "proportion"],
        s=df_final.loc[df_final["biome"]=="Proportion_Lake", "feat_imp_normal"] * 200,
        c=df_final.loc[df_final["biome"]=="Proportion_Lake", "feat_imp_normal"],
        cmap=cmap_lake,
        alpha=0.9
    )

    sc2 = ax.scatter(
        df_final.loc[df_final["biome"]=="Proportion_Oceanic", "position"],
        df_final.loc[df_final["biome"]=="Proportion_Oceanic", "proportion"],
        s=df_final.loc[df_final["biome"]=="Proportion_Oceanic", "feat_imp_normal"] * 200,
        c=df_final.loc[df_final["biome"]=="Proportion_Oceanic", "feat_imp_normal"],
        cmap=cmap_ocean,
        alpha=0.9
    )
    plt.set_xticks(np.arange(df_final["position"].min(),
                        df_final["position"].max()+1, 50), rotation=90)
    plt.colorbar(sc1, ax=ax, label="Feature importance (Lake)")
    plt.colorbar(sc2, ax=ax, label="Feature importance (Ocean)")
    
    plt.tight_layout()
    return plt

def plot_bars(df_final, cluster_num,biome1, biome2,colour1, colour2,k):
    # Define colormaps
    cmap_lake = mpl.colors.LinearSegmentedColormap.from_list("lake_cmap", ["#fafaf0", "#9dcc7e"])
    cmap_ocean = mpl.colors.LinearSegmentedColormap.from_list("ocean_cmap", ["#fafaf0", "#2B8CBF"])

    # Normalize for mapping feature importance to color range
    norm = mpl.colors.Normalize(vmin=df_final["feature_importance_vals"].min(), vmax=df_final["feature_importance_vals"].max())

    # Separate data by biome
    lake = df_final[df_final["biome"] == f"Proportion_{biome1}"].sort_values("position")
    ocean = df_final[df_final["biome"] == f"Proportion_{biome2}"].sort_values("position")

    fig, ax = plt.subplots(figsize=(30, 6))

    # Plot bars for Lake
    ax.bar(
        lake["position"],
        lake["proportion"],
        color=cmap_lake(norm(lake["feature_importance_vals"])),
        width=0.5,
        label=biome1
    )

    # Plot bars for Ocean (offset a little for side-by-side effect)
    ax.bar(
        ocean["position"].astype(float),   # slight x-shift to avoid overlap
        ocean["proportion"],
        color=cmap_ocean(norm(ocean["feature_importance_vals"])),
        width=0.5,
        label="Ocean"
    )

    # Labels & formatting
    ax.set_xticks(np.arange(df_final["position"].min(),
                    df_final["position"].max()+1, 50))
    ax.set_xlabel("MSA Position", fontsize=14)
    ax.set_ylabel(f"Proportion of Most Imp Feature in {biome1} or {biome2} Sequences(%)", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.xticks(rotation=90)

    # Add colorbars to show importance gradients
    sm_lake = mpl.cm.ScalarMappable(cmap=cmap_lake, norm=norm)
    sm_ocean = mpl.cm.ScalarMappable(cmap=cmap_ocean, norm=norm)
    cbar1 = plt.colorbar(sm_lake, ax=ax, fraction=0.03, pad=0.01)
    cbar2 = plt.colorbar(sm_ocean, ax=ax, fraction=0.03, pad=0.01)
    cbar1.set_label(f"Feature importance ({biome1})")
    cbar2.set_label(f"Feature importance ({biome2})")

    plt.tight_layout()
    # plt.show()
    return plt


def plot(path_to_rf_results,cluster_num, biome1, biome2, colour1, colour2,k):
    df_final = make_matrices(path_to_rf_results,cluster_num, biome1, biome2, k)
    plt = plot_bars(df_final, cluster_num,biome1, biome2,colour1, colour2,k)
    return plt