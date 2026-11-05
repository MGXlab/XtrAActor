import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

plt.rcParams.update({'font.size': 14})   # 🔥 global fontsize


# --- Define amino acid order by chemical properties ---
aa_groups = [
    ("Positive charged", ["R", "H", "K"], "crimson"),
    ("Negative charged", ["D", "E"], "royalblue"),
    ("Polar uncharged", ["S", "T", "N", "Q"], "seagreen"),
    ("Hydrophobic", ["A", "V", "I", "L", "M"], "goldenrod"),
    ("Aromatic", ["F", "Y", "W"], "purple"),
    ("Special cases", ["C", "U", "G", "P"], "darkorange")
]

aa_order = [aa for _, group, _ in aa_groups for aa in group]
base_colors = {aa: color for _, group, color in aa_groups for aa in group}

def lighten_color(color, factor=0.5):
    """Lighten color by blending with white"""
    c = np.array(mcolors.to_rgb(color))
    return tuple(c + (1 - c) * factor)

def darken_color(color, factor=0.5):
    """Darken color by blending with black"""
    c = np.array(mcolors.to_rgb(color))
    return tuple(c * (1 - factor))

def get_aa(idx):
    """Extract AA letter from feature index like '152_G' -> 'G'"""
    return idx.split("_")[1]

def radar_biome_plot(cluster_num, biome1, biome2, colour1, colour2,path_to_rf_results,top_k=50):
    # Get the dataframes from the files
    feat_file = path_to_rf_results.joinpath(f"{cluster_num}_Feature_importance.tsv")
    data_file = path_to_rf_results.joinpath(f"{cluster_num}_data_df.tsv")
    cluster_df = pd.read_csv(feat_file,sep='\t',header=0,index_col=False)
    biome_df = pd.read_csv(data_file, sep='\t',header=0,index_col=False)
    print(cluster_df)

    # Top-k features for this cluster
    top_feats = cluster_df.nlargest(top_k, "gini_imp").iloc[:, 0].tolist()
    feat_cols = [f for f in top_feats if f in biome_df.columns]

    if not feat_cols:
        print(f"No matching features for {cluster_num}")
        return

    sub = biome_df[feat_cols + ["biome"]]

    biome_counts = {}
    for biome in [biome1, biome2]:
        mask = sub["biome"] == biome
        counts = sub.loc[mask, feat_cols].sum()
        aa_counts = counts.rename(lambda x: get_aa(x)).groupby(level=0).sum()
        aa_counts = aa_counts.reindex(aa_order, fill_value=0)

        total = aa_counts.sum()
        if total > 0:
            aa_counts = aa_counts / total
        biome_counts[biome] = aa_counts

    # ---- Radar plot ----
    N = len(aa_order)
    dtheta = 2*np.pi / N
    angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    plt.figure(figsize=(9,9))
    ax = plt.subplot(111, polar=True)

    for biome, color in zip([biome1, biome2], [colour1, colour2]):
        if biome not in biome_counts:
            continue
        values = biome_counts[biome].tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=2, label=biome, color=color, alpha=1)
        ax.fill(angles, values, alpha=0.5, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(aa_order)

    # ---- Group coloured bands OUTSIDE the radar ----
    outer_radius = ax.get_rmax() * 1.05   # just outside the circle
    bar_height = ax.get_rmax() * 0.03     # thickness of band

    legend_patches = []
    for group_name, group_aas, group_color in aa_groups:
        start_idx = aa_order.index(group_aas[0])
        end_idx = aa_order.index(group_aas[-1])

        # Start halfway before first AA, end halfway after last AA
        start_angle = start_idx * dtheta - dtheta/2
        end_angle = (end_idx+1) * dtheta - dtheta/2
        width = end_angle - start_angle

        ax.bar(
            x=start_angle, height=bar_height, width=width,
            bottom=outer_radius * 1.05, color=group_color, alpha=0.6, align="edge"
        )
        legend_patches.append(Patch(color=group_color, alpha=0.6, label=group_name))

    # ---- Legend for groups ----
    # Legend for biomes
    biome_handles = [
        Line2D([0], [0], color="#2B8CBF", lw=3, label=biome1),
        Line2D([0], [0], color="#9dcc7e", lw=3, label=biome2),
    ]
    # Combine group and biome legends
    handles = legend_patches + biome_handles
    ax.legend(handles=handles, title="AA Groups & Biomes",
              loc="upper right", bbox_to_anchor=(1.2, 1.1))
    ax.set_title(f"Top {top_k} features - {cluster_num}", size=14, pad=20)
    return plt


def AA_bar_plot(cluster_num, biomes,path_to_rf_results,top_k=50):
    
    # Get the dataframes from the files
    feat_file = path_to_rf_results.joinpath(f"{cluster_num}_Feature_importance.tsv")
    data_file = path_to_rf_results.joinpath(f"{cluster_num}_data_df.tsv")
    cluster_df = pd.read_csv(feat_file,sep='\t')
    biome_df = pd.read_csv(data_file, sep='\t')

    # --- Select top features ---
    top_feats = cluster_df.nlargest(top_k, "gini_imp").iloc[:, 0].tolist()
    feat_cols = [f for f in top_feats if f in biome_df.columns]

    if not feat_cols:
        print(f"No matching features for {cluster_num}")
        return

    sub = biome_df[feat_cols + ["biome"]]

    # --- Count amino acid frequencies per biome ---
    biome_freqs = {}
    for biome in biomes:
        mask = sub["biome"] == biome
        counts = sub.loc[mask, feat_cols].sum()
        aa_counts = counts.rename(lambda x: get_aa(x)).groupby(level=0).sum()
        aa_counts = aa_counts.reindex(aa_order, fill_value=0)

        total = aa_counts.sum()
        if total > 0:
            aa_counts = aa_counts / total  # normalize to relative frequency
        biome_freqs[biome] = aa_counts

    df_plot = pd.DataFrame(biome_freqs).fillna(0)

    # --- Plot ---
    fig, axes = plt.subplots(2, 1, figsize=(14,12), sharex=False)

       # 1. Grouped bar chart
    bar_width = 0.4
    x = np.arange(len(aa_order))

    for biome, offset in zip(biomes, [-bar_width/2, bar_width/2]):
        values = df_plot[biome]
        if biome == biomes[0]:
            colors = [lighten_color(base_colors[aa], 0.6) for aa in aa_order]
        else:  # Oceanic
            colors = [darken_color(base_colors[aa], 0.3) for aa in aa_order]

        axes[0].bar(x + offset, values, width=bar_width, label=biome, alpha=0.9, color=colors)

    axes[0].set_xticks(x)
    axes[0].set_xticklabels(aa_order)
    axes[0].set_ylabel("Relative frequency")
    axes[0].set_title(f"Top {top_k} features - {cluster_num} (AA frequencies)")

    # --- Build legend manually ---
    biome_patches = [
        mpatches.Patch(color="lightgrey", label=f"{biomes[0]} (lighter)"),
        mpatches.Patch(color="dimgray", label=f"{biomes[1]} (darker)")
    ]
    group_patches = [mpatches.Patch(color=color, label=name) for name, _, color in aa_groups]

    axes[0].legend(handles=group_patches + biome_patches,
                   loc="upper right", bbox_to_anchor=(1.15, 1.05))


    # 2. Difference plot (Lake - Oceanic)
    diff = df_plot[biomes[1]] - df_plot[biomes[0]]
    diff_colors = [
        darken_color(base_colors[aa], 0.3) if v > 0 else lighten_color(base_colors[aa], 0.6)
        for aa, v in zip(aa_order, diff)
    ]

    axes[1].bar(x, diff, color=diff_colors, alpha=0.9)
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(aa_order)
    axes[1].set_ylabel(f"{biomes[0]} - {biomes[1]}")
    axes[1].set_title("Difference in normalized frequencies (enrichment)")

    plt.tight_layout()
    # plt.show()
    return plt
