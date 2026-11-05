import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text

plt.rcParams.update({'font.size': 14})   # 🔥 global fontsize


def plot_report(df, biome1, biome2, colour1, colour2):

    # Melt for seaborn long-form
    df_long = df.melt(id_vars=["clust_num", "biome", "support"], 
                    value_vars=["precision", "recall", "f1", "accuracy"],
                    var_name="metric", value_name="value")

    # --- Plotting ---
    palette = {biome1: colour1, biome2: colour2}

    fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=True)

    ax_dict = {
        "accuracy": axes[0, 0],
        "precision": axes[0, 1],
        "recall": axes[1, 1],
        "f1": axes[1, 0]
    }

    legend_handles = None

    for metric, ax in ax_dict.items():
        if metric == "accuracy":
            sns.barplot(
                data=df_long[df_long.metric == metric].drop_duplicates(subset=["clust_num"]),
                x="clust_num", y="value", color="dimgray", ax=ax, width=0.6
            )
        else:
            plot_data = df_long[df_long.metric == metric].copy()
            sns.barplot(
                data=plot_data,
                x="clust_num", y="value", hue="biome", palette=palette, ax=ax, width=0.6
            )
            if legend_handles is None:
                legend_handles = ax.get_legend_handles_labels()
            ax.get_legend().remove()

            # Annotate support
            for i, row in plot_data.iterrows():
                if row.metric == "f1":
                    # x-coordinate depends on cluster position and hue
                    x = list(plot_data.clust_num.unique()).index(row.clust_num)
                    offset = -0.2 if row.biome == biome1 else 0.2
                    ax.text(x + offset, row.value + 0.02, str(int(row.support)),
                            ha="center", va="bottom")

        ax.set_ylim(0, 1.15)
        ax.set_ylabel(metric.capitalize())
        ax.set_xlabel("")
        # ax.set_title(f"{metric.capitalize()}")
        ax.set_xticklabels(df["clust_num"].unique(), rotation=45, ha="right")
        ax.tick_params(axis='both', which='major', labelsize=14)

    # Single legend
    fig.legend(legend_handles[0], legend_handles[1], title="Biome", loc="upper center", ncol=2)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return plt

def plot_precision_vs_recall(df, biome1, biome2, colour1, colour2):

    fig, ax = plt.subplots(figsize=(7,5))

    texts = []  # To store label objects for adjustText

    # Plot points
    for i, row in df.iterrows():
        x = row["precision"]
        y = row["recall"]
        c = colour1 if row["biome"] == biome1 else colour2
        ax.scatter(x, y, color=c, s=150, edgecolor='k')

        # Add label object (will be adjusted later)
        texts.append(ax.text(x, y, row["clust_num"], ha="center", va="center"))

    # Adjust labels to avoid overlap
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5), expand_text=(14,14))

    ax.set_xlabel("Precision")
    ax.set_ylabel("Recall")
    ax.set_title("Precision vs Recall per Cluster")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)

    # Single legend for biomes
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colour1, markersize=10, label=biome1),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colour2, markersize=10, label=biome2)
    ]
    ax.legend(handles=handles, title="Biome", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    return plt