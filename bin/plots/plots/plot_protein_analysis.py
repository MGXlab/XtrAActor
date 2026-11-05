from Bio import SeqIO
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from itertools import combinations
from scipy.stats import ttest_ind

plt.rcParams.update({'font.size': 14})   # global fontsize

# Expected inputs:
#   <path_to_rf_results>/<cluster>_data_df.tsv : gene | <features> | biome (full biome string)
#   <fasta_path>/<cluster>.fasta (or .faa)     : unaligned protein sequences of the cluster

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# --- SCALES --- #
# charge scale at pH 7
charge_scale = {
    'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0, 'E': -1, 'Q': 0, 'G': 0,
    'H': 0.5, 'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
}


def plot_subplots(df, order, palette):
    """order = [biome1, biome2]; palette = {biome: colour}"""
    colours = [palette[b] for b in order]

    df["median_flexibility"] = df["flexibility"].apply(
        lambda x: np.nanmedian(x) if isinstance(x, (list, np.ndarray)) and len(x) > 0 else np.nan
    )

    def expand_dict_column(df, col):
        expanded = pd.json_normalize(df[col])
        expanded.columns = [f"{col}_{c}" for c in expanded.columns]
        expanded.index = df.index
        return expanded

    def expand_list_column(df, col):
        expanded = pd.DataFrame(df[col].tolist(), index=df.index)
        expanded.columns = [f"{col}_{i}" for i in range(len(expanded.columns))]
        return expanded

    def add_ttest_annotations(ax, df, feature, group_col="biome"):
        y_max = df[feature].max()
        y_min = df[feature].min()
        y_range = y_max - y_min
        offset = (y_range if y_range > 0 else 1) * 0.05   # vertical spacing for annotations
        current_y = y_max + offset

        for g1, g2 in combinations(order, 2):
            data1 = df.loc[df[group_col] == g1, feature].dropna()
            data2 = df.loc[df[group_col] == g2, feature].dropna()
            if len(data1) < 2 or len(data2) < 2:
                continue
            stat, pval = ttest_ind(data1, data2, equal_var=False)  # Welch's t-test
            if pval < 0.001:
                text = "***"
            elif pval < 0.01:
                text = "**"
            elif pval < 0.05:
                text = "*"
            else:
                text = "ns"
            x1, x2 = order.index(g1), order.index(g2)   # same order as the plotted violins
            ax.plot([x1, x1, x2, x2], [current_y, current_y + offset/2, current_y + offset/2, current_y], color='k')
            ax.text((x1 + x2) / 2, current_y + offset/2, text, ha='center', va='bottom', fontsize=8)
            current_y += offset

    def plot_pca(data, label, title, ax, show_arrows=False, n_arrows=21):
        X = data.select_dtypes(include=[np.number])

        # Drop rows with NaNs and align labels
        mask = ~X.isna().any(axis=1)
        X = X[mask]
        label = label[mask]

        if X.shape[0] < 2 or X.shape[1] < 2:
            ax.set_title(f"{title} (not enough data for PCA)")
            ax.axis("off")
            return

        X_scaled = StandardScaler().fit_transform(X)

        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X_scaled)
        pcs_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'])
        pcs_df['biome'] = label.values

        sns.scatterplot(
            data=pcs_df,
            x='PC1', y='PC2',
            hue='biome', hue_order=order, palette=palette,
            ax=ax,
            s=40, alpha=0.8
        )

        ax.set_title(
            f"{title} \n (PCA: {pca.explained_variance_ratio_[0]*100:.1f}% + {pca.explained_variance_ratio_[1]*100:.1f}%)"
        )

        if show_arrows:
            loadings = pca.components_.T  # features x PCs
            feature_names = X.columns

            arrow_scale = 3 * np.mean(np.abs(pcs)) / np.mean(np.abs(loadings))

            magnitudes = np.linalg.norm(loadings[:, :2], axis=1)
            top_idx = np.argsort(magnitudes)[-n_arrows:]

            for i in top_idx:
                x, y = loadings[i, 0] * arrow_scale, loadings[i, 1] * arrow_scale
                ax.arrow(0, 0, x, y, color='black', alpha=0.7, width=0.002, head_width=0.1)
                ax.text(x * 1.1, y * 1.1, feature_names[i][-1], fontsize=12, color='black', ha='center', va='center')

        ax.legend(fontsize=7, loc='best')

    # --- figure grid: 3 PCA plots + 7 violin plots ---
    fig, axes = plt.subplots(3, 4, figsize=(15, 10), constrained_layout=True)
    axes = axes.ravel()

    # --- PCA plots ---
    plot_pca(expand_dict_column(df, "count_amino_acids"), df["biome"],
             "Count amino acids", axes[0], show_arrows=True)
    plot_pca(expand_dict_column(df, "amino_acids_percent"), df["biome"],
             "Amino acids %", axes[1], show_arrows=True)
    plot_pca(expand_list_column(df, "secondary_structure_fraction"), df["biome"],
             "Secondary structure", axes[2], show_arrows=True)

    # --- violin plots ---
    violin_features = [
        "molecular_weight", "aromaticity", "instability_index",
        "isoelectric_point", "charge_at_pH_7", "gravy", "median_flexibility"
    ]

    for i, feature in enumerate(violin_features):
        ax = axes[3 + i]
        sns.violinplot(data=df, x="biome", y=feature, order=order, palette=colours,
                       ax=ax, inner="quartile", scale="width")
        ax.set_title(feature.replace("_", " ").capitalize())
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([b.split(";")[-1] for b in order])   # short labels on the x axis
        ax.legend_ = None
        add_ttest_annotations(ax, df, feature)

    # --- remove unused subplots ---
    for j in range(3 + len(violin_features), len(axes)):
        fig.delaxes(axes[j])

    return fig


def plot_scales(df, scale_var, order, palette, title):
    fig, ax = plt.subplots(figsize=(12, 4), constrained_layout=True)

    # Keep only rows with lists/arrays
    flex_df = df[df[scale_var].apply(lambda x: isinstance(x, (list, np.ndarray)))].copy()

    # Convert lists to long-form DataFrame
    records = []
    for idx, row in flex_df.iterrows():
        for i, val in enumerate(row[scale_var]):
            records.append({'residue_index': i, scale_var: val, 'biome': row['biome']})
    flex_long = pd.DataFrame(records)

    if flex_long.empty:
        ax.set_title(f"{title}\n(no data)")
        ax.axis("off")
        return fig

    sns.lineplot(
        data=flex_long,
        x='residue_index',
        y=scale_var,
        hue='biome', hue_order=order, palette=palette,
        estimator=np.median,
        ci=95,
        ax=ax
    )

    ax.set_xlabel("Residue index")
    ax.set_ylabel(scale_var.replace("_", " ").capitalize())
    ax.set_title(title)

    return fig


def main_plot(path_to_rf_results, fasta_path, cluster_num, colour1, colour2,
              biome1=None, biome2=None,
              flex_window=9,
              polar_window=9,
              charge_window=20):

    data_file = path_to_rf_results.joinpath(f"{cluster_num}_data_df.tsv")

    # cluster fasta: .fasta (genome entry point) or .faa (clusters entry point)
    fasta_file = None
    for ext in (".fasta", ".faa", ".fa"):
        candidate = fasta_path.joinpath(f"{cluster_num}{ext}")
        if candidate.exists():
            fasta_file = candidate
            break
    if fasta_file is None:
        raise FileNotFoundError(f"No fasta for {cluster_num} in {fasta_path} (.fasta/.faa/.fa)")

    # protein properties per sequence
    all_prot_features = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # keep only the 20 standard amino acids (drops X, *, B, Z, U, ...)
        seq = "".join(c for c in str(record.seq).upper() if c in VALID_AA)
        if len(seq) == 0:
            continue
        X = ProteinAnalysis(seq)
        all_prot_features[record.id] = [
            X.count_amino_acids(), X.amino_acids_percent, X.molecular_weight(),
            X.aromaticity(), X.instability_index(), X.isoelectric_point(), X.charge_at_pH(7), X.gravy(),
            X.secondary_structure_fraction(), X.molar_extinction_coefficient(),
            X.protein_scale(ProtParamData.Flex, window=flex_window),
            X.protein_scale(ProtParamData.kd, window=polar_window),
            X.protein_scale(charge_scale, window=charge_window)
        ]

    df_analysis = pd.DataFrame.from_dict(
        all_prot_features, orient='index',
        columns=['count_amino_acids', 'amino_acids_percent', 'molecular_weight',
                 'aromaticity', 'instability_index', 'isoelectric_point', 'charge_at_pH_7', 'gravy',
                 'secondary_structure_fraction', 'molar_extinction_coefficient',
                 'flexibility', 'polarity_profile', 'charge_profile'])

    # biome annotations from the RF data_df (gene IDs = fasta record IDs)
    data_df = pd.read_csv(data_file, usecols=['gene', 'biome'], sep='\t').set_index('gene')
    data_df['biome'] = data_df['biome'].str.replace(":", ";")

    # keep only genes with a biome label (inner join)
    df = df_analysis.join(data_df, how='inner')

    # biome order and colours: from the arguments if given, otherwise sorted
    if biome1 and biome2:
        order = [biome1, biome2]
        df = df[df['biome'].isin(order)]
    else:
        order = sorted(df['biome'].unique())[:2]
    palette = dict(zip(order, [colour1, colour2]))

    if df.empty:
        raise ValueError(f"No sequences of {cluster_num} matched the data_df genes / biomes {order}")

    plt1 = plot_subplots(df, order, palette)
    plt2 = plot_scales(df, "flexibility", order, palette,
                       title=f"Median flexibility defined by Vihinen 1994 with 95% CI,\n sliding window = {flex_window}")
    plt3 = plot_scales(df, "polarity_profile", order, palette,
                       title=f"Median GRAVY index defined by Kyte & Doolittle hydropathy with 95% CI,\n sliding window = {polar_window}")
    plt4 = plot_scales(df, "charge_profile", order, palette,
                       title=f"Median charge profile at pH 7 with 95% CI,\n sliding window = {charge_window}")
    return plt1, plt2, plt3, plt4