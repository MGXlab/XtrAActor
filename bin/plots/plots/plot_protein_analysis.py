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

plt.rcParams.update({'font.size': 14})   # 🔥 global fontsize

# --- SCALES --- #
# charge scale at pH 7
charge_scale = {
    'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0, 'E': -1, 'Q': 0, 'G': 0,
    'H': 0.5, 'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
}

# Compute the mean flexibility per row
def plot_subplots(df, colour1, colour2):
    df["median_flexibility"] = df["flexibility"].apply(
        lambda x: np.nanmedian(x) if isinstance(x, (list, np.ndarray)) else np.nan
    )

    # --- 1. Expand nested columns ---
    def expand_dict_column(df, col):
        expanded = pd.json_normalize(df[col])
        expanded.columns = [f"{col}_{c}" for c in expanded.columns]
        expanded.index = df.index
        return expanded

    def expand_list_column(df, col):
        expanded = pd.DataFrame(df[col].tolist(), index=df.index)
        expanded.columns = [f"{col}_{i}" for i in range(len(expanded.columns))]
        expanded.index = df.index
        return expanded

    dict_cols = ["count_amino_acids", "amino_acids_percent"]
    list_cols = ["flexibility", "secondary_structure_fraction"]

    expanded_parts = []

    for col in dict_cols:
        expanded = expand_dict_column(df, col)
        expanded_parts.append(expanded.reset_index(drop=True))

    for col in list_cols:
        expanded = expand_list_column(df, col)
        expanded_parts.append(expanded.reset_index(drop=True))

    # Reset df index to make sure it's unique
    df_reset = df.reset_index(drop=True)

    # Safe concatenation
    df_expanded = pd.concat(expanded_parts + [df_reset[["ecosystem_subtype"]]], axis=1)

    def add_ttest_annotations(ax, df, feature, group_col="ecosystem_subtype"):
        groups = df[group_col].unique()
        pairs = list(combinations(groups, 2))
        
        y_max = df[feature].max()
        y_min = df[feature].min()
        y_range = y_max - y_min
        offset = y_range * 0.05  # vertical spacing for annotations
        current_y = y_max + offset

        for g1, g2 in pairs:
            data1 = df.loc[df[group_col] == g1, feature].dropna()
            data2 = df.loc[df[group_col] == g2, feature].dropna()
            if len(data1) < 2 or len(data2) < 2:
                continue
            stat, pval = ttest_ind(data1, data2, equal_var=False)  # Welch's t-test
            # Add annotation
            if pval < 0.001:
                text = "***"
            elif pval < 0.01:
                text = "**"
            elif pval < 0.05:
                text = "*"
            else:
                text = "ns"
            x1, x2 = groups.tolist().index(g1), groups.tolist().index(g2)
            ax.plot([x1, x1, x2, x2], [current_y, current_y + offset/2, current_y + offset/2, current_y], color='k')
            ax.text((x1 + x2) / 2, current_y + offset/2, text, ha='center', va='bottom', fontsize=8)
            current_y += offset

    def plot_pca(data, label, title, ax, colour1,colour2,show_arrows=False, n_arrows=21):
        X = data.select_dtypes(include=[np.number])

        # Drop rows with NaNs and align labels
        mask = ~X.isna().any(axis=1)
        X = X[mask]
        label = label[mask]

        if X.shape[0] < 2 or X.shape[1] < 2:
            ax.set_title(f"{title} (not enough data for PCA)")
            ax.axis("off")
            return

        # Standardize
        X_scaled = StandardScaler().fit_transform(X)

        # PCA
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X_scaled)
        pcs_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'])
        pcs_df['ecosystem_subtype'] = label.values

        # Plot samples
        sns.scatterplot(
            data=pcs_df,
            x='PC1', y='PC2',
            hue='ecosystem_subtype',
            palette=[colour1,colour2],
            ax=ax,
            s=40, alpha=0.8
        )

        ax.set_title(
            f"{title} \n (PCA: {pca.explained_variance_ratio_[0]*100:.1f}% + {pca.explained_variance_ratio_[1]*100:.1f}%)"
        )

        # Optional: add arrows for loadings
        if show_arrows:
            loadings = pca.components_.T  # shape: features x PCs
            feature_names = X.columns

            # Scale arrows for visibility
            arrow_scale = 3 * np.mean(np.abs(pcs)) / np.mean(np.abs(loadings))

            # Pick top n features by length of vector
            magnitudes = np.linalg.norm(loadings[:, :2], axis=1)
            top_idx = np.argsort(magnitudes)[-n_arrows:]

            for i in top_idx:
                x, y = loadings[i, 0] * arrow_scale, loadings[i, 1] * arrow_scale
                ax.arrow(0, 0, x, y, color='black', alpha=0.7, width=0.002, head_width=0.1)
                ax.text(x * 1.1, y * 1.1, feature_names[i][-1], fontsize=12, color='black', ha='center', va='center')

        ax.legend(fontsize=7, loc='best')


    # --- 3. Prepare figure grid ---
    fig, axes = plt.subplots(3, 4, figsize=(15, 10), constrained_layout=True)
    axes = axes.ravel()

    # --- 4. PCA plots (4 of them) ---
    # lake: "#9dcc7e" ocean:"#2B8CBF"
    plot_pca(expand_dict_column(df, "count_amino_acids"), df["ecosystem_subtype"], 
            "Count amino acids", axes[0],colour1, colour2, show_arrows=True)
    plot_pca(expand_dict_column(df, "amino_acids_percent"), df["ecosystem_subtype"], 
            "Amino acids %", axes[1], colour1, colour2, show_arrows=True)
    # plot_pca(expand_list_column(df, "flexibility"), df["ecosystem_subtype"], "Flexibility", axes[2])
    plot_pca(expand_list_column(df, "secondary_structure_fraction"), df["ecosystem_subtype"], "Secondary structure",
            axes[3], colour1, colour2,show_arrows=True)

    # --- 5. Violin plots (6 of them) ---
    violin_features = [
        "molecular_weight", "aromaticity", "instability_index",
        "isoelectric_point", "charge_at_pH_7", "gravy", "median_flexibility"
    ]

    for i, feature in enumerate(violin_features):
        ax = axes[4 + i]
        sns.violinplot(data=df, x="ecosystem_subtype", y=feature, ax=ax, inner="quartile", scale="width", palette=[colour1, colour2])
        ax.set_title(feature.replace("_", " ").capitalize())
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.legend_ = None
        add_ttest_annotations(ax, df, feature)


    # --- 6. Final touches ---
    for j in range(len(violin_features) + 4, 10):
        fig.delaxes(axes[j])  # remove any unused subplots

    plt.tight_layout()
    # plt.show()
    return fig

def plot_scales(df, scale_var, colour1, colour2, title):
    # Create a single figure and axis
    fig, ax = plt.subplots(figsize=(12, 4), constrained_layout=True)

    # Keep only rows with lists/arrays
    flex_df = df[df[scale_var].apply(lambda x: isinstance(x, (list, np.ndarray)))].copy()

    # Convert lists to long-form DataFrame
    records = []
    for idx, row in flex_df.iterrows():
        values = row[scale_var]
        eco = row['ecosystem_subtype']
        for i, val in enumerate(values):
            records.append({'residue_index': i, scale_var: val, 'ecosystem_subtype': eco})

    flex_long = pd.DataFrame(records)

    # Plot with seaborn
    sns.lineplot(
        data=flex_long,
        x='residue_index',
        y=scale_var,
        hue='ecosystem_subtype',
        estimator=np.median,
        ci=95,
        palette=[colour1, colour2],
        ax=ax
    )

    # Labels and title
    ax.set_xlabel("Residue index")
    ax.set_ylabel(scale_var.replace("_", " ").capitalize())
    ax.set_title(title)

    return fig

def main_plot(path_to_rf_results, fasta_path, cluster_num, colour1, colour2,
              flex_window=9,
              polar_window=9,
              charge_window=20):
    
    # paths
    #data_df_path = f"RF_results/raw/{cluster_num}_data_df.tsv"
    #fasta_path = f"RF_results/cluster_fastas/{cluster_num}.faa"
    data_file = path_to_rf_results.joinpath(f"{cluster_num}_data_df.tsv")
    fasta_file = fasta_path.joinpath(f"{cluster_num}.faa")

    # read in files
    all_prot_features = {}
    with open (fasta_file, 'r') as fasta_in:
        records = SeqIO.parse(fasta_in, "fasta")
        for record in records:
            # print(record)
            seq = str(record.seq).replace("X","")
            X = ProteinAnalysis(seq)
            all_prot_features[record.id] = [X.count_amino_acids(),X.amino_acids_percent,X.molecular_weight(),
                                            X.aromaticity(),X.instability_index(),X.isoelectric_point(),X.charge_at_pH(7), X.gravy(),
                                            X.secondary_structure_fraction(),X.molar_extinction_coefficient(),X.protein_scale(ProtParamData.Flex, window=flex_window),X.protein_scale(ProtParamData.kd, window=polar_window), X.protein_scale(charge_scale, window=charge_window)]

    df_analysis = pd.DataFrame.from_dict(all_prot_features, orient='index',
                        columns=['count_amino_acids','amino_acids_percent','molecular_weight',
                                    'aromaticity','instability_index','isoelectric_point','charge_at_pH_7','gravy',
                                    'secondary_structure_fraction','molar_extinction_coefficient','flexibility', 'polarity_profile', 'charge_profile'])
    df_analysis

    # get biome annotations
    data_df = pd.read_csv(data_file,usecols=['protein','ecosystem_subtype'],sep='\t').set_index('protein')
    data_df

    # merge dataframes
    df = df_analysis.join(data_df,how='left')
    df

    # Make da Plots
    plt1 = plot_subplots(df,colour1, colour2)
    plt2 = plot_scales(df, "flexibility", colour1, colour2, title=f"Median flexibility defined by Vihinen 1994 with 95% CI,\n sliding window = {flex_window}")
    plt3 = plot_scales(df, "polarity_profile", colour1, colour2, title=f"Median GRAVY index defined by Kyte & Doolittle hydropathy with 95% CI,\n sliding window = {polar_window}")
    plt4 = plot_scales(df, "charge_profile", colour1, colour2, title=f"Median charge profile at pH 7 with 95% CI,\n sliding window = {charge_window}")
    return plt1, plt2, plt3, plt4