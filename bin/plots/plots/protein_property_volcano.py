#!/usr/bin/env python3
"""
Compare all multi-biome clusters at once: one volcano plot per protein property.

Properties (per protein, calculated as in the per-cluster protein analysis plots):
  molecular weight, aromaticity, instability index, isoelectric point, charge at pH 7,
  GRAVY index, median flexibility (Vihinen 1994, sliding window 9)

For every cluster and every property:
  x = mean of the biome1 proteins - mean of the biome2 proteins
  p = Wilcoxon rank-sum test (Mann-Whitney U, two-sided) of the per-protein values,
      biome1 proteins vs. biome2 proteins of that cluster
  padj = Benjamini-Hochberg correction over all clusters, separately for each property
  y = -log10(padj)

Every cluster is plotted; only clusters with padj < 0.05 get an annotation label.

Inputs:
  rf_results_dir : folder with <cluster>_data_df.tsv   (gene | <features> | biome)
  fasta_dir      : folder with <cluster>.fasta / .faa / .fa (unaligned proteins; header = ><gene> <annotation>)

Outputs (in --out_dir):
  protein_properties_per_protein.tsv   one row per protein: cluster, gene, biome, the 7 properties
  protein_properties_cluster_stats.tsv one row per cluster and property: n, means, diff, p_value, padj
  volcano_<property>.png               one plot per property

Usage:
  protein_property_volcano.py <rf_results_dir> <fasta_dir> --biome1 "..." --biome2 "..." \\
      --colour1 "#0bd600" --colour2 "#8e00fa" [--out_dir .]
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")            # no display on the cluster nodes
import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.transforms
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import mannwhitneyu

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
UNINFORMATIVE = ("hypothetical", "unknown function", "uncharacterized", "uncharacterised")
FLEX_WINDOW = 9
MIN_PER_BIOME = 3        # proteins per biome needed for the test
PADJ_CUTOFF = 0.05
POINT_ALPHA = 0.7        # transparency of the points

# column name -> axis label
PROPERTIES = {
    "molecular_weight":   "Molecular weight (Da)",
    "aromaticity":        "Aromaticity",
    "instability_index":  "Instability index",
    "isoelectric_point":  "Isoelectric point",
    "charge_at_pH_7":     "Charge at pH 7",
    "gravy":              "GRAVY index",
    "median_flexibility": "Median flexibility",
}


def short(label):
    return label.split(";")[-1].strip()


def match_biome(wanted, present):
    """Full biome label in the data for a --biome argument (full label or part after the last ';')."""
    wanted = wanted.replace(":", ";")
    hits = [b for b in present if b == wanted or short(b) == short(wanted)]
    return hits[0] if len(hits) == 1 else None


def bh_adjust(pvalues):
    """Benjamini-Hochberg adjusted p-values; NaNs stay NaN and are not counted as tests."""
    p = np.asarray(pvalues, dtype=float)
    padj = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = ok.sum()
    if m == 0:
        return padj
    order = np.argsort(p[ok])
    ranked = p[ok][order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    padj[ok] = out
    return padj


def protein_properties(fasta_file):
    """Per-protein properties + the cluster annotation (most frequent informative header text)."""
    rows = {}
    annotations = Counter()
    for record in SeqIO.parse(fasta_file, "fasta"):
        desc = record.description[len(record.id):].strip()
        if desc and not any(u in desc.lower() for u in UNINFORMATIVE):
            annotations[desc] += 1
        # keep only the 20 standard amino acids (drops X, *, B, Z, U, ...)
        seq = "".join(c for c in str(record.seq).upper() if c in VALID_AA)
        if len(seq) == 0:
            continue
        X = ProteinAnalysis(seq)
        flex = X.protein_scale(ProtParamData.Flex, window=FLEX_WINDOW) if len(seq) >= FLEX_WINDOW else []
        rows[record.id] = {
            "molecular_weight":   X.molecular_weight(),
            "aromaticity":        X.aromaticity(),
            "instability_index":  X.instability_index(),
            "isoelectric_point":  X.isoelectric_point(),
            "charge_at_pH_7":     X.charge_at_pH(7),
            "gravy":              X.gravy(),
            "median_flexibility": float(np.nanmedian(flex)) if len(flex) > 0 else np.nan,
        }
    annotation = annotations.most_common(1)[0][0] if annotations else "hypothetical protein"
    return pd.DataFrame.from_dict(rows, orient="index"), annotation


def cluster_label(cluster, annotation, max_len=32):
    number = cluster.replace("cluster_", "")
    text = annotation if len(annotation) <= max_len else annotation[:max_len - 1] + "…"
    return f"{text} ({number})"


def text_box(ann, renderer):
    """Bounding box of an annotation's text only (without the connecting line), in pixels."""
    ann.update_positions(renderer)
    return matplotlib.text.Text.get_window_extent(ann, renderer).padded(2)


def volcano(stats, prop, biome1, biome2, colour1, colour2, out_file):
    d = stats[(stats["property"] == prop) & stats["padj"].notna()].copy()
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    title = f"{PROPERTIES[prop]}: {short(biome1)} vs. {short(biome2)}"
    if d.empty:
        ax.set_title(f"{title}\n(no cluster with enough proteins per biome)")
        ax.axis("off")
        fig.savefig(out_file, dpi=300)
        plt.close(fig)
        return

    # padj = 0 would be infinite on the log axis
    floor = d.loc[d["padj"] > 0, "padj"].min() if (d["padj"] > 0).any() else 1e-300
    d["y"] = -np.log10(d["padj"].clip(lower=floor / 10))
    sig = d["padj"] < PADJ_CUTOFF
    d["colour"] = np.where(~sig, "#b0b0b0", np.where(d["diff"] > 0, colour1, colour2))

    ax.scatter(d["diff"], d["y"], c=d["colour"], s=70, edgecolor="black", linewidth=0.5,
               alpha=POINT_ALPHA, zorder=3)
    ax.axhline(-np.log10(PADJ_CUTOFF), linestyle="--", color="black", linewidth=1)
    ax.axvline(0, linestyle="--", color="grey", linewidth=1)

    # symmetric x axis, some head room for the labels
    x_max = max(d["diff"].abs().max() * 1.25, 1e-12)
    ax.set_xlim(-x_max, x_max)
    ax.set_ylim(bottom=0, top=max(d["y"].max(), -np.log10(PADJ_CUTOFF)) * 1.2 + 0.2)

    # one label per SIGNIFICANT cluster, moved off its point to the first free spot (no overlap with
    # other labels or any point) and connected to the point by a thin line
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    limits = fig.bbox
    points_px = ax.transData.transform(d[["diff", "y"]].to_numpy())
    taken = [matplotlib.transforms.Bbox.from_extents(x - 6, y - 6, x + 6, y + 6) for x, y in points_px]
    offsets = [(r * np.cos(t), r * np.sin(t))
               for r in (14, 28, 45, 65, 90, 120, 155, 195, 240)
               for t in np.deg2rad([45, 135, 20, 160, 70, 110, 90, 0, 180, -20, -160, -45, -135])]
    for _, r in d[sig].sort_values("y", ascending=False).iterrows():
        label = cluster_label(r["cluster"], r["annotation"])
        placed = None
        for dx, dy in offsets:
            ann = ax.annotate(label, xy=(r["diff"], r["y"]), xytext=(dx, dy), textcoords="offset points",
                              fontsize=8, ha="left" if dx >= 0 else "right", va="center",
                              arrowprops=dict(arrowstyle="-", color="grey", linewidth=0.5), zorder=4)
            box = text_box(ann, renderer)
            inside = (box.x0 >= limits.x0 and box.x1 <= limits.x1 and
                      box.y0 >= ax.bbox.y0 and box.y1 <= limits.y1)
            if inside and not any(box.overlaps(t) for t in taken):
                placed = box
                break
            ann.remove()
        if placed is None:                       # no free spot: put it at the first position anyway
            ann = ax.annotate(label, xy=(r["diff"], r["y"]), xytext=offsets[0], textcoords="offset points",
                              fontsize=8, ha="left", va="center",
                              arrowprops=dict(arrowstyle="-", color="grey", linewidth=0.5), zorder=4)
            placed = text_box(ann, renderer)
        taken.append(placed)

    ax.set_xlabel(f"Difference of the mean {PROPERTIES[prop][0].lower() + PROPERTIES[prop][1:]}\n"
                  f"({short(biome1)} − {short(biome2)})")
    ax.set_ylabel("−log10(adjusted p-value)")
    ax.set_title(title)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", markerfacecolor=colour1, markeredgecolor="black",
                   markersize=9, alpha=POINT_ALPHA, label=f"higher in {short(biome1)} (padj < {PADJ_CUTOFF})"),
        plt.Line2D([], [], marker="o", linestyle="", markerfacecolor=colour2, markeredgecolor="black",
                   markersize=9, alpha=POINT_ALPHA, label=f"higher in {short(biome2)} (padj < {PADJ_CUTOFF})"),
        plt.Line2D([], [], marker="o", linestyle="", markerfacecolor="#b0b0b0", markeredgecolor="black",
                   markersize=9, alpha=POINT_ALPHA, label="not significant"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="upper center", ncol=3)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out_file, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("rf_results_dir", help="folder with <cluster>_data_df.tsv files")
    parser.add_argument("fasta_dir", help="folder with <cluster>.fasta / .faa / .fa files")
    parser.add_argument("--biome1", required=True, help="first biome (full label or the part after the last ';')")
    parser.add_argument("--biome2", required=True, help="second biome")
    parser.add_argument("--colour1", default="#0bd600")
    parser.add_argument("--colour2", default="#8e00fa")
    parser.add_argument("--out_dir", default=".")
    args = parser.parse_args()

    rf_dir, fasta_dir, out_dir = Path(args.rf_results_dir), Path(args.fasta_dir), Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_protein, stats_rows = [], []
    biome1 = biome2 = None
    for data_file in sorted(rf_dir.glob("*_data_df.tsv")):
        cluster = data_file.name.replace("_data_df.tsv", "")
        fasta_file = next((fasta_dir / f"{cluster}{ext}" for ext in (".fasta", ".faa", ".fa")
                           if (fasta_dir / f"{cluster}{ext}").exists()), None)
        if fasta_file is None:
            print(f"WARNING: no fasta for {cluster}, skipped", file=sys.stderr)
            continue

        props, annotation = protein_properties(fasta_file)
        biomes = pd.read_csv(data_file, sep="\t", usecols=["gene", "biome"]).drop_duplicates("gene").set_index("gene")
        biomes["biome"] = biomes["biome"].str.replace(":", ";")
        df = props.join(biomes, how="inner")      # only proteins with a biome label

        b1 = match_biome(args.biome1, df["biome"].unique())
        b2 = match_biome(args.biome2, df["biome"].unique())
        if b1 is None or b2 is None:
            print(f"WARNING: {cluster} does not contain both biomes "
                  f"(found: {sorted(df['biome'].unique())}), skipped", file=sys.stderr)
            continue
        biome1, biome2 = b1, b2
        df = df[df["biome"].isin([b1, b2])]

        out = df.rename_axis("gene").reset_index()
        out.insert(0, "cluster", cluster)
        per_protein.append(out[["cluster", "gene", "biome"] + list(PROPERTIES)])

        for prop in PROPERTIES:
            v1 = df.loc[df["biome"] == b1, prop].dropna()
            v2 = df.loc[df["biome"] == b2, prop].dropna()
            row = {"cluster": cluster, "annotation": annotation, "property": prop,
                   "n_biome1": len(v1), "n_biome2": len(v2),
                   "mean_biome1": v1.mean(), "mean_biome2": v2.mean(),
                   "diff": v1.mean() - v2.mean(), "U": np.nan, "p_value": np.nan}
            if len(v1) >= MIN_PER_BIOME and len(v2) >= MIN_PER_BIOME:
                row["U"], row["p_value"] = mannwhitneyu(v1, v2, alternative="two-sided")
            stats_rows.append(row)

    if not stats_rows:
        sys.exit("ERROR: no cluster with a data_df, a fasta and both biomes")

    pd.concat(per_protein).to_csv(out_dir / "protein_properties_per_protein.tsv", sep="\t", index=False)

    stats = pd.DataFrame(stats_rows)
    stats["padj"] = np.nan
    for prop in PROPERTIES:                      # BH over the clusters, separately for each property
        sel = stats["property"] == prop
        stats.loc[sel, "padj"] = bh_adjust(stats.loc[sel, "p_value"])
    stats.insert(2, "biome1", biome1)
    stats.insert(3, "biome2", biome2)
    stats.to_csv(out_dir / "protein_properties_cluster_stats.tsv", sep="\t", index=False)

    for prop in PROPERTIES:
        volcano(stats, prop, biome1, biome2, args.colour1, args.colour2, out_dir / f"volcano_{prop}.png")
        n_sig = int(((stats["property"] == prop) & (stats["padj"] < PADJ_CUTOFF)).sum())
        print(f"{prop}: {int((stats['property'] == prop).sum())} clusters, {n_sig} with padj < {PADJ_CUTOFF}",
              file=sys.stderr)


if __name__ == "__main__":
    main()