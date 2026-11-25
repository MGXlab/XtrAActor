#!/usr/bin/env python3
"""
Protein enrichment report from pyseer results.

Significance: Bonferroni over the number of distinct presence/absence patterns tested
(clusters with identical patterns count once), as recommended by pyseer:
    threshold = alpha / n_unique_patterns
A Benjamini-Hochberg adjusted p-value over all tested clusters is added for reference.

Direction: biome2 is phenotype 1, so beta > 0 = enriched in biome2, beta < 0 = enriched in biome1.

Outputs:
  enrichment_all.tsv      every tested cluster
  enrichment_report.tsv   significant clusters only
  enrichment_summary.txt  counts

Usage:
  enrichment_report.py --pyseer pyseer.assoc --patterns patterns.txt --info cluster_info.tsv \\
      --biome1 "..." --biome2 "..." [--alpha 0.05]
"""
import argparse
import sys

import numpy as np
import pandas as pd


def short(label):
    return label.split(";")[-1].strip()


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = ok.sum()
    if m:
        order = np.argsort(p[ok])
        ranked = p[ok][order] * m / np.arange(1, m + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        tmp = np.empty(m)
        tmp[order] = np.minimum(ranked, 1.0)
        out[ok] = tmp
    return out


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--pyseer", required=True, help="pyseer output (tab-separated)")
parser.add_argument("--patterns", required=True, help="pyseer --output-patterns file")
parser.add_argument("--info", required=True, help="cluster_info.tsv from enrichment_inputs.py presence")
parser.add_argument("--biome1", required=True)
parser.add_argument("--biome2", required=True)
parser.add_argument("--alpha", type=float, default=0.05)
parser.add_argument("--out_all", default="enrichment_all.tsv")
parser.add_argument("--out_report", default="enrichment_report.tsv")
parser.add_argument("--out_summary", default="enrichment_summary.txt")
a = parser.parse_args()

res = pd.read_csv(a.pyseer, sep="\t")
info = pd.read_csv(a.info, sep="\t")
with open(a.patterns) as fh:
    n_patterns = len({l.strip() for l in fh if l.strip()})
threshold = a.alpha / n_patterns if n_patterns else np.nan

res = res.rename(columns={"variant": "cluster", "lrt-pvalue": "p_value", "beta-std-err": "beta_se",
                          "filter-pvalue": "filter_p_value"})
df = info.merge(res, on="cluster", how="inner")
df["padj_bh"] = bh_adjust(df["p_value"])
df["significant"] = df["p_value"] < threshold
b1, b2 = short(a.biome1), short(a.biome2)
df["enriched_in"] = np.where(df["beta"] > 0, b2, b1)
df.loc[df["beta"].isna(), "enriched_in"] = ""
df = df.sort_values(["p_value", "cluster"], na_position="last")

first = ["cluster", "annotation", "enriched_in", "beta", "beta_se", "p_value", "padj_bh", "significant",
         "af", "n_genomes", "n_genomes_biome1", "n_genomes_biome2", "single_biome", "in_rf_analysis",
         "n_proteins", "n_proteins_with_derep", "representative"]
cols = [c for c in first if c in df.columns] + [c for c in df.columns if c not in first]
df = df[cols]
for c in ["beta", "beta_se", "p_value", "padj_bh", "filter_p_value", "af"]:
    if c in df.columns:
        df[c] = df[c].apply(lambda x: float(f"{x:.4g}") if pd.notna(x) else x)

df.to_csv(a.out_all, sep="\t", index=False)
sig = df[df["significant"]]
sig.to_csv(a.out_report, sep="\t", index=False)

n_tested = int(df["p_value"].notna().sum())
lines = [
    f"clusters in the presence/absence matrix: {len(info)}",
    f"clusters tested by pyseer (after its frequency filter): {n_tested}",
    f"distinct presence/absence patterns: {n_patterns}",
    f"significance threshold: p < {a.alpha} / {n_patterns} = {threshold:.3g}",
    f"significant clusters: {len(sig)}",
    f"  enriched in {b1} (biome1): {int((sig['enriched_in'] == b1).sum())}",
    f"  enriched in {b2} (biome2): {int((sig['enriched_in'] == b2).sum())}",
    f"  of these single-biome: {int(sig['single_biome'].astype(str).eq('True').sum())}",
    f"  of these also in the RF analysis: {int(sig['in_rf_analysis'].astype(str).eq('True').sum())}",
    f"clusters with BH padj < {a.alpha} (for reference): {int((df['padj_bh'] < a.alpha).sum())}",
]
if "notes" in df.columns:
    notes = df["notes"].dropna().astype(str).str.split(",").explode().value_counts()
    if len(notes):
        lines.append("pyseer notes: " + ", ".join(f"{k} ({v})" for k, v in notes.items()))
open(a.out_summary, "w").write("\n".join(lines) + "\n")
print("\n".join(lines), file=sys.stderr)
