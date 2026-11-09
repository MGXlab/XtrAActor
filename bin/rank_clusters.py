#!/usr/bin/env python3
"""
Rank clusters by RF accuracy (highest first) and attach a functional annotation
taken from the headers of each cluster fasta.

One output row per cluster AND biome (the RF report has one row per biome):
  cluster | annotation | biome | accuracy | accuracy_std | precision | recall | f1 | support | ...

Header format of the cluster fastas: ><protein_id> <annotation>, e.g.
  >NODE_5731_length_13416_cov_14.9073||full_CDS_0010 glycoside hydrolase
  >IMGVR_UViG_..._CDS_0010 hypothetical protein

Annotation per cluster = the most frequent informative annotation among its sequences;
'hypothetical protein' only if no sequence has a real annotation.

Usage:
  rank_clusters.py <reports_dir> <fasta_dir> [-o cluster_ranking.tsv]
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO

UNINFORMATIVE = ("hypothetical", "unknown function", "uncharacterized", "uncharacterised")
METRICS = ["accuracy", "accuracy_std", "precision", "recall", "f1", "support"]


def is_informative(desc):
    d = desc.strip().lower()
    return bool(d) and not any(u in d for u in UNINFORMATIVE)


def cluster_annotation(fasta):
    n_seqs = 0
    annotations = Counter()
    for rec in SeqIO.parse(fasta, "fasta"):
        n_seqs += 1
        desc = rec.description[len(rec.id):].strip()   # text after the ID
        if is_informative(desc):
            annotations[desc] += 1

    n_annotated = sum(annotations.values())
    if annotations:
        top, top_n = annotations.most_common(1)[0]
        share = top_n / n_annotated
    else:
        top, share = "hypothetical protein", 0.0

    other = "; ".join(f"{a} ({n})" for a, n in annotations.most_common(4)[1:])
    return {
        "annotation": top,
        "n_sequences": n_seqs,
        "n_annotated": n_annotated,
        "annotation_share": round(share, 3),   # fraction of annotated sequences with the top annotation
        "other_annotations": other,
    }


def read_report(report):
    """All rows (one per biome) of one RF report, with biome + metric columns."""
    df = pd.read_csv(report, sep="\t")      # quoted multi-line fields (RF_parameters, conf_matrix) are fine
    missing = [c for c in ["biome"] + METRICS if c not in df.columns]
    if missing:
        raise ValueError(f"missing columns {missing}")
    return df[["biome"] + METRICS]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reports_dir", help="folder with <cluster>_RandomForest_Report.tsv files")
    parser.add_argument("fasta_dir", help="folder with <cluster>.fasta / .faa files")
    parser.add_argument("-o", "--output", default="cluster_ranking.tsv")
    args = parser.parse_args()

    reports_dir, fasta_dir = Path(args.reports_dir), Path(args.fasta_dir)

    # RF metrics per cluster (cluster name from the file name)
    reports = {}
    for f in reports_dir.glob("*_RandomForest_Report.tsv"):
        name = f.name.replace("_RandomForest_Report.tsv", "")
        try:
            reports[name] = read_report(f)
        except Exception as e:
            print(f"WARNING: could not read {f.name}: {e}", file=sys.stderr)

    # one row per cluster and biome
    rows = []
    fastas = [f for ext in ("*.fasta", "*.faa", "*.fa") for f in fasta_dir.glob(ext)]
    for f in fastas:
        name = f.stem
        ann = cluster_annotation(f)
        if name in reports:
            for _, r in reports[name].iterrows():
                rows.append({"cluster": name, "annotation": ann["annotation"], **r.to_dict(),
                             **{k: v for k, v in ann.items() if k != "annotation"}})
        else:
            # no RF result for this cluster: one row with empty metrics
            rows.append({"cluster": name, "annotation": ann["annotation"], "biome": None,
                         **{m: None for m in METRICS},
                         **{k: v for k, v in ann.items() if k != "annotation"}})

    if not rows:
        sys.exit("ERROR: no cluster fasta files found")

    ranking = pd.DataFrame(rows)
    ranking[METRICS] = ranking[METRICS].apply(pd.to_numeric, errors="coerce")

    # sort: clusters by accuracy (highest first); the biome rows of a cluster stay together
    ranking = ranking.sort_values(["accuracy", "cluster", "biome"],
                                  ascending=[False, True, True], na_position="last")

    for m in ["accuracy", "accuracy_std", "precision", "recall", "f1"]:
        ranking[m] = ranking[m].round(3)
    ranking["support"] = ranking["support"].astype("Int64")   # 173.0 -> 173, empty stays empty

    ranking = ranking[["cluster", "annotation", "biome"] + METRICS +
                      ["n_sequences", "n_annotated", "annotation_share", "other_annotations"]]
    ranking.to_csv(args.output, sep="\t", index=False)

    n_clusters = ranking["cluster"].nunique()
    n_no_rf = ranking.loc[ranking["accuracy"].isna(), "cluster"].nunique()
    print(f"Clusters ranked: {n_clusters} ({n_no_rf} without RF result, listed last)", file=sys.stderr)
    print("Top 10 clusters:", file=sys.stderr)
    for cluster in ranking["cluster"].drop_duplicates().head(10):
        sub = ranking[ranking["cluster"] == cluster]
        per_biome = ", ".join(
            f"{str(r['biome']).split(';')[-1]}: P={r['precision']} R={r['recall']} F1={r['f1']} n={r['support']}"
            for _, r in sub.iterrows())
        print(f"{cluster} | {sub['annotation'].iloc[0]} | acc={sub['accuracy'].iloc[0]} | {per_biome}",
              file=sys.stderr)


if __name__ == "__main__":
    main()