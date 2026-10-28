#!/usr/bin/env python3
"""
Create the per-cluster input table for phyloglm.

Inputs (both written by rf_model_v3.py):
  feature_df : <cluster>_Feature_importance.tsv   columns: feature, gini_imp, cluster_num
  data_df    : <cluster>_data_df.tsv               columns: gene, <one-hot MSA features...>, biome

Outputs:
  -o <cluster>_phyloglm_dataset.tsv : gene | ecosystem_subtype_<target> (0/1) | features (by gini, descending)
  <cluster>_gini.tsv                : feature | gini
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser(description="Create dataset for PhyloGLM")
parser.add_argument("feature_df", help="<cluster>_Feature_importance.tsv from rf_model_v3.py")
parser.add_argument("data_df", help="<cluster>_data_df.tsv from rf_model_v3.py")
parser.add_argument("-o", "--output", default="output_phyloglm_dataset.tsv", help="Output dataset path")
parser.add_argument("--target_biome", default=None,
                    help="Biome coded as 1 (full label or the part after the last ';'). "
                         "Default: the alphabetically last of the two biomes")
parser.add_argument("--gini_output", default=None,
                    help="Gini lookup path (default: derived from --output)")
args = parser.parse_args()


def short(label):
    return label.split(";")[-1].strip()


def safe(name):
    # column-name-safe version of a biome label, e.g. "Hot (42-90C)" -> "Hot_42_90C"
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


# ---------------------------------------------------------------
# Feature importances (already sorted by rf_model_v3.py; sort again to be safe)
# ---------------------------------------------------------------
feat = pd.read_csv(args.feature_df, sep="\t")
feat = feat.rename(columns={feat.columns[0]: "feature", feat.columns[1]: "gini"})[["feature", "gini"]]
feat["feature"] = feat["feature"].astype(str)
feat = feat.drop_duplicates("feature").sort_values("gini", ascending=False)

# ---------------------------------------------------------------
# Data table: one row per gene, one-hot features + biome label
# ---------------------------------------------------------------
data = pd.read_csv(args.data_df, sep="\t")
if "gene" not in data.columns:                      # index column written without a name
    data = data.rename(columns={data.columns[0]: "gene"})
if "biome" not in data.columns:
    sys.exit(f"ERROR: no 'biome' column in {args.data_df}. Columns start with: {list(data.columns[:5])}")

data = data.drop_duplicates("gene").set_index("gene")
data = data[data["biome"].notna() & (data["biome"] != "unknown")]

biomes = sorted(data["biome"].unique())
if len(biomes) != 2:
    sys.exit(f"ERROR: phyloglm needs exactly 2 biomes, found {len(biomes)}: "
             f"{data['biome'].value_counts().to_dict()}")

# Which biome is 1 (target) and which is 0 (reference)
if args.target_biome:
    matches = [b for b in biomes
               if b == args.target_biome or short(b) == short(args.target_biome)]
    if len(matches) != 1:
        sys.exit(f"ERROR: --target_biome '{args.target_biome}' matches {len(matches)} of {biomes}")
    target = matches[0]
else:
    target = biomes[1]
reference = biomes[0] if target == biomes[1] else biomes[1]
target_col = f"ecosystem_subtype_{safe(short(target))}"

# ---------------------------------------------------------------
# Keep features that exist in the data, in gini order
# ---------------------------------------------------------------
features = [f for f in feat["feature"] if f in data.columns]
n_missing = len(feat) - len(features)
if n_missing:
    print(f"WARNING: {n_missing} feature(s) from {args.feature_df} not found in {args.data_df}",
          file=sys.stderr)
if not features:
    sys.exit("ERROR: none of the features in the importance table are columns of the data table")

dataset = pd.DataFrame({target_col: (data["biome"] == target).astype(int)}, index=data.index)
dataset = dataset.join(data[features])
dataset.reset_index().to_csv(args.output, sep="\t", index=False)

# ---------------------------------------------------------------
# Gini lookup table (feature | gini), same features as the dataset
# ---------------------------------------------------------------
if args.gini_output:
    gini_out = Path(args.gini_output)
else:
    out = Path(args.output)
    stem = out.name.replace("_phyloglm_dataset.tsv", "")
    if stem == out.name:                            # output didn't follow the naming scheme
        stem = out.stem
    gini_out = out.with_name(f"{stem}_gini.tsv")

feat[feat["feature"].isin(features)].to_csv(gini_out, sep="\t", index=False)

counts = dataset[target_col].value_counts()
print(f"Target (1): {target}  n={counts.get(1, 0)}", file=sys.stderr)
print(f"Reference (0): {reference}  n={counts.get(0, 0)}", file=sys.stderr)
print(f"Features written: {len(features)}  -> {args.output}, {gini_out}", file=sys.stderr)