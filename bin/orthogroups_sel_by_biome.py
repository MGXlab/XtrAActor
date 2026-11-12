import argparse, re, sys
from collections import Counter
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import matplotlib
matplotlib.use("Agg")            # no display on the cluster nodes
import matplotlib.pyplot as plt

MIN_PER_BIOME = 30      # a cluster is written for RF only if EACH biome in the metadata has at least this many members
MIN_SINGLE_BIOME = 30   # a single-biome cluster is reported only if its one biome has at least this many members
SINGLE_BIOME_OUT = "single_biome_clusters.tsv"

UNINFORMATIVE = ("hypothetical", "unknown function", "uncharacterized", "uncharacterised")

parser = argparse.ArgumentParser()
parser.add_argument("cluster_tsv")
parser.add_argument("proteins_fasta")
parser.add_argument("metadata_file")
parser.add_argument("--fasta_out_dir", default="cluster_fastas")
parser.add_argument("--biome_out_dir", default="cluster_biomes")
args = parser.parse_args()

fasta_out_dir = Path(args.fasta_out_dir)
biome_out_dir = Path(args.biome_out_dir)
fasta_out_dir.mkdir(parents=True, exist_ok=True)
biome_out_dir.mkdir(parents=True, exist_ok=True)

clusters = pd.read_csv(args.cluster_tsv, sep="\t", names=["representative", "member"])
seqs = SeqIO.index(args.proteins_fasta, "fasta")  # indexes on disk instead of loading all into memory

# Metadata: split only at the first whitespace, since biome names contain spaces
biome_lookup = {}
with open(args.metadata_file) as fh:
    next(fh)  # skip header line
    for line in fh:
        parts = line.strip().split(None, 1)
        if len(parts) == 2:
            biome_lookup.setdefault(parts[0], parts[1].strip())  # keep first entry per genome

# Biomes to compare = the distinct biomes in the metadata
biomes = sorted(set(biome_lookup.values()))
print(f"Biomes in metadata: {biomes}", file=sys.stderr)
if len(biomes) != 2:
    print(f"WARNING: expected 2 biomes in metadata, found {len(biomes)}; "
          f"each of them must reach {MIN_PER_BIOME} members", file=sys.stderr)

min_size = MIN_PER_BIOME * len(biomes)        # smallest cluster size that can pass the RF filter
min_check = min(min_size, MIN_SINGLE_BIOME)   # smallest cluster size worth mapping at all

# ---------------------------------------------------------------
# Cluster size histogram (all clusters, before any filter)
# ---------------------------------------------------------------
cluster_sizes = clusters.groupby("representative").size()         # members per cluster
size_counts = cluster_sizes.value_counts().sort_index()           # size -> number of clusters

fig, ax = plt.subplots(figsize=(8, 5))
ax.vlines(size_counts.index, 0, size_counts.values, color="#002df5", linewidth=1.5)
ax.scatter(size_counts.index, size_counts.values, color="#002df5", s=8, zorder=3)

ax.axvline(min_size, color="#d62728", linestyle="--", linewidth=1,
           label=f"min. size to pass filter ({len(biomes)} × {MIN_PER_BIOME} = {min_size})")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(bottom=0.8)
ax.set_xlabel("Cluster size (number of sequences)")
ax.set_ylabel("Number of clusters")
ax.set_title(f"Cluster size distribution ({len(cluster_sizes):,} clusters, "
             f"{int(cluster_sizes.sum()):,} sequences)")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig("cluster_size_histogram.png", dpi=200)
plt.close(fig)

print(f"Cluster sizes: min={cluster_sizes.min()}, median={cluster_sizes.median():.0f}, "
      f"max={cluster_sizes.max()}, singletons={int((cluster_sizes == 1).sum())}", file=sys.stderr)
# ---------------------------------------------------------------


def extract_genome(member):
    # strip the CDS suffix: <contig>_CDS_0003 -> <contig>   (Pharokka, and Bakta after renaming)
    if re.search(r'_CDS_\d+$', member):
        contig = re.sub(r'_CDS_\d+$', '', member)
    else:
        # Prodigal-style / renamed: <id>_<n> -> <id>
        contig = member.rsplit('_', 1)[0]
    # IMG/VR and assembly headers carry extra '|' fields: the metadata uses only the part before the first '|'
    #   IMGVR_UViG_3300032238_001132|3300032238|Ga0335266_1155248 -> IMGVR_UViG_3300032238_001132
    #   NODE_231_length_74457_cov_12.8046||0_partial               -> NODE_231_length_74457_cov_12.8046
    return contig.split('|')[0]


def cluster_annotation(members):
    """Most frequent informative annotation (text after the ID in the fasta header).
    If no member has one, the most frequent uninformative label instead
    (e.g. 'hypothetical protein', 'unknown function'); 'hypothetical protein' if headers are empty."""
    informative = Counter()
    uninformative = Counter()
    for m in members:
        rec = seqs[m]
        desc = rec.description[len(rec.id):].strip()
        d = desc.lower()
        if d and not any(u in d for u in UNINFORMATIVE):
            informative[desc] += 1
        elif d:
            uninformative[desc] += 1

    n_annotated = sum(informative.values())
    if informative:
        top, top_n = informative.most_common(1)[0]
        return top, n_annotated, round(top_n / n_annotated, 3)
    if uninformative:
        return uninformative.most_common(1)[0][0], 0, 0.0
    return "hypothetical protein", 0, 0.0


n_clusters = 0
n_checked = 0
n_written = 0
n_missing_seqs = 0
n_unknown_biome = 0
single_biome_rows = []

for i, (rep, group) in enumerate(clusters.groupby("representative")):
    n_clusters += 1

    members = [m for m in group["member"] if m in seqs]
    n_missing_seqs += len(group) - len(members)

    # Quick skip: too small for either filter (saves the mapping for ~all singletons)
    if len(members) < min_check:
        continue
    n_checked += 1

    # Biome mapping
    rows = []
    for member in members:
        genome = extract_genome(member)
        biome = biome_lookup.get(genome, "unknown")
        if biome == "unknown":
            n_unknown_biome += 1
        rows.append({"gene": member, "representative": rep, "genome": genome, "biome": biome})

    counts = Counter(r["biome"] for r in rows)
    known = {b: n for b, n in counts.items() if b != "unknown"}

    # --- single-biome clusters: all members with a known biome come from ONE biome,
    #     or no member has a known biome at all (reported as biome 'unknown') ---
    if len(known) <= 1:
        if known:
            only_biome, n_biome = next(iter(known.items()))
        else:
            only_biome, n_biome = "unknown", counts.get("unknown", 0)
        if n_biome >= MIN_SINGLE_BIOME:
            annotation, n_annotated, share = cluster_annotation(members)
            single_biome_rows.append({
                "cluster": f"cluster_{i}",
                "annotation": annotation,
                "biome": only_biome,
                "n_members": len(members),
                "n_biome": n_biome,
                "n_unknown": counts.get("unknown", 0),
                "n_annotated": n_annotated,
                "annotation_share": share,
            })
        continue   # one biome (or none) -> can never pass the per-biome filter for RF

    # --- per-biome filter for RF: every metadata biome needs at least MIN_PER_BIOME members ---
    if len(members) < min_size or any(counts[b] < MIN_PER_BIOME for b in biomes):
        continue

    SeqIO.write((seqs[m] for m in members), fasta_out_dir / f"cluster_{i}.fasta", "fasta")
    pd.DataFrame(rows).to_csv(biome_out_dir / f"cluster_{i}_biome.tsv", sep="\t", index=False)
    n_written += 1

# --- single-biome report (always written, header only if none found) ---
single_cols = ["cluster", "annotation", "biome", "n_members", "n_biome",
               "n_unknown", "n_annotated", "annotation_share"]
single_df = pd.DataFrame(single_biome_rows, columns=single_cols)
single_df = single_df.sort_values(["biome", "n_biome"], ascending=[True, False])
single_df.to_csv(SINGLE_BIOME_OUT, sep="\t", index=False)

print(f"Clusters total: {n_clusters}", file=sys.stderr)
print(f"Clusters checked for biomes (>= {min_check} members): {n_checked}", file=sys.stderr)
print(f"Clusters written (>= {MIN_PER_BIOME} members per biome): {n_written}", file=sys.stderr)
print(f"Single-biome clusters (>= {MIN_SINGLE_BIOME} members of one biome): {len(single_df)}", file=sys.stderr)
for b, n in single_df["biome"].value_counts().items():
    print(f"    {b}: {n}", file=sys.stderr)
print(f"Members missing from FASTA: {n_missing_seqs}", file=sys.stderr)
print(f"Genes with no biome match (in checked clusters): {n_unknown_biome}", file=sys.stderr)

if n_written == 0:
    sys.exit("ERROR: no cluster passed the biome filter")