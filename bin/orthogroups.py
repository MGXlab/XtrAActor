import argparse, re, sys
from collections import Counter
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import matplotlib
matplotlib.use("Agg")            # no display on the cluster nodes
import matplotlib.pyplot as plt

MIN_PER_BIOME = 10   # a cluster is written only if EACH biome in the metadata has at least this many members

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

# ---------------------------------------------------------------
# Cluster size histogram (all clusters, before the biome filter)
# ---------------------------------------------------------------
cluster_sizes = clusters.groupby("representative").size()         # members per cluster
size_counts = cluster_sizes.value_counts().sort_index()           # size -> number of clusters

fig, ax = plt.subplots(figsize=(8, 5))
ax.vlines(size_counts.index, 0, size_counts.values, color="#1f77b4", linewidth=1.5)
ax.scatter(size_counts.index, size_counts.values, color="#1f77b4", s=8, zorder=3)

min_pass = MIN_PER_BIOME * len(biomes)                            # smallest size that can pass
ax.axvline(min_pass, color="#d62728", linestyle="--", linewidth=1,
           label=f"min. size to pass filter ({len(biomes)} × {MIN_PER_BIOME} = {min_pass})")

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
    # Pharokka: EAF_Gu_TropicPondC2c0115_c1_CDS_0001 -> EAF_Gu_TropicPondC2c0115_c1
    if re.search(r'_CDS_\d+$', member):
        return re.sub(r'_CDS_\d+$', '', member)
    # Prodigal-style / renamed: <id>_<n> -> <id>
    return member.rsplit('_', 1)[0]


n_clusters = 0
n_written = 0
n_missing_seqs = 0
n_unknown_biome = 0

for i, (rep, group) in enumerate(clusters.groupby("representative")):
    n_clusters += 1

    members = [m for m in group["member"] if m in seqs]
    n_missing_seqs += len(group) - len(members)

    # Biome mapping
    rows = []
    for member in members:
        genome = extract_genome(member)
        biome = biome_lookup.get(genome, "unknown")
        if biome == "unknown":
            n_unknown_biome += 1
        rows.append({"gene": member, "representative": rep, "genome": genome, "biome": biome})

    # Filter: every metadata biome needs at least MIN_PER_BIOME members
    counts = Counter(r["biome"] for r in rows)
    if any(counts[b] < MIN_PER_BIOME for b in biomes):
        continue

    SeqIO.write((seqs[m] for m in members), fasta_out_dir / f"cluster_{i}.fasta", "fasta")
    pd.DataFrame(rows).to_csv(biome_out_dir / f"cluster_{i}_biome.tsv", sep="\t", index=False)
    n_written += 1

print(f"Clusters total: {n_clusters}", file=sys.stderr)
print(f"Clusters written (>= {MIN_PER_BIOME} members per biome): {n_written}", file=sys.stderr)
print(f"Members missing from FASTA: {n_missing_seqs}", file=sys.stderr)
print(f"Genes with no biome match: {n_unknown_biome}", file=sys.stderr)

if n_written == 0:
    sys.exit("ERROR: no cluster passed the biome filter")