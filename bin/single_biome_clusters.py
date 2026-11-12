#!/usr/bin/env python3
"""
Select the "single-biome clusters" of the MMseqs2 50/50 clustering: clusters whose
members all come from genomes of ONE biome. (orthogroups_sel_by_biome.py keeps the
clusters with enough members in EVERY biome; this script keeps the opposite set.)

Inputs
  cluster_tsv    mmseq2_cluster.tsv            representative <TAB> member (no header)
  metadata_file  genome <whitespace> biome     (header line; same file as for orthogroups_sel_by_biome.py)
  --derep_dir    folder with <biome_set>_derep_cluster.tsv files (optional). The 50/50 clustering
                 runs on dereplicated proteins, so a genome whose protein was removed as a
                 near-identical copy would look like it lacks the cluster. With these files every
                 cluster member is expanded back to all proteins it represents.
  --proteins_fasta  all predicted proteins (optional); only used for the annotation column

Output (-o): one row per single-biome cluster found in at least --min_genomes genomes
  cluster | representative | biome | n_genomes | n_proteins | annotation | genomes (comma separated)

'cluster' is cluster_<i> with the same numbering as orthogroups_sel_by_biome.py.
"""
import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

UNINFORMATIVE = ("hypothetical", "unknown function", "uncharacterized", "uncharacterised")

parser = argparse.ArgumentParser()
parser.add_argument("cluster_tsv")
parser.add_argument("metadata_file")
parser.add_argument("--derep_dir", default=None)
parser.add_argument("--proteins_fasta", default=None)
parser.add_argument("--min_genomes", type=int, default=3,
                    help="minimum number of genomes carrying the cluster (default: 3)")
parser.add_argument("-o", "--output", default="single_biome_clusters.tsv")
args = parser.parse_args()


def extract_genome(member):
    # identical to orthogroups_sel_by_biome.py
    if re.search(r'_CDS_\d+$', member):
        contig = re.sub(r'_CDS_\d+$', '', member)
    else:
        contig = member.rsplit('_', 1)[0]
    return contig.split('|')[0]


# Metadata: split only at the first whitespace, since biome names contain spaces
biome_lookup = {}
with open(args.metadata_file) as fh:
    next(fh)
    for line in fh:
        parts = line.strip().split(None, 1)
        if len(parts) == 2:
            biome_lookup.setdefault(parts[0], parts[1].strip())
biomes = sorted(set(biome_lookup.values()))
print(f"Biomes in metadata: {biomes}", file=sys.stderr)

# Dereplication map: representative -> all proteins it stands for
derep = defaultdict(list)
if args.derep_dir and Path(args.derep_dir).is_dir():
    for f in sorted(Path(args.derep_dir).glob("*.tsv")):
        with open(f) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0] != parts[1]:
                    derep[parts[0]].append(parts[1])
print(f"Dereplicated proteins added back: {sum(len(v) for v in derep.values())}", file=sys.stderr)

# Clusters: representative -> member proteins (expanded)
clusters = defaultdict(list)
with open(args.cluster_tsv) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        clusters[parts[0]].append(parts[1])
        clusters[parts[0]].extend(derep.get(parts[1], []))

n_multi = n_unknown = n_small = 0
kept = []
for i, rep in enumerate(sorted(clusters)):          # same order as pandas groupby("representative")
    members = clusters[rep]
    genomes = {extract_genome(m) for m in members}
    member_biomes = {biome_lookup.get(g, "unknown") for g in genomes}
    if "unknown" in member_biomes:                  # cannot be called single-biome
        n_unknown += 1
        continue
    if len(member_biomes) != 1:
        n_multi += 1
        continue
    if len(genomes) < args.min_genomes:
        n_small += 1
        continue
    kept.append({"cluster": f"cluster_{i}", "representative": rep, "biome": member_biomes.pop(),
                 "genomes": sorted(genomes), "members": set(members)})

# Annotation = most frequent informative header text among the cluster's proteins
annotation = {}
if args.proteins_fasta and kept:
    prot2idx = {m: k for k, c in enumerate(kept) for m in c["members"]}
    counts = [Counter() for _ in kept]
    with open(args.proteins_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                pid, _, desc = line[1:].strip().partition(" ")
                k = prot2idx.get(pid)
                desc = desc.strip()
                if k is not None and desc and not any(u in desc.lower() for u in UNINFORMATIVE):
                    counts[k][desc] += 1
    annotation = {k: c.most_common(1)[0][0] for k, c in enumerate(counts) if c}

with open(args.output, "w") as out:
    out.write("cluster\trepresentative\tbiome\tn_genomes\tn_proteins\tannotation\tgenomes\n")
    for k, c in enumerate(kept):
        out.write("\t".join([c["cluster"], c["representative"], c["biome"], str(len(c["genomes"])),
                             str(len(c["members"])), annotation.get(k, "hypothetical protein"),
                             ",".join(c["genomes"])]) + "\n")

print(f"Clusters total: {len(clusters)}", file=sys.stderr)
print(f"Clusters with members from more than one biome: {n_multi}", file=sys.stderr)
print(f"Clusters with members of unknown biome (skipped): {n_unknown}", file=sys.stderr)
print(f"Single-biome clusters in < {args.min_genomes} genomes (skipped): {n_small}", file=sys.stderr)
print(f"Single-biome clusters written: {len(kept)}", file=sys.stderr)
for b, n in sorted(Counter(c["biome"] for c in kept).items()):
    print(f"  {b}: {n}", file=sys.stderr)
