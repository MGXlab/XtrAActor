#!/usr/bin/env python3
"""
Inputs for the protein enrichment analysis (pyseer).

Subcommands
-----------
presence   MMseqs2 clusters -> presence/absence matrix (clusters x genomes, pyseer --pres),
           phenotype file (genome -> biome as 0/1) and a cluster info table.
genomes    genome fastas -> one fasta with one record per genome, named by genome ID
           (input for ViPTreeGen and mash -i), plus a table of the names used in it.
mash       `mash dist` output (stdin) -> square distance matrix named by genome ID.
viptree    ViPTreeGen all.dist.matrix -> square distance matrix named by genome ID.

Genome IDs
----------
  protein ID  <genome>_CDS_<n>    -> genome = part before "_CDS_<n>", then before the first '|'
                                     (same rule as orthogroups_sel_by_biome.py)
  fasta file  <genome>.fna[.gz]   -> genome = file name without extension (non-phage: one file per genome)
  contig      ><genome>|... desc  -> genome = header up to the first whitespace, then before the
                                     first '|' (phage: one contig per genome)

All matrices are tab-separated with genome IDs as header row and first column, which is the
format pyseer reads with --distances / --similarity.
"""
import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

UNINFORMATIVE = ("hypothetical", "unknown function", "uncharacterized", "uncharacterised")


def extract_genome(member):
    """Same rule as orthogroups_sel_by_biome.py."""
    if re.search(r'_CDS_\d+$', member):
        contig = re.sub(r'_CDS_\d+$', '', member)
    else:
        contig = member.rsplit('_', 1)[0]
    return contig.split('|')[0]


def short(label):
    return label.split(";")[-1].strip()


def read_biomes(metadata_file):
    """genome -> biome; header line skipped, split at the first whitespace (biome names contain spaces)."""
    lookup = {}
    with open(metadata_file) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                lookup.setdefault(parts[0], parts[1].strip())
    return lookup


def match_biome(wanted, present):
    hits = [b for b in present if b == wanted or short(b) == short(wanted)]
    if len(hits) != 1:
        sys.exit(f"ERROR: biome '{wanted}' matches {len(hits)} of the metadata biomes {sorted(present)}")
    return hits[0]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def genome_id_from_file(path):
    name = Path(path).name
    name = re.sub(r"\.gz$", "", name)
    return re.sub(r"\.[^.]+$", "", name)


# =====================================================================
# presence
# =====================================================================
def run_presence(a):
    genome2biome = read_biomes(a.metadata)
    biomes = set(genome2biome.values())
    b1, b2 = match_biome(a.biome1, biomes), match_biome(a.biome2, biomes)
    print(f"biome1 (phenotype 0): {b1}\nbiome2 (phenotype 1): {b2}", file=sys.stderr)

    # dereplication: representative -> the proteins it stands for
    derep = defaultdict(list)
    if a.derep_dir and Path(a.derep_dir).is_dir():
        for f in sorted(Path(a.derep_dir).glob("*.tsv")):
            with open(f) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 2 and p[0] != p[1]:
                        derep[p[0]].append(p[1])

    clusters = defaultdict(list)
    with open(a.clusters) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                clusters[p[0]].append(p[1])

    rf_clusters = set()
    if a.rf_clusters and Path(a.rf_clusters).exists():
        rf_clusters = {l.strip() for l in open(a.rf_clusters) if l.strip()}

    # samples: genomes of the two biomes that carry at least one clustered protein
    kept, carriers = [], {}
    all_genomes = set()
    for i, rep in enumerate(sorted(clusters)):      # same numbering as orthogroups_sel_by_biome.py
        members = clusters[rep]
        expanded = members + [d for m in members for d in derep.get(m, [])]
        genomes = {extract_genome(m) for m in expanded}
        all_genomes |= genomes
        if len(members) < a.min_proteins:
            continue
        name = f"cluster_{i}"
        kept.append((name, rep, len(members), len(expanded), set(members)))
        carriers[name] = genomes

    samples = sorted(g for g in all_genomes if genome2biome.get(g) in (b1, b2))
    n_unknown = len([g for g in all_genomes if g not in genome2biome])
    print(f"genomes with clustered proteins: {len(all_genomes)} | in biome1/biome2: {len(samples)} | "
          f"not in metadata: {n_unknown}", file=sys.stderr)
    if n_unknown:
        ex = sorted(g for g in all_genomes if g not in genome2biome)[:3]
        print(f"  e.g. {ex} (check that the metadata genome IDs match the protein IDs)", file=sys.stderr)
    if len(samples) < 10:
        sys.exit("ERROR: fewer than 10 genomes with a biome; check the metadata genome IDs")

    # annotation per cluster: most frequent informative header text among its proteins
    annotation = {}
    if a.proteins_fasta:
        prot2cl = {m: name for name, _, _, _, mem in kept for m in mem}
        counts = defaultdict(Counter)
        with open_text(a.proteins_fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    pid, _, desc = line[1:].strip().partition(" ")
                    name = prot2cl.get(pid)
                    desc = desc.strip()
                    if name and desc and not any(u in desc.lower() for u in UNINFORMATIVE):
                        counts[name][desc] += 1
        annotation = {n: c.most_common(1)[0][0] for n, c in counts.items() if c}

    sample_set = set(samples)
    with open(a.out_pres, "w") as pres, open(a.out_info, "w") as info:
        pres.write("Gene\t" + "\t".join(samples) + "\n")
        info.write("cluster\trepresentative\tannotation\tn_proteins\tn_proteins_with_derep\t"
                   "n_genomes\tn_genomes_biome1\tn_genomes_biome2\tsingle_biome\tin_rf_analysis\n")
        n_written = 0
        for name, rep, n_prot, n_exp, _ in kept:
            g = carriers[name] & sample_set
            if not g:
                continue
            n1 = sum(genome2biome[x] == b1 for x in g)
            n2 = len(g) - n1
            pres.write(name + "\t" + "\t".join("1" if s in g else "0" for s in samples) + "\n")
            info.write(f"{name}\t{rep}\t{annotation.get(name, 'hypothetical protein')}\t{n_prot}\t{n_exp}\t"
                       f"{len(g)}\t{n1}\t{n2}\t{n1 == 0 or n2 == 0}\t{name in rf_clusters}\n")
            n_written += 1

    with open(a.out_pheno, "w") as fh:
        fh.write("samples\tbiome\n")
        for s in samples:
            fh.write(f"{s}\t{1 if genome2biome[s] == b2 else 0}\n")

    n1 = sum(genome2biome[s] == b1 for s in samples)
    print(f"clusters total: {len(clusters)} | with >= {a.min_proteins} proteins: {len(kept)} | "
          f"written: {n_written}", file=sys.stderr)
    print(f"phenotype: {n1} x 0 ({short(b1)}), {len(samples) - n1} x 1 ({short(b2)})", file=sys.stderr)


# =====================================================================
# genomes (ViPTreeGen input)
# =====================================================================
def viptree_name(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def run_genomes(a):
    """One record per genome, named by genome ID (ViPTreeGen / mash -i input).
    phage (default): every contig is one genome, ID = header before the first whitespace and '|'
    --one_per_file:  every file is one genome, ID = file name; its contigs are joined with 100 N"""
    seen = {}
    with open(a.out, "w") as out:
        for path in a.input:
            if a.one_per_file:
                gid = genome_id_from_file(path)
                if gid in seen:
                    sys.exit(f"ERROR: two files give the genome ID {gid}")
                seen[gid] = viptree_name(gid)
                out.write(f">{seen[gid]}\n")
                first = True
                with open_text(path) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            if not first:
                                out.write("N" * 100 + "\n")
                            first = False
                        else:
                            out.write(line)
                continue
            with open_text(path) as fh:
                write = False
                for line in fh:
                    if line.startswith(">"):
                        gid = line[1:].split()[0].split("|")[0]
                        write = gid not in seen
                        if not write:
                            print(f"WARNING: duplicate genome {gid} in {path}, skipped", file=sys.stderr)
                            continue
                        seen[gid] = viptree_name(gid)
                        out.write(f">{seen[gid]}\n")
                    elif write:
                        out.write(line)
    clash = [v for v, n in Counter(seen.values()).items() if n > 1]
    if clash:
        sys.exit(f"ERROR: {len(clash)} genome IDs become identical after renaming, e.g. {clash[:3]}")
    with open(a.names, "w") as names:
        names.write("genome\tviptree_name\n")
        for gid, v in seen.items():
            names.write(f"{gid}\t{v}\n")
    print(f"{len(seen)} genomes from {len(a.input)} file(s) written to {a.out}", file=sys.stderr)


# =====================================================================
# square distance matrices
# =====================================================================
def write_square(names, get, out_file):
    import numpy as np
    n = len(names)
    with open(out_file, "w") as fh:
        fh.write("\t" + "\t".join(names) + "\n")
        for i in range(n):
            row = get(i)
            fh.write(names[i] + "\t" + "\t".join(f"{x:.6g}" for x in row) + "\n")


def run_mash(a):
    import numpy as np
    v2g = {}
    if a.names:
        with open(a.names) as fh:
            next(fh)
            for line in fh:
                g, v = line.rstrip("\n").split("\t")
                v2g[v] = g
    idx, pairs = {}, []
    src = open(a.input) if a.input != "-" else sys.stdin
    for line in src:
        p = line.split("\t")
        if len(p) < 3:
            continue
        if v2g:          # sequences of genomes.fna (mash sketch -i)
            g1, g2 = v2g[p[0]], v2g[p[1]]
        else:            # one file per genome
            g1, g2 = genome_id_from_file(p[0]), genome_id_from_file(p[1])
        for g in (g1, g2):
            if g not in idx:
                idx[g] = len(idx)
        pairs.append((idx[g1], idx[g2], float(p[2])))
    names = sorted(idx, key=idx.get)
    D = np.ones((len(names), len(names)), dtype=np.float32)
    for i, j, d in pairs:
        D[i, j] = D[j, i] = d
    np.fill_diagonal(D, 0.0)
    order = np.argsort(names)
    names = [names[k] for k in order]
    D = D[np.ix_(order, order)]
    write_square(names, lambda i: D[i], a.out)
    print(f"mash distances for {len(names)} genomes written to {a.out}", file=sys.stderr)


def run_viptree(a):
    v2g = {}
    with open(a.names) as fh:
        next(fh)
        for line in fh:
            g, v = line.rstrip("\n").split("\t")
            v2g[v] = g
    with open(a.input) as fh, open(a.out, "w") as out:
        header = fh.readline().rstrip("\n").split("\t")
        missing = [h for h in header[1:] if h not in v2g]
        if missing:
            sys.exit(f"ERROR: {len(missing)} ViPTreeGen names not in {a.names}, e.g. {missing[:3]}")
        out.write("\t" + "\t".join(v2g[h] for h in header[1:]) + "\n")
        n = 0
        for line in fh:
            p = line.rstrip("\n").split("\t")
            out.write(v2g[p[0]] + "\t" + "\t".join(p[1:]) + "\n")
            n += 1
    print(f"ViPTreeGen distances for {n} genomes written to {a.out}", file=sys.stderr)


# =====================================================================
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = parser.add_subparsers(dest="command", required=True)

p = sub.add_parser("presence", help="clusters -> presence/absence matrix, phenotypes, cluster info")
p.add_argument("--clusters", required=True, help="mmseq2_cluster.tsv (representative <TAB> member)")
p.add_argument("--metadata", required=True, help="genome <whitespace> biome, with a header line")
p.add_argument("--biome1", required=True, help="coded 0 (full label or the part after the last ';')")
p.add_argument("--biome2", required=True, help="coded 1")
p.add_argument("--derep_dir", default=None, help="folder with *_derep_cluster.tsv (adds dereplicated proteins back)")
p.add_argument("--proteins_fasta", default=None, help="clustered proteins, for the annotation column")
p.add_argument("--rf_clusters", default=None, help="file with the cluster names of the multi-biome (RF) analysis")
p.add_argument("--min_proteins", type=int, default=10, help="minimum cluster size (members in the cluster table)")
p.add_argument("--out_pres", default="clusters.Rtab")
p.add_argument("--out_pheno", default="phenotypes.tsv")
p.add_argument("--out_info", default="cluster_info.tsv")
p.set_defaults(func=run_presence)

p = sub.add_parser("genomes", help="genome fastas -> one record per genome (ViPTreeGen, mash -i)")
p.add_argument("-i", "--input", nargs="+", required=True)
p.add_argument("--one_per_file", action="store_true",
               help="each file is one genome (non-phage): contigs are joined with a 100 N spacer")
p.add_argument("-o", "--out", default="genomes.fna")
p.add_argument("--names", default="genome_names.tsv")
p.set_defaults(func=run_genomes)

p = sub.add_parser("mash", help="`mash dist` output -> square matrix")
p.add_argument("-i", "--input", default="-", help="mash dist output (default: stdin)")
p.add_argument("--names", default=None, help="genome_names.tsv from 'genomes' when mash ran on genomes.fna with -i")
p.add_argument("-o", "--out", default="distances.tsv")
p.set_defaults(func=run_mash)

p = sub.add_parser("viptree", help="ViPTreeGen all.dist.matrix -> square matrix with genome IDs")
p.add_argument("-i", "--input", required=True)
p.add_argument("--names", required=True, help="genome_names.tsv from 'genomes'")
p.add_argument("-o", "--out", default="distances.tsv")
p.set_defaults(func=run_viptree)

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
args.func(args)
