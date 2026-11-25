#!/usr/bin/env python3
"""
Fast version of dice_tree.py for many genomes (tested idea: ~20,000 genomes).

Same method as dice_tree.py:
  1. keep hits passing the evalue, identity and (optional) cover thresholds
  2. sum bitscores for every genome pair -> S[a,b]
  3. Dice distance: d(a,b) = 1 - 2*max(S[a,b], S[b,a]) / (S[a,a] + S[b,b])
     (pairs with no hits get distance 1)
  4. neighbour-joining tree from the distance matrix

What is different:
  - the blast file is read in pieces and reduced to bitscore sums per genome pair straight
    away, so the full hit table is never held in memory. This step can run once per blast
    chunk ('sum') and the small results are added up afterwards ('tree').
  - the tree is built with rapidNJ instead of Biopython (same algorithm, minutes instead of
    weeks for 20,000 genomes)
  - the distance table only lists each genome pair once (seq1 < seq2) and only pairs that
    share hits (dice < 1). Every pair that is not listed has distance 1.

Usage:
  # step 1, once per blast file / chunk (can run in parallel)
  dice_tree_fast.py sum  -i chunk_1_diamond.out -o chunk_1_pairsums.tsv [-e 1e-5 -p 40 -qc 70 -tc 70]

  # step 2, once
  dice_tree_fast.py tree -i *_pairsums.tsv -out_csv dice_distances.tsv -out_tree dice_tree.nwk

Blast input format (as for dice_tree.py), tab-separated, no header:
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore [qcov tcov]
Sequence IDs must look like <genome>_CDS_<number> (see --sep).
"""
import argparse
import re
import subprocess
import sys

import numpy as np
import pandas as pd

BASE_COLS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
             'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']


# =====================================================================
# step 1: blast hits -> bitscore sums per genome pair
# =====================================================================
def genome_of(ids, sep, col):
    """Genome name = everything before the last <sep> in the sequence ID."""
    codes, uniques = pd.factorize(ids)                  # split each distinct ID only once
    uniques = pd.Series(uniques.astype(str))
    bad = uniques[~uniques.str.contains(sep, regex=False)]
    if not bad.empty:
        sys.exit(f"ERROR: {len(bad)} ID(s) in {col} do not contain '{sep}', e.g. {bad.iloc[0]}. "
                 f"Use --sep to set the genome/protein separator.")
    return uniques.str.rsplit(sep, n=1).str[0].to_numpy()[codes]


def run_sum(args):
    use = ['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']
    usecols = [BASE_COLS.index(c) for c in use]
    if args.qcover is not None:
        use.append('qcov'); usecols.append(12)
    if args.tcover is not None:
        use.append('tcov'); usecols.append(13)

    totals = None
    n_hits = n_kept = 0
    for path in args.input:
        try:
            reader = pd.read_csv(path, sep='\t', header=None, usecols=usecols, names=None,
                                 chunksize=args.chunksize)
            for df in reader:
                df.columns = use                         # usecols are returned in file order = order of `use`
                n_hits += len(df)
                keep = (df['evalue'] <= args.evalue) & (df['pident'] >= args.pident)
                if args.qcover is not None:
                    keep &= df['qcov'] >= args.qcover
                if args.tcover is not None:
                    keep &= df['tcov'] >= args.tcover
                df = df.loc[keep, ['qseqid', 'sseqid', 'bitscore']].dropna()
                if df.empty:
                    continue
                n_kept += len(df)
                df['qcontig'] = genome_of(df['qseqid'], args.sep, 'qseqid')
                df['scontig'] = genome_of(df['sseqid'], args.sep, 'sseqid')
                part = df.groupby(['qcontig', 'scontig'])['bitscore'].sum()
                totals = part if totals is None else totals.add(part, fill_value=0.0)
        except pd.errors.EmptyDataError:
            print(f"WARNING: {path} is empty", file=sys.stderr)
        except ValueError as e:
            sys.exit(f"ERROR: could not read {path}: {e}\n(-qc / -tc need columns 13 / 14 of the blast file)")

    if totals is None:
        totals = pd.Series([], index=pd.MultiIndex.from_arrays([[], []], names=['qcontig', 'scontig']),
                           name='bitscore', dtype=float)
    totals.rename('bitscore_sum').reset_index().to_csv(args.output, sep='\t', index=False)
    print(f"hits read: {n_hits} | kept after filtering: {n_kept} | genome pairs: {len(totals)}")
    print(f"bitscore sums written to {args.output}")


# =====================================================================
# step 2: bitscore sums -> Dice distances -> tree
# =====================================================================
def newick_label(name):
    if re.search(r"[\s()\[\]':;,]", name):
        return "'" + name.replace("'", "''") + "'"
    return name


def run_tree(args):
    print("reading bitscore sums...")
    parts = [pd.read_csv(f, sep='\t', dtype={'qcontig': str, 'scontig': str}) for f in args.input]
    bs = pd.concat(parts, ignore_index=True)
    if bs.empty:
        sys.exit("ERROR: no blast hits left after filtering.")
    bs = bs.groupby(['qcontig', 'scontig'], sort=False)['bitscore_sum'].sum().reset_index()

    names = sorted(set(bs['qcontig']) | set(bs['scontig']))
    n = len(names)
    idx = {g: i for i, g in enumerate(names)}
    qi = bs['qcontig'].map(idx).to_numpy()
    si = bs['scontig'].map(idx).to_numpy()
    print(f"genomes: {n} | genome pairs with hits: {len(bs)}")

    # square genome x genome matrix of bitscore sums; pairs without any hit stay 0
    S = np.zeros((n, n), dtype=np.float64)
    S[qi, si] = bs['bitscore_sum'].to_numpy()
    del bs, qi, si

    # d(i,j) = 1 - 2 * max(S_ij, S_ji) / (S_ii + S_jj), calculated in blocks of rows to save memory
    print("calculating Dice distances...")
    self_bs = np.diag(S).copy()
    dice = np.empty((n, n), dtype=np.float32)
    step = 2000
    for a in range(0, n, step):
        b = min(a + step, n)
        best = np.maximum(S[a:b, :], S[:, a:b].T)
        denom = self_bs[a:b, None] + self_bs[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            d = np.where(denom > 0, 1.0 - 2.0 * best / denom, 1.0)
        dice[a:b, :] = np.clip(d, 0.0, 1.0)             # paralogs can push the raw value slightly below 0
    del S
    np.fill_diagonal(dice, 0.0)

    no_self = [g for g, s in zip(names, self_bs) if s == 0]
    if no_self:
        print(f"WARNING: {len(no_self)} genome(s) have no self-hits after filtering; "
              f"their distances are set to 1: {', '.join(no_self[:20])}"
              f"{' ...' if len(no_self) > 20 else ''}", file=sys.stderr)

    # ------------------------------------------------------- write distances
    # one row per genome pair (seq1 < seq2) with dice < 1; all other pairs have distance 1
    names_arr = np.array(names, dtype=object)
    n_pairs = 0
    with open(args.out_csv, 'w') as fh:
        fh.write("seq1\tseq2\tdice\n")
        for a in range(0, n, step):
            b = min(a + step, n)
            block = dice[a:b, :]
            r, c = np.nonzero(block < 1.0)
            r = r + a
            keep = c > r
            r, c = r[keep], c[keep]
            if len(r):
                pd.DataFrame({'seq1': names_arr[r], 'seq2': names_arr[c], 'dice': dice[r, c]}) \
                    .to_csv(fh, sep='\t', header=False, index=False, float_format='%.6f')
                n_pairs += len(r)
    print(f"distances written to {args.out_csv} ({n_pairs} genome pairs with dice < 1; "
          f"the other {n * (n - 1) // 2 - n_pairs} pairs have distance 1)")

    # ------------------------------------------------------- full square matrix (optional)
    if args.out_matrix:
        with open(args.out_matrix, 'w') as fh:
            fh.write("\t" + "\t".join(names) + "\n")
            for a in range(0, n, step):
                b = min(a + step, n)
                pd.DataFrame(dice[a:b, :], index=names[a:b]) \
                    .to_csv(fh, sep='\t', header=False, float_format='%.6f')
        print(f"square distance matrix written to {args.out_matrix}")
    if args.no_tree:
        return

    # ------------------------------------------------------- build tree
    if n < 3:
        sys.exit("Need at least 3 genomes to build a tree; distances were still written.")

    # PHYLIP matrix with neutral labels (T0, T1, ...), so genome names can contain any character
    phylip = args.out_tree + ".phylip"
    print("writing distance matrix for rapidNJ...")
    with open(phylip, 'w') as fh:
        fh.write(f"{n}\n")
        for a in range(0, n, step):
            b = min(a + step, n)
            pd.DataFrame(dice[a:b, :], index=[f"T{i}" for i in range(a, b)]) \
                .to_csv(fh, sep=' ', header=False, float_format='%.6f')
    del dice

    print("building neighbour-joining tree with rapidNJ...")
    tmp_tree = args.out_tree + ".tmp"
    cmd = [args.rapidnj, phylip, '-i', 'pd', '-o', 't', '-x', tmp_tree] + args.rapidnj_args.split()
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        sys.exit(f"ERROR: '{args.rapidnj}' not found. Install rapidNJ (conda install -c bioconda rapidnj) "
                 f"or give its path with --rapidnj. Distances were still written.")
    except subprocess.CalledProcessError as e:
        sys.exit(f"ERROR: rapidNJ failed with exit status {e.returncode}. Distances were still written.")

    # put the genome names back
    with open(tmp_tree) as fh:
        newick = fh.read()
    newick, n_sub = re.subn(r"'?T(\d+)'?(?=[:,)])", lambda m: newick_label(names[int(m.group(1))]), newick)
    if n_sub != n:
        sys.exit(f"ERROR: expected {n} tips in the rapidNJ tree, found {n_sub}.")
    with open(args.out_tree, 'w') as fh:
        fh.write(newick.strip() + "\n")
    if not args.keep_matrix:
        import os
        os.remove(phylip)
        os.remove(tmp_tree)
    print(f"tree written to {args.out_tree}")


# =====================================================================
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = parser.add_subparsers(dest='command', required=True)

p_sum = sub.add_parser('sum', help='blast file(s) -> bitscore sums per genome pair')
p_sum.add_argument('-i', '--input', required=True, nargs='+', help='all-v-all blastp file(s) or chunk(s)')
p_sum.add_argument('-o', '--output', default='pairsums.tsv', help='output table: qcontig, scontig, bitscore_sum')
p_sum.add_argument('-e', '--evalue', type=float, default=1e-5, help='evalue threshold (default: 1e-5)')
p_sum.add_argument('-p', '--pident', type=float, default=40, help='minimum percent identity (default: 40)')
p_sum.add_argument('-tc', '--tcover', type=float, default=None, help='minimum target cover in %% (column 14)')
p_sum.add_argument('-qc', '--qcover', type=float, default=None, help='minimum query cover in %% (column 13)')
p_sum.add_argument('--sep', default='_CDS_', help='text separating genome name from protein number (default: _CDS_)')
p_sum.add_argument('--chunksize', type=int, default=5_000_000, help='blast lines read at a time (default: 5000000)')
p_sum.set_defaults(func=run_sum)

p_tree = sub.add_parser('tree', help='bitscore sums -> Dice distances -> neighbour-joining tree')
p_tree.add_argument('-i', '--input', required=True, nargs='+', help="output file(s) of 'sum'")
p_tree.add_argument('-out_csv', '--out_csv', default='dice_distances.tsv',
                    help='pairwise distances, pairs with dice < 1 only (default: dice_distances.tsv)')
p_tree.add_argument('-out_tree', '--out_tree', default='dice_tree.nwk', help='tree in newick format')
p_tree.add_argument('--rapidnj', default='rapidnj', help='rapidNJ executable (default: rapidnj)')
p_tree.add_argument('--rapidnj_args', default='', help='extra rapidNJ arguments, in quotes')
p_tree.add_argument('--keep_matrix', action='store_true', help='keep the PHYLIP distance matrix (large)')
p_tree.add_argument('--out_matrix', default=None,
                    help='also write the full square distance matrix (tab-separated, genome IDs as header and '
                         'first column), e.g. for pyseer --distances')
p_tree.add_argument('--no_tree', action='store_true', help='stop after writing the distances (no rapidNJ)')
p_tree.set_defaults(func=run_tree)

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
args.func(args)
