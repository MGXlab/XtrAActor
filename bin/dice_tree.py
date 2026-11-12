#!/usr/bin/env python3
"""
All-v-all blastp output -> pairwise Dice distances between genomes (based on
summed bitscores of shared genes) -> neighbour-joining tree.

Expected input: BLAST tabular format (outfmt 6) where query and subject IDs
look like  <genome>_CDS_<number>, e.g. genome_A_icnght_CDS_0001
(genome = genome_A_icnght, protein = CDS_0001):

    qseqid  sseqid  pident  length  mismatch  gapopen
    qstart  qend  sstart  send  evalue  bitscore  [qcov  tcov]

The two optional trailing columns (query cover, target cover, in %) are only
needed if -qc / -tc are used.
"""
import argparse
import sys

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

EPILOG = """\
input format:
  Tab-separated BLAST output (-outfmt 6), one hit per line:

    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore [qcov tcov]

  Sequence IDs must look like <genome>_CDS_<number>, e.g. genome_A_icnght_CDS_0001
  (genome = genome_A_icnght, protein = CDS_0001). Use --sep if your IDs use a
  different separator. Columns 13 and 14 (query / target cover in %%) are optional
  and only needed for -qc / -tc. For example, DIAMOND writes them with:

    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

  (NCBI blastp has qcovs for query cover but no target-cover field.)

method:
  1. keep hits passing the evalue, identity and (optional) cover thresholds
  2. sum bitscores for every genome pair -> S[a,b]
  3. Dice distance: d(a,b) = 1 - 2*max(S[a,b], S[b,a]) / (S[a,a] + S[b,b])
     (pairs with no hits get distance 1)
  4. neighbour-joining tree from the distance matrix

output:
  out_csv   tab-separated table with columns seq1, seq2, dice
  out_tree  neighbour-joining tree in newick format

examples:
  %(prog)s -i allvall.blastp
  %(prog)s -i allvall.blastp -out_csv dist.tsv -out_tree tree.nwk -e 1e-10 -p 50
  %(prog)s -i allvall.blastp -qc 70 -tc 70
"""

parser = argparse.ArgumentParser(
    description='All-v-all blastp output -> Dice distances between genomes (from summed '
                'bitscores of shared genes) -> neighbour-joining tree.',
    epilog=EPILOG,
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--input', required=True, help='path to all-v-all blastp file')
parser.add_argument('-out_csv', '--out_csv', default='dice_distances.tsv',
                    help='tab-separated file of pairwise distances between genomes (default: dice_distances.tsv)')
parser.add_argument('-out_tree', '--out_tree', default='tree.nwk',
                    help='tree file in newick format (default: dice_tree.nwk)')
parser.add_argument('-e', '--evalue', type=float, default=1e-5,
                    help='evalue threshold to use (default: 1e-5)')
parser.add_argument('-p', '--pident', type=float, default=40,
                    help='minimum percent identity (default: 40)')
parser.add_argument('-tc', '--tcover', type=float, default=None,
                    help='minimum target cover in %% (default: none; needs column 14 of the blast file)')
parser.add_argument('-qc', '--qcover', type=float, default=None,
                    help='minimum query cover in %% (default: none; needs column 13 of the blast file)')
parser.add_argument('--sep', default='_CDS_',
                    help="text separating genome name from protein number in the IDs (default: _CDS_)")
if len(sys.argv) == 1:      # no arguments: show the full help instead of a short error
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# ---------------------------------------------------------------- read input
df = pd.read_csv(args.input, sep='\t', header=None)

base_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
             'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
if df.shape[1] < len(base_cols):
    sys.exit(f"ERROR: expected at least {len(base_cols)} tab-separated columns, "
             f"found {df.shape[1]}.")
extra_cols = ['qcov', 'tcov'][:df.shape[1] - len(base_cols)]
df = df.iloc[:, :len(base_cols) + len(extra_cols)]
df.columns = base_cols + extra_cols

# Genome name = everything before the last "_CDS_" in the sequence ID,
# e.g. genome_A_icnght_CDS_0001 -> genome_A_icnght
for id_col, genome_col in (('qseqid', 'qcontig'), ('sseqid', 'scontig')):
    ids = df[id_col].astype(str)
    bad = ids[~ids.str.contains(args.sep, regex=False)]
    if not bad.empty:
        sys.exit(f"ERROR: {len(bad)} ID(s) in {id_col} do not contain '{args.sep}', "
                 f"e.g. {bad.iloc[0]}. Use --sep to set the genome/protein separator.")
    df[genome_col] = ids.str.rsplit(args.sep, n=1).str[0]

# ------------------------------------------------------------------- filter
print("filtering blast results...")
keep = (df['evalue'] <= args.evalue) & (df['pident'] >= args.pident)
for value, col, flag in ((args.qcover, 'qcov', '-qc'), (args.tcover, 'tcov', '-tc')):
    if value is not None:
        if col not in df.columns:
            sys.exit(f"ERROR: {flag} was given but the blast file has no {col} column "
                     f"(expected as an extra column after bitscore).")
        keep &= df[col] >= value
df = df.loc[keep, ['qcontig', 'scontig', 'bitscore']].dropna()
if df.empty:
    sys.exit("ERROR: no blast hits left after filtering.")

# ------------------------------------------------------- sum bitscores
print("calculating bitscore sums...")
bs = df.groupby(['qcontig', 'scontig'])['bitscore'].sum()

# Square genome x genome matrix; pairs without any hit get a bitscore sum of 0.
names = sorted(set(df['qcontig']) | set(df['scontig']))
S = bs.unstack().reindex(index=names, columns=names).fillna(0.0).to_numpy()

# ------------------------------------------------------- Dice distance
# d(i,j) = 1 - 2 * max(S_ij, S_ji) / (S_ii + S_jj)
# originally was using a dictionary and reverse dictionary but now can do it all at once using numpy 
best = np.maximum(S, S.T)
self_bs = np.diag(S)
denom = self_bs[:, None] + self_bs[None, :]
with np.errstate(divide='ignore', invalid='ignore'):
    dice = np.where(denom > 0, 1.0 - 2.0 * best / denom, 1.0)
dice = np.clip(dice, 0.0, 1.0)   # paralogs can push the raw value slightly below 0
np.fill_diagonal(dice, 0.0)

no_self = [n for n, s in zip(names, self_bs) if s == 0]
if no_self:
    print(f"WARNING: {len(no_self)} genome(s) have no self-hits after filtering; "
          f"their distances are set to 1: {', '.join(no_self)}", file=sys.stderr)

# ------------------------------------------------------- write distances
dist_df = pd.DataFrame(dice, index=names, columns=names)
df_final = (dist_df.rename_axis(index='seq1', columns='seq2')
            .stack().rename('dice').reset_index())
df_final.to_csv(args.out_csv, sep='\t', header=True, index=False)
print(df_final)
print(f"distances written to {args.out_csv}")

# ------------------------------------------------------- build tree
if len(names) < 3:
    sys.exit("Need at least 3 genomes to build a tree; distances were still written.")

# Biopython wants a lower-triangular matrix including the diagonal
matrix = [[float(x) for x in dice[i, :i + 1]] for i in range(len(names))]
tree = DistanceTreeConstructor().nj(DistanceMatrix(names=names, matrix=matrix))
Phylo.write(tree, args.out_tree, "newick")
print(f"tree written to {args.out_tree}")