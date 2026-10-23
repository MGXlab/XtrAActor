#!/usr/bin/env python3
"""
Build per-biome consensus sequences from a one-hot encoded dataframe.

The dataframe columns are:
  - 'protein'           : sequence identifier
  - 'ecosystem_subtype' : biome label (e.g. 'Lake', 'Oceanic')
  - '{pos}_{AA}'        : one-hot encoded amino acid columns (e.g. '0_D', '0_E', ...)

Usage:
    python consensus_from_msa.py dataframe.csv [output_prefix] [threshold]

    dataframe.csv  : path to CSV file
    output_prefix  : prefix for output files (default: "consensus")
    threshold      : min frequency to call a residue instead of 'X' (default: 0.5)

Examples:
    python consensus_from_msa.py data.csv
    python consensus_from_msa.py data.csv my_output 0.3
"""

import sys
import re
import pandas as pd
import numpy as np
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

def diagnose_position_frequencies(df, pos_to_cols):
    """Print distribution of best-AA frequencies across all positions."""
    import numpy as np

    n = len(df)
    freqs = []

    for pos, col_aa_pairs in sorted(pos_to_cols.items()):
        counts = {}
        for col, aa in col_aa_pairs:
            if col in df.columns:
                s = df[col].sum()
                if s > 0:
                    counts[aa] = s

        if not counts:
            continue

        n_present = sum(counts.values())
        if n_present / n < 0.5:
            continue

        best_freq = counts[max(counts, key=counts.get)] / n_present
        freqs.append(best_freq)

    freqs = np.array(freqs)
    print(f"\n>>> Per-position best-AA frequency distribution:")
    print(f"    N positions assessed : {len(freqs)}")
    print(f"    Mean                 : {freqs.mean():.3f}")
    print(f"    Median               : {np.median(freqs):.3f}")
    print(f"    Min                  : {freqs.min():.3f}")
    print(f"    Max                  : {freqs.max():.3f}")
    print(f"\n    Positions where best AA freq is:")
    for cutoff in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        n_above = (freqs >= cutoff).sum()
        print(f"      >= {cutoff:.1f} : {n_above:4d} / {len(freqs)} "
              f"({100*n_above/len(freqs):.1f}%)")

    print(f"\n    Positions that would be X at threshold:")
    for cutoff in [0.1, 0.2, 0.3, 0.4, 0.5]:
        n_x = (freqs < cutoff).sum()
        print(f"      threshold={cutoff:.1f} → {n_x:4d} X's "
              f"({100*n_x/len(freqs):.1f}%)")
        
def parse_ohe_columns(columns):
    """
    Parse OHE column names like '0_D', '10_A', '565_X' into
    a dict: {position (int): [list of (col_name, aa_letter)]}
    sorted by position.
    """
    pattern = re.compile(r'^(\d+)_([A-Z])$')
    pos_to_cols = defaultdict(list)

    for col in columns:
        m = pattern.match(col)
        if m:
            pos  = int(m.group(1))
            aa   = m.group(2)
            pos_to_cols[pos].append((col, aa))

    return dict(sorted(pos_to_cols.items()))


def build_consensus_from_ohe(df, pos_to_cols, threshold=0.1):
    n = len(df)
    if n == 0:
        return ""

    consensus = []
    for pos, col_aa_pairs in sorted(pos_to_cols.items()):
        counts = {}
        for col, aa in col_aa_pairs:
            if col in df.columns:
                s = df[col].sum()
                if s > 0:
                    counts[aa] = s

        if not counts:
            consensus.append("-")
            continue

        # Number of sequences that have ANY amino acid at this position
        n_present = sum(counts.values())

        # If majority have nothing encoded here, emit gap
        if n_present / n < 0.5:
            consensus.append("-")
            continue

        best_aa   = max(counts, key=counts.get)
        best_freq = counts[best_aa] / n_present  # <-- divide by present, not total

        if best_freq >= threshold:
            consensus.append(best_aa)
        else:
            consensus.append("X")

    return "".join(consensus)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    csv_file     = sys.argv[1]
    out_prefix   = sys.argv[2] if len(sys.argv) > 2 else "consensus"
    threshold    = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5

    print(f">>> Reading: {csv_file}")
    df = pd.read_csv(csv_file, index_col=0, low_memory=False, sep='\t')
    print(f">>> Shape: {df.shape}")

    # Parse OHE columns
    pos_to_cols = parse_ohe_columns(df.columns)
    n_positions = len(pos_to_cols)
    print(f">>> Found {n_positions} alignment positions")

    # Get biomes
    if "ecosystem_subtype" not in df.columns:
        print("ERROR: 'ecosystem_subtype' column not found")
        sys.exit(1)

    biomes = df["ecosystem_subtype"].dropna().unique()
    print(f">>> Biomes found: {list(biomes)}")
    print(f">>> Threshold: {threshold}\n")

    output_records = []
    stats_rows     = []

    for biome in sorted(biomes):
        subset = df[df["ecosystem_subtype"] == biome]
        n_seqs = len(subset)
        print(f">>> Building consensus for '{biome}' (n={n_seqs})...")
        diagnose_position_frequencies(subset, pos_to_cols)
        consensus = build_consensus_from_ohe(subset, pos_to_cols, threshold)

        x_count  = consensus.count("X")
        ungapped = consensus  # OHE-derived consensus has no gaps by default

        print(f"    Consensus length : {len(consensus)}")
        print(f"    Ambiguous (X)    : {x_count} ({100*x_count/max(len(consensus),1):.1f}%)\n")

        # Sanitize biome name for use in sequence ID
        biome_id = biome.replace(" ", "_").replace("/", "_")

        output_records.append(SeqRecord(
            Seq(consensus),
            id=f"consensus_{biome_id}",
            description=f"n={n_seqs} threshold={threshold} len={len(consensus)}"
        ))

        stats_rows.append({
            "biome"          : biome,
            "n_seqs"         : n_seqs,
            "n_positions"    : n_positions,
            "consensus_len"  : len(consensus),
            "ambiguous_X"    : x_count,
            "pct_ambiguous"  : round(100 * x_count / max(len(consensus), 1), 2),
        })

    # Write FASTA
    out_fasta = f"{out_prefix}_biome_consensus.fasta"
    SeqIO.write(output_records, out_fasta, "fasta")
    print(f">>> Written: {out_fasta}")

    # Write stats TSV
    out_stats = f"{out_prefix}_biome_consensus_stats.tsv"
    pd.DataFrame(stats_rows).to_csv(out_stats, sep="\t", index=False)
    print(f">>> Written: {out_stats}")

    print("\n>>> Done!")


if __name__ == "__main__":
    main()