#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
similar_sequences.py

All-vs-all global pairwise identity between the sequences of a FASTA file,
plus each sequence's identity to its closest paralog (best hit).

Identity = number of identical aligned positions / length of the longer
sequence * 100 (global alignment, match = 1, mismatch = 0, gaps = 0).

Author: Aldana Cepeda (2024)
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")  # save figures without needing a display
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Align import PairwiseAligner


def build_aligner():
    """Global aligner equivalent to the old pairwise2.align.globalxx."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    return aligner


def calculate_identity(aligner, seq1, seq2):
    matches = aligner.score(seq1, seq2)  # score == number of matches
    length = max(len(seq1), len(seq2))
    return matches / length * 100


def plot_identity_distribution(identities, out_png):
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    sns.histplot(identities, bins=20, kde=True, color="skyblue")
    plt.xlabel("Identity (%)")
    plt.ylabel("Frequency")
    plt.title("Distribution of pairwise sequence identity")
    plt.grid(True)

    plt.subplot(1, 2, 2)
    sns.violinplot(y=identities, inner="quartile", color="lightblue")
    plt.ylabel("Identity (%)")
    plt.title("Pairwise sequence identity")

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="All-vs-all pairwise identity and closest-paralog identity from a FASTA file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA (nucleotide or protein).")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory (default: current).")
    parser.add_argument("-p", "--prefix", default=None,
                        help="Prefix for output files (default: input file name).")
    parser.add_argument("--top", type=int, default=10, help="Number of most similar pairs to print.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.prefix or os.path.splitext(os.path.basename(args.input))[0]

    sequences = list(SeqIO.parse(args.input, "fasta"))
    if len(sequences) < 2:
        raise SystemExit("The FASTA file needs at least two sequences.")

    aligner = build_aligner()
    pairs = []
    n = len(sequences)
    for i in range(n):
        for j in range(i + 1, n):
            s1, s2 = sequences[i], sequences[j]
            pairs.append((s1.id, s2.id, calculate_identity(aligner, s1.seq, s2.seq)))

    # 1) All pairwise identities
    df = pd.DataFrame(pairs, columns=["Sequence 1", "Sequence 2", "Identity (%)"])
    df = df.sort_values("Identity (%)", ascending=False)
    pairwise_csv = os.path.join(args.outdir, f"{prefix}_pairwise_identity.csv")
    df.to_csv(pairwise_csv, index=False)

    # 2) Closest paralog per sequence
    long_df = pd.concat([
        df.rename(columns={"Sequence 1": "Sequence", "Sequence 2": "Closest paralog"}),
        df.rename(columns={"Sequence 2": "Sequence", "Sequence 1": "Closest paralog"}),
    ])
    best = (long_df.sort_values("Identity (%)", ascending=False)
                   .drop_duplicates("Sequence")
                   .sort_values("Identity (%)", ascending=False))
    best = best[["Sequence", "Closest paralog", "Identity (%)"]]
    best_csv = os.path.join(args.outdir, f"{prefix}_closest_paralog.csv")
    best.to_csv(best_csv, index=False)

    # 3) Plot
    png = os.path.join(args.outdir, f"{prefix}_identity_distribution.png")
    plot_identity_distribution(df["Identity (%)"].tolist(), png)

    print(f"Sequences: {n} | Pairs: {len(df)}")
    print(f"\nTop {args.top} most similar pairs:")
    print(df.head(args.top).to_string(index=False))
    print(f"\nOutputs:\n  {pairwise_csv}\n  {best_csv}\n  {png}")


if __name__ == "__main__":
    main()
