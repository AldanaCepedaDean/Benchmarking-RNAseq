# Identifying highly similar genes and sequences

This step estimates how similar each gene of interest is to the other genes in a set (e.g., all members of a multigene family). Genes with very high identity to a paralog are the ones most likely to suffer from ambiguous read assignment during RNA-seq quantification, so their expression estimates should be interpreted with caution.

## Approach

For every pair of sequences in a FASTA file, the script computes a **global pairwise alignment** and reports the percent identity:

```
Identity (%) = identical aligned positions / length of the longer sequence × 100
```

The alignment uses a simple scoring scheme (match = 1, mismatch = 0, gaps = 0), so the alignment score equals the number of identical positions. Dividing by the length of the longer sequence penalizes length differences: a short fragment that is identical to part of a longer gene will not reach 100%.

From the all-vs-all table, the script also reports, for each sequence, its **closest paralog** (the other sequence with the highest identity). This per-gene value is the most useful one to cross with expression or quantification-error results.

## Requirements

- Python ≥ 3.8
- Biopython, pandas, matplotlib, seaborn

```bash
pip install biopython pandas matplotlib seaborn
```

## Input

A FASTA file with the sequences to compare (nucleotide or protein, but not mixed). For example, the CDS of all members of a gene family:

```
>gene_A
ATGGCT...
>gene_B
ATGGCC...
```

Use the same type of sequence (e.g., CDS vs. full transcripts) that reads will be mapped to, so identity reflects what the quantifier actually sees.

## Usage

```bash
python scripts/similar_sequences.py \
    -i data/DGF1_YC6.fasta \
    -o results/similar_sequences \
    -p DGF1_YC6
```

| Argument | Description | Default |
|---|---|---|
| `-i`, `--input` | Input FASTA file | required |
| `-o`, `--outdir` | Output directory | current directory |
| `-p`, `--prefix` | Prefix for output files | input file name |
| `--top` | Number of most similar pairs printed to screen | 10 |

## Outputs

| File | Content |
|---|---|
| `<prefix>_pairwise_identity.csv` | All pairs: `Sequence 1`, `Sequence 2`, `Identity (%)`, sorted from most to least similar |
| `<prefix>_closest_paralog.csv` | One row per sequence: `Sequence`, `Closest paralog`, `Identity (%)` |
| `<prefix>_identity_distribution.png` | Histogram (with KDE) and violin plot of all pairwise identities |

The top pairs and the output paths are also printed to the terminal.

## Interpreting the results

- **Closest-paralog identity** is the key metric: a gene whose closest paralog is, for example, >95% identical over its full length is likely to share many reads with it and to have less reliable expression estimates.
- The **distribution plot** gives an overview of the family: a long tail at high identity indicates groups of near-identical copies.
- These values can be combined with the quantification benchmark (simulated reads) to relate identity to the quantification error of each method.

## Notes and limitations

- **Runtime grows quadratically** with the number of sequences (n·(n−1)/2 alignments). It is fine for hundreds of sequences; for thousands, consider a faster pre-filter (e.g., BLAST/DIAMOND all-vs-all or MMseqs2) and align only the candidate pairs.
- Identity is computed over the full length of the sequences; it does not detect **local** stretches of high identity (e.g., shared conserved domains or repeated motifs) that can still cause ambiguous mapping of short reads.
- Nucleotide identity is the relevant measure for read mapping; protein identity can be used to describe the family but does not directly reflect mapping ambiguity.
