# Benchmarking RNA-seq

This repository accompanies the manuscript and provides a practical workflow to explore and benchmark RNA-seq data analysis, with a focus on organisms with highly repetitive genomes and large multigene families.

It includes scripts and guidelines to:

1. 📊 **Understand multi-mapping rates in your dataset** (STAR, Salmon and Bowtie2).
2. 🔍 **Understand the presence of highly similar genes and sequences** in your organism.
3. 🧬 **Choose a quantification strategy**, with parameters selected for repetitive genomes:
   - `STAR + Salmon`
   - `Salmon`
   - `Kallisto`
4. ✨ **Annotate UTRs**, if they are not available in the reference genome.
5. 🧪 *(Optional)* **Run simulation experiments** to further evaluate quantification accuracy and biases.

The goal is to provide a reproducible framework that can be adapted to different datasets, helping users understand the impact of repetitiveness, mapping and quantification strategies on RNA-seq analysis.

---

## 1. Understanding multi-mapping rates in your dataset

Multi-mapping reads (i.e., reads that align equally well to more than one location) are a major source of uncertainty in RNA-seq quantification. The amount of multi-mapping in a dataset is not a fixed property of the organism: it results from the combination of **genome repetitiveness**, **read length** and **annotation quality** (e.g., whether UTRs are annotated). Estimating it before choosing a pipeline helps anticipate how much each tool's strategy for handling ambiguous reads will affect your results.

We provide approaches to estimate global and per-gene/per-transcript multi-mapping rates with **STAR**, **Salmon** and **Bowtie2**.

➡️ **How to run it:** [docs/multimapping.md](docs/multimapping.md)

---

## 2. Understanding the presence of highly similar genes and sequences

<!-- TODO: brief description of the approach (e.g., pairwise identity to the closest paralog) and link to the corresponding script/doc. -->

Genes with very high sequence identity to other genes in the genome (e.g., members of multigene families) are especially hard to quantify, since short reads often cannot be unambiguously assigned to a single copy. Knowing which genes of interest fall in this category helps interpret their expression estimates with the appropriate caution.

➡️ **How to run it:** [docs/similar_sequences.md](docs/similar_sequences.md)

---

## 3. Quantification strategies

There is an abundance of strategies for handling multi-mapping reads, both for bulk and single-cell RNA-seq (see our [comparative overview of tools](docs/strategies_comparison.md)).

As discussed in our manuscript, we found that the best-performing aligner–quantifier combinations, especially for highly repetitive genomes and for genes/transcripts with high sequence similarity, are those that assign multi-mapping reads probabilistically: **STAR + Salmon**, **Salmon** and **Kallisto**. In contrast, count-based approaches such as Bowtie2 + featureCounts strongly overestimated read counts for multigene families. Quantification error increased sharply for genes sharing more than ~98–99% identity with their closest paralog, and Salmon showed consistently lower error than Kallisto in that range.

<!-- TODO: replace/adjust this paragraph with the wording from the final abstract/results of the manuscript. -->

All reads were mapped or quantified against the *T. cruzi* YC6 reference genome and transcriptome (Wang et al., 2021). Parameters were adjusted to account for the complex genome architecture of *T. cruzi*.

### 3.1 STAR + Salmon (alignment-based mode)

This strategy combines **STAR** for alignment and **Salmon** for transcript quantification. It is a convenient choice when you need to inspect the actual alignments or require them for other downstream analyses (e.g., QC, coverage, variant calling), while still benefiting from Salmon's probabilistic quantification.

Paired-end reads are first aligned to the reference genome with STAR, which outputs coordinate-sorted genomic BAM files. Using the gene annotation, STAR also projects these genomic alignments onto transcript coordinates, producing a separate transcriptome BAM file. This transcriptome-projected BAM is then used as input for Salmon in alignment-based mode.

**Key parameters (STAR):**
- Allow multiple mappings per read: `--outFilterMultimapNmax 50`
- Adapt splice junction detection:
  `--outSJfilterReads Unique`
  `--outSJfilterOverhangMin 150 150 150 150`
  `--outFilterType BySJout`
- Mismatch filtering: `--outFilterMismatchNoverReadLmax 0.01`
- End-to-end alignment: `--alignEndsType EndToEnd`
- Output sorted BAM files: `--outSAMtype BAM SortedByCoordinate`
- Output transcriptome alignments: `--quantMode TranscriptomeSAM`

> [!IMPORTANT]
> **Organism-specific parameterization**
>
> These parameters were optimized for *Trypanosoma cruzi* and must be adjusted for other organisms.
>
> - *T. cruzi* has a highly repetitive genome, particularly due to large multigene families → `--outFilterMultimapNmax 50` was set to allow a high number of valid multi-mapping alignments per read.
> - *T. cruzi* does not contain canonical introns → the splice-related parameters were configured to prevent artificial splice junction detection and avoid spliced alignments.
>
> For other organisms:
> - **Complex exon–intron structures:** STAR's defaults are optimized for mammalian genomes (mainly human and mouse), so they will probably work better than the parameters above unless you plan a careful, organism-specific optimization.
> - **Lower repetitiveness:** a smaller `--outFilterMultimapNmax` value may be more appropriate.

**Salmon (alignment-based mode):**
- Input: transcriptome-projected BAM files from STAR
- Run with default parameters

### 3.2 Salmon (quasi-mapping mode)

This strategy uses **Salmon** directly with its quasi-mapping algorithm, without prior genome alignment.

**Steps:**
1. Index the genome and GTF annotation.
2. Run Salmon with paired-end reads.

**Key parameters:**
- Automatic library type detection: `-l A`
- Improved mapping validation: `--validateMappings`
- Default parameters for quantification

### 3.3 Kallisto (pseudo-alignment mode)

This strategy uses **Kallisto** for fast quantification through pseudo-alignment against the *T. cruzi* transcriptome FASTA.

**Steps:**
1. Build the Kallisto index from the transcript FASTA.
2. Run Kallisto with paired-end reads.

**Key parameters:**
- Default pseudo-alignment settings
- Input: transcriptome FASTA

> [!NOTE]
> All commands should be adapted to your file paths and computational environment. The provided parameters are optimized for the *T. cruzi* genome, which is highly repetitive and requires careful handling of multi-mapping reads.

<!-- Optional: case study with Florencia's data showing how each strategy performs on a multigene family -> docs/case_study.md -->

---

## 4. UTR annotation

In many organisms, annotation files such as GFFs typically include untranslated regions (UTRs). However, in some cases these annotations are absent. Incorporating UTRs during read mapping can be particularly useful to reduce ambiguity and, hence, multi-mapping in highly conserved regions. To address this limitation, we explored different computational strategies to infer UTRs.

Several tools exist for UTR prediction, including **GETUTR** (Kim et al., 2015), [`UTRme`](https://github.com/sradiouy/UTRme) (Radío et al., 2018), **ExUTR** (Huang and Teeling, 2017) and **F3UTER** (Sethi et al., 2022). However, many of these are limited to specific organisms or depend on input formats that may not be available for non-model organisms. More recently, [`peaks2utr`](https://github.com/haessar/peaks2utr) (Haese-Hill et al., 2023) has been introduced as an organism-independent option for UTR prediction.

#### General requirements
- GFF/GTF annotation file (with or without existing 3′ UTR annotations)
- Either RNA-seq FASTQ files or a BAM file with aligned reads, depending on the selected tool

<sub>Some organisms, such as *T. cruzi*, do not provide GTF files. In such cases, [`gffread`](https://github.com/gpertea/gffread) can be used to generate a GTF from the reference genome and its corresponding GFF.</sub>

#### Next steps
After selecting the appropriate tool, integrate the predicted UTR coordinates into the original GFF annotation file.
<!-- TODO: specify whether this was done with a tool or with custom scripts (and link them if so). -->

---

## 5. Simulations

To perform RNA-seq simulations, we highly recommend the [**Polyester**](https://github.com/alyssafrazee/polyester) R package. It only requires a FASTA file with the transcript/CDS sequences and a table with the expected number of reads per gene, and lets you simulate different sequencing depths by scaling those counts.

➡️ **Example script:** [scripts/Simulation_example.R](scripts/Simulation_example.R)
