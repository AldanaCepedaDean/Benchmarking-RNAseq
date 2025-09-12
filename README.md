# Benchmarking RNA-seq

This repository provides a practical workflow to explore and benchmark RNA-seq data analysis.  
It includes scripts and guidelines for key steps such as:

- ✨ **UTR annotation** (if not available in the reference genome).  
- 🧬 **Quantification strategies** with selected parameters:  
  - `STAR + Salmon`  
  - `Salmon`  
  - `Kallisto`  
- 📊 **Multi-mapping analysis** (for STAR+Salmon and Salmon).  
- 🧪 *(Optional, second stage)* Scripts to perform **simulation experiments** to further evaluate quantification accuracy and biases.  

The goal is to provide a reproducible framework that can be adapted to different datasets,  
helping users understand the impact of mapping and quantification strategies on RNA-seq analysis.  

---
## The case of UTRs 

In many organisms, annotation files such as GFFs typically include untranslated regions (UTRs). However, in some cases these annotations are absent. Incorporating UTRs during read mapping can be particularly useful to reduce multi-mapping in highly conserved regions. To address this limitation, we explore different computational strategies to infer UTRs.

Several tools exist for UTR prediction, including **GETUTR** (Kim et al., 2015), **UTRme** (Radío et al., 2018), **ExUTR** (Huang and Teeling, 2017), and **F3UTER** (Sethi et al., 2022). However, many of these are limited to specific organisms or depend on input formats that may not be available for non-model organisms. More recently, **peaks2utr** (Haese-Hill et al., 2023) has been introduced as an organism-independent option for UTR prediction.

### Requirements
- **GFF/GTF annotation file** (with or without existing 3′ UTR annotations)  
- Either:
  - **RNA-seq FASTQ files** (for quality assessment), or  
  - **BAM file** with aligned reads (depending on the selected tool)  

> **Note:** Some organisms, such as *Trypanosoma cruzi*, do not provide GTF files. In such cases, [`gffread`](https://github.com/gpertea/gffread) can be used to generate a GTF from the reference genome and its corresponding GFF.

### Next steps
After selecting the appropriate tool, integrate the predicted UTR coordinates into the original GFF annotation file.

---
## How much multi-mapping does my organism possess?

During read mapping, multi-mapping reads (i.e., reads that align to multiple genomic locations) can pose a challenge, depending on how each program handles them. Therefore, it is crucial to assess how difficult our organism of interest may be in this regard, particularly when working with short reads.  

To address this, we applied different strategies with two widely used tools in the community:  

### STAR
- Multi-mapping reads were identified using **secondary alignments**.  
- Unique and multi-mapped reads were separated based on alignment flags and mapping quality scores.  
- Multi-mapping percentages were calculated as the ratio of multi-mapped reads to total mapped reads, both globally and at the gene level.  

### Salmon
- Transcript-level alignments were generated with the `--writeMappings` option to produce SAM files.  
- SAM files were filtered to remove unmapped reads.  
- Transcript-level counts were extracted by summarizing mappings per transcript (using custom `awk` commands in bash).  
- Multi-mapping estimates were derived from transcript-level distributions.  

---

### Workflow Overview

```mermaid
flowchart TD
    A[Input RNA-seq reads] --> B[STAR]
    A --> C[Salmon]

    B --> B1[Secondary alignments]
    B1 --> B2[Separate unique vs. multi-mapped reads]
    B2 --> B3[Calculate multi-mapping % globally and per gene]

    C --> C1[--writeMappings → SAM output]
    C1 --> C2[Filter unmapped reads]
    C2 --> C3[Summarize transcript-level counts]
    C3 --> C4[Estimate multi-mapping % at transcript level]

## Workflow overview

```text
         ┌─────────────────────┐
         │   RNA-seq reads     │
         └─────────┬───────────┘
                   │
                   ▼
         ┌─────────────────────┐
         │  UTR annotation     │
         │ (if not available)  │
         └─────────┬───────────┘
                   │
       ┌───────────┼───────────────────┐
       │           │                   │
       ▼           ▼                   ▼
┌─────────────┐ ┌─────────────┐ ┌─────────────┐
│ STAR+Salmon │ │   Salmon    │ │  Kallisto   │
└───────┬─────┘ └──────┬──────┘ └──────┬──────┘
        │              │               │
        ▼              ▼               ▼
┌─────────────┐ ┌─────────────┐ ┌─────────────┐
│ Multimapping│ │ Multimapping│ │   (N/A)     │
│   analysis  │ │   analysis  │ │             │
└─────────────┘ └─────────────┘ └─────────────┘
                   │
                   ▼
         ┌─────────────────────┐
         │  (Optional) Sims    │
         │   & evaluation      │
         └─────────────────────┘

