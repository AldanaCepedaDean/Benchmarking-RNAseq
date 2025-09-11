# Benchmarking RNA-seq

This repository provides a practical workflow to explore and benchmark RNA-seq data analysis.  
It includes scripts and guidelines for key steps such as:

- ✨ **UTR annotation** (if not available in the reference genome).  
- 🧬 **Quantification strategies** with selected parameters:  
  - `STAR + Salmon`  
  - `Salmon`  
  - `Kallisto`  
- 📊 **Multimapping analysis** (for STAR+Salmon and Salmon).  
- 🧪 *(Optional, second stage)* Scripts to perform **simulation experiments** to further evaluate quantification accuracy and biases.  

The goal is to provide a reproducible framework that can be adapted to different datasets,  
helping users understand the impact of mapping and quantification strategies on RNA-seq analysis.  

---
## The case of UTRs 

In many organisms, annotation files such as GFFs typically include untranslated regions (UTRs). However, in some cases these annotations are absent. Incorporating UTRs during read mapping can be particularly useful to reduce multi-mapping in highly conserved regions. To address this limitation, we explore different computational strategies to infer UTRs.
There are several existing tools to predict UTRs including GETUTR (Kim et al. 2015), UTRme (Radío et al. 2018), ExUTR (Huang and Teeling, 2017), and F3UTER (Sethi et al. 2022). However, many of these tools are limited to specific organisms or require specific input that may not be available for nonmodel organisms. peaks2utr (Haese-Hill et al. 2023) is now another option for the prediction of UTRs. 

---
## How much multi-mapping does my organism possess?

During read mapping, multi-mapping reads (i.e., reads that align to multiple genomic locations) can pose a challenge, depending on how each program handles them. Therefore, it is crucial to assess how difficult our organism of interest may be in this regard, particularly when working with short reads. To address this, we propose two different approaches to estimate the extent of multi-mapping using two programs that are widely employed in the scientific community.
 - STAR
 - SALMON

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

