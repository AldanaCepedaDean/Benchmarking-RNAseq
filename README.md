# Benchmarking RNA-seq

This repository provides a practical workflow to explore and benchmark RNA-seq data analysis.  
It includes scripts and guidelines for key steps such as:

- ✨ **UTR annotation** (if not available in the reference).  
- 🧬 **Quantification strategies** with selected parameters:  
  - `STAR + Salmon`  
  - `Salmon`  
  - `Kallisto`  
- 📊 **Multimapping analysis** (for STAR+Salmon and Salmon).  
- 🧪 *(Optional, second stage)* Scripts to perform **simulation experiments** to further evaluate quantification accuracy and biases.  

The goal is to provide a reproducible framework that can be adapted to different datasets,  
helping users understand the impact of mapping and quantification strategies on RNA-seq analysis.  

---

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
