# Benchmarking-RNAseq
RNA-seq benchmarking toolkit including UTR annotation, quantification (STAR+Salmon, Salmon, Kallisto), multimapping evaluation, and optional simulation scripts.


This repository provides a practical workflow to explore and benchmark RNA-seq data analysis. It includes scripts and guidelines for key steps such as:

UTR annotation (if not available in the reference).

Running three quantification strategies with the selected parameters:

STAR + Salmon

Salmon

Kallisto

Multimapping analysis for the first two strategies.

(Optional, second stage) Scripts to perform simulation experiments to further evaluate quantification accuracy and biases.

The goal is to provide a reproducible framework that can be adapted to different datasets, helping users understand the impact of mapping and quantification strategies on RNA-seq analysis.
