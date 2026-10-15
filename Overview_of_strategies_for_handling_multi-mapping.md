### Comparative overview of strategies for handling multi-mapping reads in RNA-seq

Summary of commonly used alignment, pseudoalignment and quantification approaches for handling multi-mapping reads in bulk and single-cell RNA-seq, including their main advantages, limitations and recommended applications (Supplementary Table 4 of the manuscript). The table is not intended to be exhaustive; tools were selected as representative examples to serve as a practical guide.

<details>
<summary><strong>Bulk RNA-seq approaches</strong></summary>

| Tool / Approach | Type | Pros | Cons | Recommended scenarios |
|---|---|---|---|---|
| **Bowtie2 + featureCounts** | Alignment-based | Mature, widely used; generates alignments; flexible parameterization | Retains many multimappers → inflated counts in repetitive regions; poor paralog resolution | Need for BAM alignments; exploratory analyses |
| **STAR + featureCounts** | Alignment-based | Fast splice-aware aligner; widely supported; extended use; flexible parameterization | Handling of multimappers depends on settings; sensitive to incomplete annotation | Downstream analyses requiring BAMs or junctions |
| **STAR + Salmon** | Hybrid alignment + quantifier | Keeps alignments while benefiting from probabilistic quantification; good compromise between speed/accuracy | Sensitive to annotation quality; heavier computational burden than pure Salmon | When alignments are required (QC, variant calling) but accurate quantification is needed |
| **Salmon** | Alignment-free (quasi-mapping) | High accuracy; probabilistic assignment of multimappers; fast and memory-efficient | No full BAM alignments; cannot resolve perfect duplicates; sensitive to incomplete transcript models | General-purpose quantification; repetitive genomes with extended annotations |
| **Kallisto** | Alignment-free (pseudoalignment) | Extremely fast; good accuracy at transcript level; low resource requirements | Less flexible than Salmon; fewer options for inspecting ambiguous read assignments | Transcriptome-wide quantification; exploratory analyses at scale |
| **MMquant** | Alignment-based read reassignment | Redistributes multimappers across compatible features → maximizes data usage | Risk of overestimation in pseudogenes/poorly resolved paralogs; less adopted, partly due to its more complex installation and configuration | Multigene family studies where coverage retention is critical |
| **STARconsensus** | Alignment-based (haplotype-aware) | Allows haplotype-specific reference transformation; supports consensus genomes | Complex setup; not widely adopted; limited tools downstream | Studies with divergent haplotypes; allele-specific expression |
| **Graph-based pangenomes (vg, Giraffe, minigraph)** | Graph aligners | Reduces reference bias; integrates haplotypes/structural variants; handles diversity | Still experimental in RNA-seq; complex to implement; computationally intensive | Promising for highly polymorphic/repetitive species; pangenome projects |

</details>

---

<details>
<summary><strong>Single-cell RNA-seq approaches</strong></summary>

| Tool / Approach | Type | Pros | Cons | Recommended scenarios |
|---|---|---|---|---|
| **Cell Ranger (10x)** | Single-cell pipeline (STAR-based) | Industry standard; robust UMI handling; reproducible outputs; discards multimappers conservatively | Underestimates expression of repetitive gene families; STAR dependency; less customizable | Routine 10x scRNA-seq; industry-grade analyses |
| **STARsolo** | Single-cell aligner (STAR-based) | Transparent STAR-based pipeline; customizable multimapper handling | Same multimapping limitations as STAR; heavier resource use | Alternative to Cell Ranger; when flexibility is needed in scRNA-seq |
| **Salmon-Alevin / Alevin-fry** | Single-cell quantifier (alignment-free) | Efficient UMI-aware quantification; probabilistic handling of multimappers; lower memory footprint | Less extended use | scRNA-seq with focus on quantification accuracy; repetitive transcriptomes |
| **Kallisto + bustools** | Single-cell quantifier (pseudoalignment) | Extremely fast; modular and memory-efficient; widely benchmarked | Requires careful pre-processing; less user-friendly | scRNA-seq at scale where speed and modularity are critical |

</details>
