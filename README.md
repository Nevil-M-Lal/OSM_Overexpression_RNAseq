# OSM Overexpression RNA-seq Analysis

## Project overview

This repository contains a bulk RNA-seq analysis of **Oncostatin M (OSM) overexpression in HCT116 cells**, compared with empty-vector control cells.

The workflow covers SRA data retrieval, read quality control, adapter trimming, reference-based alignment, gene-level quantification, differential gene expression analysis, gene annotation, and functional enrichment analysis.

**Dataset:** [GEO GSE350441](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE350441)

## Experimental design

The dataset contains six paired-end RNA-seq samples:

| Condition | Samples | Biological replicates |
|---|---|---:|
| OSM overexpression | OSM_1, OSM_2, OSM_3 | 3 |
| Empty-vector control | Control_1, Control_2, Control_3 | 3 |

- Organism: *Homo sapiens*
- Cell line: HCT116
- Sequencing: Paired-end, 2 × 150 bp
- Comparison: OSM overexpression vs Control

Sample accessions and condition assignments are provided in `docs/samples.tsv`.

## Workflow

1. **Data retrieval:** Download SRA files and convert them into paired-end FASTQ files using SRA Toolkit.
2. **Initial quality control:** Assess raw reads using FastQC and summarize the reports using MultiQC.
3. **Read preprocessing:** Perform quality filtering and adapter trimming using fastp.
4. **Post-trimming quality control:** Reassess processed reads using FastQC and MultiQC.
5. **Reference preparation:** Obtain the Ensembl GRCh38 release 92 gene annotation and use a prebuilt HISAT2 genome index.
6. **Alignment:** Align paired-end reads using HISAT2 and sort BAM files using samtools.
7. **Alignment quality assessment:** Evaluate BAM integrity, alignment statistics and library strandedness.
8. **Gene-level quantification:** Generate gene-level raw counts using featureCounts.
9. **Differential expression:** Perform differential expression analysis using DESeq2.
10. **Gene annotation:** Map Ensembl gene identifiers to gene symbols and Entrez identifiers using `org.Hs.eg.db`.
11. **Functional enrichment:** Perform separate enrichment analyses for upregulated and downregulated genes using Enrichr.

## Differential expression analysis

DESeq2 was used to compare OSM-overexpressing samples against empty-vector controls.

- Design formula: `~ condition`
- Reference condition: Control
- Contrast: OSM_overexpression vs Control
- Low-count filtering: Genes with a total raw count below 10 across samples were removed.
- Statistical significance: Adjusted p-value (Benjamini–Hochberg) < 0.05.

An additional subset was defined using adjusted p-value < 0.05 and absolute log2 fold change ≥ 1.

### Summary

- Genes tested: 23,500
- Significant differentially expressed genes (adjusted p-value < 0.05): 1,307
- Significant genes with absolute log2 fold change ≥ 1: 490

The fold-change threshold was applied as an additional selection criterion; it was not part of the DESeq2 model-fitting procedure.

## Functional enrichment analysis

Functional enrichment was performed separately for upregulated and downregulated gene sets using Enrichr.

The input gene sets were derived from significant DEGs satisfying:

- Upregulated: adjusted p-value < 0.05 and log2 fold change ≥ 1.
- Downregulated: adjusted p-value < 0.05 and log2 fold change ≤ −1.

A total of 432 mapped gene symbols were used: 188 upregulated and 244 downregulated.

The following libraries were examined:

- Gene Ontology Biological Process (GO BP)
- Gene Ontology Molecular Function (GO MF)
- Gene Ontology Cellular Component (GO CC)
- KEGG
- Reactome

Enrichment terms were considered significant at adjusted p-value ≤ 0.05.

### Enrichment summary

| Library | Upregulated | Downregulated |
|---|---:|---:|
| Reactome 2024 | 96 | 37 |
| KEGG 2026 | 27 | 2 |
| GO BP 2026 | 80 | 13 |
| GO CC 2026 | 26 | 0 |
| GO MF 2026 | 5 | 1 |
| **Total** | **234** | **53** |

The upregulated gene set showed enrichment involving immune responses, antigen presentation, interferon and cytokine-related processes, cell-cycle regulation, DNA replication and repair, and epithelial organization.

The downregulated gene set showed enrichment involving extracellular matrix organization, chromatin and histone regulation, cell migration, development, and epithelial proliferation regulation.

These findings represent functional overrepresentation among the selected gene sets and do not independently establish pathway activation or inhibition.

## Repository structure

```text
OSM_Overexpression_RNAseq/
├── .gitignore
├── LICENSE
├── README.md
├── docs/
│   └── samples.tsv
├── envs/
│   ├── environment.yml
│   └── explicit-spec.txt
├── logs/
│   ├── Control_1_hisat2.log
│   ├── Control_2_hisat2.log
│   ├── Control_3_hisat2.log
│   ├── OSM_1_hisat2.log
│   ├── OSM_2_hisat2.log
│   └── OSM_3_hisat2.log
├── results/
│   ├── alignment/
│   │   └── qc/
│   │       ├── Control_1.flagstat.txt
│   │       ├── Control_2.flagstat.txt
│   │       ├── Control_3.flagstat.txt
│   │       ├── OSM_1.flagstat.txt
│   │       ├── OSM_2.flagstat.txt
│   │       └── OSM_3.flagstat.txt
│   ├── counts/
│   │   ├── gene_counts.txt.gz
│   │   └── gene_counts.txt.summary
│   ├── deseq2/
│   │   ├── DESeq2_OSM_vs_Control.csv
│   │   ├── DESeq2_dataset.rds
│   │   ├── DESeq2_summary.txt
│   │   ├── normalized_counts.csv
│   │   ├── plots/
│   │   │   ├── DEG_heatmap.png
│   │   │   ├── MA_plot.png
│   │   │   ├── PCA.png
│   │   │   ├── sample_correlation_heatmap.png
│   │   │   ├── sessionInfo.txt
│   │   │   └── volcano_plot.png
│   │   ├── sessionInfo.txt
│   │   ├── significant_DEGs_padj05.csv
│   │   └── significant_DEGs_padj05_lfc1.csv
│   ├── enrichr/
│   │   ├── DOWN/
│   │   │   ├── DOWN_Enrichment_Results.xlsx
│   │   │   ├── DOWN_GO_BP_2026_bar.png
│   │   │   ├── DOWN_GO_CC_2026_bar.png
│   │   │   ├── DOWN_GO_MF_2026_bar.png
│   │   │   ├── DOWN_KEGG_2026_bar.png
│   │   │   └── DOWN_Reactome_2024_bar.png
│   │   └── UP/
│   │       ├── UP_Enrichment_Results.xlsx
│   │       ├── UP_GO_BP_2026_bar.png
│   │       ├── UP_GO_CC_2026_bar.png
│   │       ├── UP_GO_MF_2026_bar.png
│   │       ├── UP_KEGG_2026_bar.png
│   │       └── UP_Reactome_2024_bar.png
│   ├── fastp/
│   │   ├── SRR41141592_fastp.html
│   │   ├── SRR41141593_fastp.html
│   │   ├── SRR41141594_fastp.html
│   │   ├── SRR41141595_fastp.html
│   │   ├── SRR41141596_fastp.html
│   │   └── SRR41141597_fastp.html
│   └── qc/
│       ├── multiqc_raw/
│       │   └── multiqc_report.html
│       └── multiqc_trimmed/
│           └── multiqc_report.html
└── scripts/
    ├── deseq2_analysis.R
    ├── deseq2_visualization.R
    ├── gtf_to_bed12.py
    └── preprocessing.sh
```

## Scripts

| Script | Purpose |
|---|---|
| `preprocessing.sh` | SRA retrieval, FASTQ conversion, read QC, trimming, alignment, alignment QC and gene quantification |
| `gtf_to_bed12.py` | Converts exon records in the GTF annotation into transcript-level BED12 records for RSeQC |
| `deseq2_analysis.R` | Differential expression analysis, gene annotation and result generation |
| `deseq2_visualization.R` | PCA, sample correlation heatmap, MA plot, volcano plot and DEG heatmap |

## Reproducibility

The Conda environment specification is provided in `envs/environment.yml`. The explicit environment specification is available in `envs/explicit-spec.txt`.

The analysis was performed using the following principal software versions:

- HISAT2 2.2.3
- samtools 1.24
- fastp 1.4.0
- FastQC 0.13.0
- MultiQC 1.35
- Subread/featureCounts 2.1.1
- R 4.4.3
- DESeq2 1.46.0

The R analysis scripts save session information alongside their outputs.

## Limitations

- The analysis includes three biological replicates per condition.
- A prebuilt transcript-aware HISAT2 GRCh38 index was used instead of building a new index locally. Exact nucleotide sequence identity between the prebuilt index and the downloaded Ensembl release 92 FASTA was not independently verified.
- Library strandedness was assessed using representative samples (OSM_1 and Control_1). The observed results supported the use of unstranded counting (`-s 0`).
- Of the 490 significant DEGs selected using the additional fold-change threshold, 432 were mapped to gene symbols. The remaining 58 lacked corresponding SYMBOL and ENTREZID annotations in the installed `org.Hs.eg.db` database and were excluded from functional enrichment analysis.
- Functional enrichment results represent overrepresentation of annotated genes and should not be interpreted as direct experimental confirmation of pathway activity.
- No independent experimental validation of the differential expression or enrichment findings was performed as part of this computational workflow.

## Acknowledgements and declaration

This is an individual computational biology project carried out using publicly available RNA-seq data and the computational resources available on my personal laptop. The workflow was completed within the practical limitations of the available hardware and computing environment.

AI tools were used as supplementary assistance during code development, debugging, and documentation writing, including improving the grammar, clarity, and presentation of written content. The analytical decisions, selection of parameters, evaluation of results, and biological interpretation were reviewed and undertaken by me.

Some methodological and computational limitations remain due to the available resources and scope of the project. These limitations are acknowledged in this repository. This project is intended as an independent learning and portfolio project demonstrating the application of a bulk RNA-seq analysis workflow, rather than a production-grade or clinically validated pipeline.

## License

This repository is distributed under the MIT License. The license applies to the original code and documentation in this repository. External datasets, reference annotations and software retain their respective terms and licenses; users should consult the original sources for those conditions.
