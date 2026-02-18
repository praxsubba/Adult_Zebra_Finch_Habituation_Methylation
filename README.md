# DNA Methylation of Cell Adhesion Genes is Associated with Habituation of the Song Response in Zebra Finches
## Supplementary Code Repository — RRBS Analysis Workflow

This repository contains the full bioinformatics workflow used in the manuscript, including read alignment, CpG methylation extraction, differential methylation analysis, and genomic annotation. All steps are described in sequential order.

**Reference Genome:** `GCF_003957565.2_bTaeGut1.4.pri_genomic.fna`  
Available at: https://www.ncbi.nlm.nih.gov/assembly/GCF_003957565.2/

---

## Table of Contents

1. [BSMAP: RRBS Read Alignment](#1-bsmap-rrbs-read-alignment)
2. [Minimap2: cDNA-to-Transcriptome Mapping](#2-minimap2-cdna-to-transcriptome-mapping)
3. [Gene Interval Definition (R)](#3-gene-interval-definition-r)
4. [samtools: Extracting Target Alignments](#4-samtools-extracting-target-alignments)
5. [BSMAPz: CpG Methylation Extraction](#5-bsmapz-cpg-methylation-extraction)
6. [methylKit: Differential Methylation Analysis](#6-methylkit-differential-methylation-analysis)
7. [genomation: Annotation](#7-genomation-annotation)
8. [Session Info](#8-session-info)

---

## 1. BSMAP: RRBS Read Alignment

Map RRBS reads to the zebra finch reference genome using [BSMAP](https://github.com/genome-vendor/bsmap).

**Input:**
- `<sample>.fastq` — Raw RRBS reads
- `GCF_003957565.2_bTaeGut1.4.pri_genomic.fna` — Reference genome

**Command:**

```bash
bsmap \
  -a <input_fastq> \
  -d GCF_003957565.2_bTaeGut1.4.pri_genomic.fna \
  -o <output_bam> \
  -D C-CGG \
  -D T-CGA \
  -w 100 \
  -v 0.08 \
  -r 0 \
  -p 4 \
  -n 0 \
  -s 12 \
  -S 0 \
  -f 5 \
  -q 0 \
  -u \
  -V 2
```

**Parameter Notes:**
| Flag | Value | Description |
|------|-------|-------------|
| `-D` | `C-CGG`, `T-CGA` | Restriction enzyme cut sites (RRBS mode) |
| `-w` | `100` | Maximum number of equal best hits |
| `-v` | `0.08` | Maximum mismatch rate |
| `-r` | `0` | Report all best hits |
| `-p` | `4` | Number of threads |
| `-s` | `12` | Seed size |
| `-f` | `5` | Filter reads with >5 Ns |
| `-u` | — | Report unmapped reads |

---

## 2. Minimap2: cDNA-to-Transcriptome Mapping

Map cDNA clone sequences from Dong et al. (2009) to the zebra finch transcriptome to obtain clone ID–to–transcript ID mappings.

**Reference Transcriptome:** `GCF_003957565.2_bTaeGut1.4.pri_rna.fna`  
Available at: https://www.ncbi.nlm.nih.gov/assembly/GCF_003957565.2/

### 2.1 Alignment

```bash
minimap2 -I 13G \
  GCF_003957565.2_bTaeGut1.4.pri_rna.fna \
  sb_array_seq.FASTA \
  > minimap2_RNA_output_clone_transcript.paf
```

### 2.2 Post-Processing in R

Parse minimap2 output and join with gene name annotations extracted from the reference FASTA.

```r
library(dplyr)

# Load minimap2 output
minimap_output <- read.csv(
  file = "minimap2_RNA_output_clone_transcript.csv",
  sep  = "\t"
)
colnames(minimap_output) <- c("Clone_ID", "Transcript_ID")

# Load transcript-to-gene name table
# (extracted from GCF_003957565.2_bTaeGut1.4.pri_rna.fna headers using awk)
transcript_gene <- read.csv(
  file = "transcript_gene_names_zebra_finch.csv",
  sep  = "\t"
)
colnames(transcript_gene) <- c("Transcript_ID", "gene")

# Join to get Clone_ID -> Transcript_ID -> Gene mapping
clone_transcript_gene <- inner_join(transcript_gene, minimap_output, by = "Transcript_ID")
head(clone_transcript_gene)
```

---

## 3. Gene Interval Definition (R)

Define gene body intervals extended 2 kb upstream to capture potential regulatory regions. Filter to significant differentially expressed genes (Habituated vs. Novel; AASA analysis).

```r
library(rtracklayer)
library(tidyverse)

# --- 3.1 Load significant DE genes ---
hab_vs_nov_sig_genes <- read.csv("AASA_Manuscript_plot.csv")

hab_vs_nov_sig_genes_final <- hab_vs_nov_sig_genes %>%
  dplyr::select(gene, clone_id, AASA_fold_clone, AASA_fdr_clone, log2FoldChange, fdr) %>%
  filter(AASA_fdr_clone < 0.05)

# --- 3.2 Parse gene coordinates from GFF ---
# Select only "gene"-type features (one row per gene, not per transcript)
my_tags    <- c("Name", "Dbxref", "gene")
my_columns <- c("seqid", "start", "end", "strand", "type")
my_filter  <- list(type = "gene")

dat <- readGFF(
  "GCF_003957565.2_bTaeGut1.4.pri_genomic.gff",
  tags    = my_tags,
  columns = my_columns,
  filter  = my_filter
)

# --- 3.3 Extend intervals 2 kb upstream of transcription start site ---
# For genes on the (+) strand: interval_start = start - 2000
# For genes on the (-) strand: interval_stop  = end   + 2000
dat_intervals <- as.data.frame(dat) %>%
  mutate(
    interval_start = ifelse(strand == "+", start - 2000, start),
    interval_start = ifelse(interval_start < 1, 1, interval_start),  # clamp to chromosome start
    interval_stop  = ifelse(strand == "-", end + 2000, end)
  )

saveRDS(dat_intervals, "bTaeGut_v2p_gene_intervals_plus_2kb.RDS")

# --- 3.4 Filter intervals to significant DE genes only ---
sig_gene_intervals <- dplyr::semi_join(dat_intervals, hab_vs_nov_sig_genes_final, by = "gene")

saveRDS(sig_gene_intervals, "sig_Habituated_vs_Novel_gene_intervals.RDS")

write.table(
  sig_gene_intervals,
  file      = "sig_Habituated_vs_Novel_gene_intervals.txt",
  sep       = "\t",
  row.names = FALSE
)
# Note: The output .txt file (interval_start, interval_stop columns) is used
# as a BED file in steps 4 and 5.
```

---

## 4. samtools: Extracting Target Alignments

Subset BSMAP BAM files to reads overlapping the significant gene intervals defined above.

```bash
samtools view -b \
  -L sig_Habituated_vs_Novel_gene_intervals.bed \
  <sample>_output.bam \
  > <sample>_interval.bam
```

---

## 5. BSMAPz: CpG Methylation Extraction

Extract per-CpG methylation ratios from the interval-subset BAM files using [BSMAPz](https://github.com/dohlee/bsmapz).

```bash
python methratio.py \
  -o <sample>_methratio.txt \
  -d GCF_003957565.2_bTaeGut1.4.pri_genomic.fna \
  -z \
  -x CG \
  <sample>_interval.bam
```

---

## 6. methylKit: Differential Methylation Analysis

Differential methylation analysis comparing **Habituated (FA)** vs. **Silence (SI)** conditions using [methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html).

### 6.1 Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylKit")
library(methylKit)
```

### 6.2 Load Methylation Data

```r
# Input files: per-sample methratio output from BSMAPz (Step 5)
# Conditions: FA = Habituated (treatment = 1), SI = Silence (treatment = 0)

file.list <- list(
  "G1FA160_methratio.txt", "G1FA211_methratio.txt", "G1FA237_methratio.txt",
  "G2FA149_methratio.txt", "G2FA209_methratio.txt", "G2FA222_methratio.txt",
  "G1SI152_methratio.txt", "G1SI199_methratio.txt", "G1SI246_methratio.txt",
  "G2SI146_methratio.txt", "G2SI188_methratio.txt", "G2SI218_methratio.txt"
)

myobj <- methRead(
  file.list,
  pipeline = list(
    fraction    = TRUE,
    chr.col     = 1,
    start.col   = 2,
    end.col     = 2,
    coverage.col = 6,
    strand.col  = 3,
    freqC.col   = 5
  ),
  sample.id = list(
    "G1FA160", "G1FA211", "G1FA237",
    "G2FA149", "G2FA209", "G2FA222",
    "G1SI152", "G1SI199", "G1SI246",
    "G2SI146", "G2SI188", "G2SI218"
  ),
  assembly  = "taeGut1",
  treatment = c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
)
```

### 6.3 Quality Control

```r
# Methylation statistics per sample
for (i in 1:12) {
  getMethylationStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
}
```

### 6.4 Coverage Filtering

```r
# Retain CpGs with a minimum read depth of 10
filtered.myobj <- filterByCoverage(
  myobj,
  lo.count = 10,
  lo.perc  = NULL,
  hi.count = NULL,
  hi.perc  = NULL
)

# QC plots post-filtering
for (i in 1:12) {
  getCoverageStats(filtered.myobj[[i]], plot = TRUE, both.strands = FALSE)
  getMethylationStats(filtered.myobj[[i]], plot = TRUE, both.strands = FALSE)
}
```

### 6.5 Unite Samples

```r
# Unite CpG sites across all samples (strand-specific)
meth <- unite(filtered.myobj, destrand = FALSE)

# Unite with destranding (merges CpGs on opposite strands; increases coverage)
meth_destrand <- unite(filtered.myobj, destrand = TRUE)
```

### 6.6 Sample Clustering and PCA

```r
# Hierarchical clustering
clusterSamples(meth,          dist = "correlation", method = "ward", plot = TRUE, sd.threshold = 0.90)
clusterSamples(meth_destrand, dist = "correlation", method = "ward", plot = TRUE, sd.threshold = 0.90)

# PCA
PCASamples(meth,          adj.lim = c(0.4, 1), sd.threshold = 0.90)
PCASamples(meth_destrand, adj.lim = c(0.0, 1), sd.threshold = 0.50)
```

### 6.7 Differential Methylation

```r
# Calculate differential methylation (strand-specific and destranded)
myDiff   <- calculateDiffMeth(meth)
myDiff_d <- calculateDiffMeth(meth_destrand)

# Export full destranded results
getData(myDiff_d) %>%
  dplyr::arrange(qvalue) %>%
  write.csv(file = "differential_methylation_all_CpGs.csv")

# Filter: q-value < 0.01, methylation difference > 25%
myDiff10p_destranded <- getMethylDiff(myDiff_d, qvalue = 0.01, difference = 25)

getData(myDiff10p_destranded) %>%
  dplyr::arrange(qvalue) %>%
  write.csv("differential_methylation_sig_CpGs_unannotated.csv")
```

---

## 7. genomation: Annotation

Annotate differentially methylated CpGs with promoter, exon, and intron features using [genomation](https://bioconductor.org/packages/release/bioc/html/genomation.html).

```r
library(genomation)

# --- 7.1 Load gene model BED file ---
# BED file derived from UCSC genePred table for taeGut1
gene.obj <- readTranscriptFeatures(
  "genePred.bed.txt",
  up.flank      = 2000,
  down.flank    = 0,
  remove.unusual = FALSE
)

# --- 7.2 Annotate significant DMCs ---
diffAnn          <- annotateWithGeneParts(as(myDiff10p_destranded, "GRanges"), gene.obj)
diffAnn_features <- annotateWithFeatures(as(myDiff10p_destranded, "GRanges"), gene.obj, intersect.chr = TRUE)

plotTargetAnnotation(
  diffAnn_features,
  precedence = TRUE,
  main       = "Habituated vs. Silence — Significant DMC Annotation"
)
plotTargetAnnotation(diffAnn, precedence = TRUE)

head(getAssociationWithTSS(diffAnn))
getTargetAnnotationStats(diffAnn, percentage = TRUE, precedence = TRUE)

# --- 7.3 Annotate all CpGs (background) ---
diffAnn_ALL          <- annotateWithGeneParts(as(myDiff_d, "GRanges"), gene.obj, intersect.chr = TRUE)
diffAnn_features_ALL <- annotateWithFeatures(as(myDiff_d, "GRanges"), gene.obj, intersect.chr = TRUE)

plotTargetAnnotation(
  diffAnn_features_ALL,
  precedence = TRUE,
  main       = "Habituated vs. Silence — All CpG Annotation"
)

head(getAssociationWithTSS(diffAnn_ALL))
getTargetAnnotationStats(diffAnn_ALL, percentage = TRUE, precedence = TRUE)

# --- 7.4 Promoter region counts ---
promoters <- regionCounts(filtered.myobj, gene.obj$promoters)
head(promoters[])[1]
```

---

## 8. Session Info

```r
sessionInfo()
```

---

## Dependencies

| Tool | Version Used | Source |
|------|-------------|--------|
| BSMAP | 2.9 | https://github.com/genome-vendor/bsmap |
| minimap2 | — | https://github.com/lh3/minimap2 |
| BSMAPz (methratio.py) | — | https://github.com/dohlee/bsmapz |
| samtools | — | http://www.htslib.org/ |
| R / methylKit | Bioconductor | https://bioconductor.org/packages/methylKit |
| R / genomation | Bioconductor | https://bioconductor.org/packages/genomation |
| R / rtracklayer | Bioconductor | https://bioconductor.org/packages/rtracklayer |
| R / tidyverse | CRAN | https://www.tidyverse.org/ |

---

## Data Availability

Raw sequencing data and processed methylation ratio files are deposited at [repository link]. The zebra finch reference genome assembly `GCF_003957565.2` (bTaeGut1.4.pri) is available through NCBI at https://www.ncbi.nlm.nih.gov/assembly/GCF_003957565.2/.

---

## Citation

If you use this workflow, please cite:  
> [Author list]. DNA Methylation of Cell Adhesion Genes is Associated with Habituation of the Song Response in Zebra Finches. *[Journal]*, [Year]. DOI: [DOI]
