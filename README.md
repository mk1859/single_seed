Here we provide code for preparation of figures for "Transcriptional heterogeneity of single seeds during the secondary dormancy induction".

#### Load required R packages.

```
library (tidyverse)
library (ggthemes)
library (ggbeeswarm)
library (viridis)
library (sctransform)
library (Seurat)
library (rlist)
library (cowplot)
library (patchwork)
library (gprofiler2)
library (rrvgo)
library (scran)
library (scater)
library (graph)
library (RBGL)
library (data.table)
library (eulerr)
library (DESeq2)
library (VISION)
library (org.At.tair.db)
library (biomaRt)
library (GO.db)
```

## Import the data

On GEO we deposited results of three different experiments. 
1) single seed RNA-seq during time course of secondary dormancy induction. There was 6 time points (plus control treatment for one of them) each consisting of 3 libraries. We provide matrix of read counts for each libarary in which rows are genes and columns are one of 32 seeds used for library preparation.

2) single seed RNA-seq for *dog1-4* and Col-0 seeds in two time points of secondary dormancy induction. Again, each of four tratments (genotype + time point) consist of 3 libraries of the same structure as above.

3) RNA-seq results for analysis of *dog1-4* and Col-0 dry seed pools (3 biological replicas). Here we provide read counts as separate files for each of replicas. 

We created the function called import_counts to upload the data and combine matrices into single matrix.

Data for time-course experiment.
```
data_timecourse <- import_counts ("matrix/timecourse/", header = TRUE)
```

Data for *dog1-4* and Col-0 single seed experiment.
```
data_dog1 <- import_counts ("matrix/dog1/", header = TRUE)
```

Data for *dog1-4* and Col-0 seed pool experiment.
```
data_dry_dog1 <- import_counts ("matrix/dog1_htseq/", header = FALSE)
```

Our library preparation protocol is design to detect mRNAs. To filter out non-protein coding genes we need reference file with information about gene types.

```
head (Araport)
  chr    source feature start   end score strand frame      gene           type
1   1 Araport11    gene  3631  5899     .      +     . AT1G01010 protein_coding
2   1 Araport11    gene  6788  9130     .      -     . AT1G01020 protein_coding
3   1 Araport11    gene 11649 13714     .      -     . AT1G01030 protein_coding
4   1 Araport11    gene 23121 31227     .      +     . AT1G01040 protein_coding
5   1 Araport11    gene 31170 33171     .      -     . AT1G01050 protein_coding
6   1 Araport11    gene 33365 37871     .      -     . AT1G01060 protein_coding
```
