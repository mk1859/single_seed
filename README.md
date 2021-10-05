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
``` R
data_timecourse <- import_counts ("/matrix/timecourse/", header = TRUE)
```

Data for *dog1-4* and Col-0 single seed experiment.
``` R
data_dog1 <- import_counts ("/matrix/dog1/", header = TRUE)
```

Data for *dog1-4* and Col-0 seed pool experiment.
``` R
data_dry_dog1 <- import_counts ("/matrix/dog1_htseq/", header = FALSE)
```

Our library preparation protocol is design to detect mRNAs. To filter out non-protein coding genes we need reference file with information about gene types.

``` R
Araport <- read.csv ("Araport.txt", sep = "\t", header = TRUE)

head (Araport)
```
```
       gene           type
1 AT1G01010 protein_coding
2 AT1G01020 protein_coding
3 AT1G01030 protein_coding
4 AT1G01040 protein_coding
5 AT1G01050 protein_coding
6 AT1G01060 protein_coding
```

## Pre-filter single seed data

Simirally to single cell experiments our count data is quite sparse. We need to clean our data by:
1) removal of non-protein coding genes
2) removel of genes encoded in organells
3) removal of summary lines at last rows of count matrix
4) filtering out genes with low count number
5) filtering seeds with few reads

To do that we creted finction prefilter_matrix and appllayed it to our single seed matrices.
We require mean expression of gene to be at least 1 read per seed for gene to be included and at least 5,000 reads per seed for seed to be included.

``` R
# time-course experiment
filtered_timecourse <- prefilter_matrix (data_timecourse, mean_exp=1, n_reads=5000)

# dog1 experiment
filtered_dog1 <- prefilter_matrix (data_dog1, mean_exp=1, n_reads=5000)
```
``` R
# time-course experiment
dim (filtered_timecourse) # genes / seeds remaining
[1] 8687  659

# dog1 experiment
dim (filtered_dog1) # genes / seeds remaining
[1] 10989   382
```

To plot number of sequenced reads and identified genes per seed in matrix, we created function that require also defined order of treatments as input.

time-course experiment
``` R
nreads_plot (filtered_timecourse, order= c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"))
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/nreads_timecourse.png" width=50% height=50%>


dog1 experiment
``` R
nreads_plot (filtered_dog1, order = c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h"))
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/nreads_dog1.png" width=50% height=50%>

## Filter seeds for number of background reads

As it is visible above our libraries differ in number of identified genic reads. One source of that may be different quality of libraries.
To estimate their quality we decied to count fraction of off-target reads (not in protein coding genes) for each seed by using created function.

Function uses raw and pre-filtered matices as input.

``` R
background_timecourse <- background_reads (data_timecourse, filtered_timecourse)

background_dog1 <- background_reads (data_dog1, filtered_dog1)
```

Then we made plots for number of sequenced reads and background level. Again we created custome function for that.

time-course experiment
``` R
background_plot (filtered_timecourse, order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"), 
                 background = background_timecourse)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/background_timecourse.png" width=50% height=50%>

dog1 experiment
``` R
background_plot (filtered_dog1, order = c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h"), 
                 background = background_dog1)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/background_dog1.png" width=50% height=50%>

Presence of background reads may imply that some counts atrributed to genes may be generated by technical problems.
Closer examination of read tracks in browser showed that distribution of background reads is not random and they tend to create hot spots laying both between and in genes.
In addition strength of "true" peaks is negatievly correlated with number of background reads.
Based on this observations we decided to remove from our analysis genes which read count number is strongly postively correlated with number of background reads among the seeds.
As gene expression pattern is different between treatments, we calculated that correlations for each treatment separetly as well as for all seeds. To do that we created function called correlation_table.

``` R
correlation_timecourse <- correlation_table (filtered_timecourse, background_timecourse)

correlation_dog1 <- correlation_table (filtered_dog1, background_dog1)
```

We created histograms to show fraction of genes showing correlation to background reads above treshold of 0.3 using created function. 

time-course experiment
``` R
corr_hist (correlation_timecourse, treshold = 0.3, order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS", "all"))
```


dog1 experiment
``` R
corr_hist (correlation_dog1, treshold = 0.3, order = c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h", "all"))
```
