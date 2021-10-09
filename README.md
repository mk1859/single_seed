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
<img src="https://github.com/mk1859/single_seed/blob/main/images/histograms_timepoint.png" width=50% height=50%>

dog1 experiment
``` R
corr_hist (correlation_dog1, treshold = 0.3, order = c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h", "all"))
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/histograms_dog1.png" width=50% height=50%>

As genes can be correlated to the background in severla treatments we created a plot showing how many genes were found correlated.

time-course experiment
``` R
corr_number (correlation_timecourse, treshold = 0.3)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/area_timecourse.png" width=10% height=10%>

dog1 experiment
``` R
corr_number (correlation_dog1, treshold = 0.3)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/area_dog1.png" width=10% height=10%>

Fianally we remove from our pre-filtered matrices genes showing read count correlation to background reads higher than 0.3 in any treatment or when all seeds considered.

``` R
library (matrixStats)
filtered_timecourse <- filtered_timecourse [-which (rowMaxs (correlation_timecourse) > 0.3),]
nrow (filtered_timecourse) # genes remaining
```
```
[1] 6430
```
``` R
filtered_dog1 <- filtered_dog1 [-which (rowMaxs (correlation_dog1) > 0.3),]
nrow (filtered_dog1) # genes remaining
```
```
[1] 9110
```

## Seurat object

We obtained final matrices of counts for both single seed experiments. Now we prepare Seurat object (REF) with sctransform normalization (REF).
To do that we prepared wraper function that takes matrix and extracts information about seeds from their names. Optionally fraction of background reads may be used as varaible to regress.

Time-course experiment. Here, due to high background content in some libraries we regressed fraction of background reads.
``` R
seurat_timecourse <- seurat_object (filtered_timecourse, background = background_timecourse)
```

DOG1 experiment. Here, we did not regress fraction of background reads.
``` R
seurat_dog1 <- seurat_object (filtered_dog1, background = background_dog1, include_background = FALSE)
```

We calculated PCA during preparation of Seurat objects. Now we ploted it to show tratments and batches (libraries).
To do this we prepared function exporting dimention reduction and metadata from Seurat object. 
It is possible to chose color pallet from ggthemes and exclude some treatments from the plot.

time-course experiment
``` R
pca_discrete (seurat_timecourse, type = "timepoint", tableu = "Tableau 10")

pca_discrete (seurat_timecourse, type = "batch", tableu = "Tableau 20", excluded ="SD7dPS")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/pca_timepoint_timecourse.png" width=40% height=40%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_batch_timecourse.png" width=40% height=40%>

dog1 experiment
``` R
pca_discrete (seurat_dog1, type = "timepoint", tableu = "Tableau 10")

pca_discrete (seurat_dog1, type = "batch", tableu = "Tableau 20")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/pca_treatment_dog1.png" width=40% height=40%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_batch_dog1.png" width=40% height=40%>

Some technical parameters like number of reads, number of identified genes and fraction of backround reads per seed may affect position of seeds on the map.
To check continous values we crate another plotting function.

time-course experiment
``` R
pca_continuous (seurat_timecourse, column = "log10_reads")

pca_continuous (seurat_timecourse, column = "n_gene")

pca_continuous (seurat_timecourse, column = "background")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/pca_reads_timecourse.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_genes_timecourse.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_background_timecourse.png" width=33% height=33%>

dog1 experiment
``` R
pca_continuous (seurat_dog1, column = "log10_reads")

pca_continuous (seurat_dog1, column = "n_gene")

pca_continuous (seurat_dog1, column = "background")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/pca_reads_dog1.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_genes_dog1.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_background_dog1.png" width=33% height=33%>

## Differential gene expression

To compare gene expression changes between two sets of conditions we created wrapper function to Seurat FindMarkers.

For time-course experiment we compared sequential time points.
``` R
deg_timecourse <- deg_list (seurat_timecourse, vector1 = c ("SD1h","SD1d","SD3d","SD5d","SD7d"), 
                           vector2 = c ("SD1d","SD3d","SD5d","SD7d","SD7d24h"), 
                           column = "timepoint", padj = 0.05, log2FC_threshold = 1)
```                           

For dog1 experiment we compared mutant and wild type in two time points.
``` R
deg_dog1 <- deg_list (seurat_dog1, vector1 = c ("SD_Col0_3d","SD_Col0_7d24h"), 
                           vector2 = c ("SD_dog1_3d","SD_dog1_7d24h"), 
                           column = "timepoint", padj = 0.05, log2FC_threshold = 1)
``` 

We plotted number of DEGs with division for upregulated and downregulated genes using custom function.
It takes parameters of the plot as input.

``` R
deg_plot (deg_timecourse, direction = TRUE, limits = c(-400,800), by = 200)
``` 
<img src="https://github.com/mk1859/single_seed/blob/main/images/number_degs.png" width=50% height=50%>

To analyze GO terms enriched among affected genes and plot them we created set of interdepended functions:
go_res.R, multiple_category.R, multiple_genelists.R, go_heatmap.R, go_bubble.R and go_pca_map.R

We noticed that in three comparisons in time-course experiment genes involved in translation are significantly enriched.

``` R
translation_list <- list()
translation_list$SD1d_up <- rownames (deg_timecourse$SD1h_SD1d [deg_timecourse$SD1h_SD1d$avg_log2FC > 0,])
translation_list$SD3d_down <- rownames (deg_timecourse$SD1d_SD3d [deg_timecourse$SD1d_SD3d$avg_log2FC < 0,])
translation_list$SD7d24h_up <- rownames (deg_timecourse$SD7d_SD7d24h [deg_timecourse$SD7d_SD7d24h$avg_log2FC > 0,])

# function to create data frame with enriched go terms in lists of genes, 
# using rrvgo package it also generalizes results by finding parent GO terms

translation_go <- multiple_genelists (translation_list, background = rownames(filtered_timecourse), 
                                      p_value = 0.05, rrvgo_threshold=0.99)
                                      
# plotting heatmap of GO terms enrichment p-values, with GO id instead of GO term names 
# with term ontology and parent GO terms indicated                                     

go_heatmap (translation_go, term_name =FALSE, term_category = TRUE, parent_term = TRUE)
``` 
<img src="https://github.com/mk1859/single_seed/blob/main/images/translation_heatmap.png" width=50% height=50%>

Due to fluctuation of translation related genes we decided to plot expression of genes belonging to translation GO term on PCA map.
Function allows to exclude treatament and select column with treatemnts that will be plotted as violin plot in indicated order.

``` R
go_pca_map (seurat_timecourse, go_id = "GO:0006412", excluded = NULL, column = "timepoint",
            order =  c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"))
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/translation_pca.png" width=50% height=50%>

We also observed intresting set of enriched GO terms for genes downregulated upon transition from 1h to 1d time point.

``` R
dry_list <- list()
dry_list$SD1d_down <- rownames (deg_timecourse$SD1h_SD1d [deg_timecourse$SD1h_SD1d$avg_log2FC < 0,])

dry_go <- multiple_genelists (dry_list, background = rownames(filtered_timecourse), 
                                   p_value = 0.05, rrvgo_threshold=0.99)

go_heatmap (dry_go, term_name =TRUE, term_category = FALSE, parent_term = FALSE)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/dry_heatmap.png" width=33% height=33%>


## Gene expression patterns

We created custome function to plot normalized expression of gene on PCA plot with violin plot inset to show its expression in treatments.
``` R
gene_exp (seurat_timecourse, gene = "AT2G36270", order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"), column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/abi5_plot.png" width=50% height=50%>

We identified 500 most varaibly expressed genes in time-course experiment and looked for GO terms enriched among them in BP ontology.
``` R
hvg_timecourse <- hvg (seurat_timecourse, top = 500)

hvg_go <- go_res (hvg_timecourse$gene, background = rownames(filtered_timecourse), 
                  p_value = 0.05, category = "BP", rrvgo_threshold=0.8)

go_bubble (hvg_go)
```


We also identified co-expressed gene groups using custome function. It calculates pairwise gene expression correlations and filters them. Next creates graph object, looks for highly connected groups in it and output groups with gene number above set treshold.

Time-course experiment
``` R
seurat_timecourse@active.ident <- as.factor(seurat_timecourse$timepoint)

clusters_timepoint <- coexpressed (subset(seurat_timecourse, idents = "SD7dPS", invert = TRUE), 
                                   threshold = 0.5, n_genes = 10)
               
lengths (clusters_timepoint)
```
```
cluster_1 cluster_2 cluster_3 cluster_4 cluster_5 
      266       123        74        36        13
```

We noticed that cluster_1 is enriched in translation-related GO terms
``` R
cluster1_go <- go_res (clusters_timepoint$cluster_1, background = rownames(filtered_timecourse), 
                   p_value = 0.05, category = "BP", rrvgo_threshold=0.8)

go_bubble (cluster1_go)
```


dog1 experiment
``` R
clusters_dog1 <- coexpressed (seurat_dog1, threshold = 0.6, n_genes = 10)

lengths (clusters_dog1)
```
```
cluster_1 cluster_2 cluster_3 cluster_4 cluster_5 
      385       154        39        27        21
```

Identified co-expressed gene groups were used to create signatures that were plotted on PCA maps.

Time-course experiment
``` R
seurat_timecourse <- AddModuleScore(seurat_timecourse, features = clusters_timepoint, name = "cluster_")

signature_map (seurat_timecourse, signature = "cluster_1", 
               order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"), column = "timepoint")

signature_map (seurat_timecourse, signature = "cluster_2", 
               order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"), column = "timepoint")

signature_map (seurat_timecourse, signature = "cluster_3", 
               order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"), column = "timepoint")

signature_map (seurat_timecourse, signature = "cluster_4", 
               order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"), column = "timepoint")

signature_map (seurat_timecourse, signature = "cluster_5", 
               order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"), column = "timepoint")
```

dog1 experiment
``` R
seurat_dog1 <- AddModuleScore(seurat_dog1, features = clusters_timepoint, name = "old_")

seurat_dog1 <- AddModuleScore(seurat_dog1, features = clusters_dog1, name = "new_")

signature_map (seurat_dog1, signature = "new_1", 
               order = c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h"), column = "timepoint")

signature_map (seurat_dog1, signature = "new_2", 
               order = c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h"), column = "timepoint")

signature_map (seurat_dog1, signature = "old_2", 
               order = c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h"), column = "timepoint")
```

We look for gene overlaps between identified groups in two expriments using Venn diagrams.
``` R
# plotting Venn diagrams
plot <- list(clusters_timepoint [[1]], clusters_dog1 [[1]])
names(plot) <- c("time course cluster 1","dog1 cluster 1")
plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF"))

plot <- list(clusters_timepoint [[2]], clusters_dog1 [[2]])
names(plot) <- c("time course cluster 2","dog1 cluster 2")
plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF"))
```


We finally plotted levels of two signatures using custom function.

Time-course experiment
``` R
sig_vs_sig (seurat_timecourse, "cluster_1", "cluster_2", exclude= "SD7dPS",
            order = c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS"),
            column = "timepoint")
```

dog1 experiment
``` R
sig_vs_sig (seurat_dog1, "new_1", "new_2",
            order = c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h"),
            column = "timepoint")
```


## Gene expression variability

We looked for gene expression variability in each of time points of our time-course experiment.

First we create Seourat object for each of time points. To avoid including lowly expressed genes we filtered them for each time point.


calculated varience of top variable genes in each time point.






