Here we provide code for the preparation of figures for "Transcriptional heterogeneity of single seeds during the secondary dormancy induction".

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
library (scales)
library (matrixStats)
```

## Import the data

On GEO we deposited results of three different experiments. 
1) single seed RNA-seq during the time course of secondary dormancy induction. There were 6 time points (plus control treatment for one of them) each consisting of 3 libraries. We provide a matrix of read counts for each library in which rows are genes and columns are one of 32 seeds used for library preparation.

2) single seed RNA-seq for *dog1-4* mutant and Col-0 in two time points of secondary dormancy induction. Each of the four treatments (genotype + time point) consist of 3 libraries and count matrices have the same structure as above.

3) RNA-seq results for the analysis of *dog1-4* and Col-0 dry seed pools (3 biological replicas). Here we provide read counts as separate files for each of the replicas. 

We created the function called import_counts to upload the data and combine matrices into a single matrix.

Data for the time-course experiment.
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
colnames (data_dry_dog1) <- gsub( ".txt","", list.files("/home/rstudio/matrix/dog1_htseq/"))
data_dry_dog1 <- data_dry_dog1[-grep ("__", rownames (data_dry_dog1)),]
```

Our library preparation protocol was designed to detect mRNAs. To filter out non-protein-coding genes we needed a reference file with information about gene types.

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

## Pre-filtering single seed data

Similarly to single-cell experiments, our count data is sparse. We needed to clean it by:
1) removing of non-protein-coding genes
2) removing of genes encoded in organelles
3) removing of summary lines at last rows of the count matrix
4) filtering out genes with a low count number
5) filtering seeds with not enough reads

To do that we created function prefilter_matrix and applied it to our single seed matrices. By default it uses Araport data frame with columns described as above.
We require the mean expression of a gene to be at least 1 read per seed for a gene to remain and at least 5,000 reads per seed for a seed to remain.

``` R
# time-course experiment
filtered_timecourse <- prefilter_matrix (data_timecourse, mean_exp=1, n_reads=5000)

# dog1-4 experiment
filtered_dog1 <- prefilter_matrix (data_dog1, mean_exp=1, n_reads=5000)
```
``` R
# time-course experiment
dim (filtered_timecourse) # genes / seeds remaining
[1] 8687  659

# dog1-4 experiment
dim (filtered_dog1) # genes / seeds remaining
[1] 10989   382
```

We wrote a function to plot the number of sequenced reads and identified genes per seed. We wanted to show treatments in the specified order.

time-course experiment
``` R
timepoints <- c ("SD1h","SD1d","SD3d","SD5d","SD7d","SD7d24h","SD7dPS")

nreads_plot (filtered_timecourse, order= timepoints)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/nreads_timecourse.png" width=50% height=50%>


*dog1-4* experiment
``` R
treatments <- c ("SD_Col0_3d","SD_dog1_3d","SD_Col0_7d24h","SD_dog1_7d24h")

nreads_plot (filtered_dog1, order = treatments)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/nreads_dog1.png" width=50% height=50%>

## Filter seeds for a number of background reads

As visible on the plots above, our libraries vary in the number of identified genic reads. One source of this may be the different quality of sequenced libraries reflected by the ratio of target and off-target reads.
We counted the fraction of off-target reads (not in protein-coding genes) for each seed with the background_reads function which uses raw and pre-filtered matrices as input.

``` R
background_timecourse <- background_reads (data_timecourse, filtered_timecourse)

background_dog1 <- background_reads (data_dog1, filtered_dog1)
```

Then, we made plots for the number of target reads and fraction of background reads with the background_plot function.

time-course experiment
``` R
background_plot (filtered_timecourse, order = timepoints, background = background_timecourse)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/background_timecourse.png" width=50% height=50%>

*dog1-4* experiment
``` R
background_plot (filtered_dog1, order = treatments, background = background_dog1)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/background_dog1.png" width=50% height=50%>

The abundance of background reads may imply that some counts attributed to genes may not reflect their expression.
Closer examination of read tracks in the browser showed that the distribution of background reads is not random and they tend to create hot spots laying both between genes and partially overlapping with them. In addition, the strength of genic peaks is negatively correlated with the number of background reads.
Based on these observations, we decided to remove from our analysis genes whose read count is strongly positively correlated with the number of background reads. As gene expression patterns are different between treatments, we calculated these correlations for each of them separately as well as for all seeds combined. To do that we wrote the function called correlation_table.

``` R
correlation_timecourse <- correlation_table (filtered_timecourse, background_timecourse)

correlation_dog1 <- correlation_table (filtered_dog1, background_dog1)
```

We created histograms to show the fraction of genes with correlation to background reads above a threshold of 0.3. 

time-course experiment
``` R
corr_hist (correlation_timecourse, threshold = 0.3, order = timepoints)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/histograms_timepoint.png" width=50% height=50%>

*dog1-4* experiment
``` R
corr_hist (correlation_dog1, threshold = 0.3, order = treatments)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/histograms_dog1.png" width=50% height=50%>

Some genes could be correlated to the background fraction in several treatments while others only if all seeds were considered. We created plots showing how many times each gene was found significantly correlated.

time-course experiment
``` R
corr_number (correlation_timecourse, threshold = 0.3)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/area_timecourse.png" width=10% height=10%>

*dog1-4* experiment
``` R
corr_number (correlation_dog1, threshold = 0.3)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/area_dog1.png" width=10% height=10%>

From pre-filtered matrices, we removed genes showing read count correlation to background reads fraction higher than 0.3 in any treatment or when all seeds were considered.

``` R
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

After we obtained filtered matrices of counts for both single seed experiments, we created Seurat objects (REF) with sctransform normalization (REF). To do this, we prepared a wrapper function that takes the count matrix and extracts information about seeds from their names. Optionally fraction of background reads may be used as a variable to regress.

Time-course experiment. Here, due to high background content in some libraries, we regressed the fraction of background reads.
``` R
seurat_timecourse <- seurat_object (filtered_timecourse, background = background_timecourse)
```

*dog1-4* experiment. Here, we did not regress the fraction of background reads.
``` R
seurat_dog1 <- seurat_object (filtered_dog1, background = background_dog1, include_background = FALSE)
```

We calculated PCA during the preparation of Seurat objects. Now, we plotted it to show treatments and batches (libraries) with the pca_discrete function.
This function exports dimension reduction and metadata from the Seurat object. It is possible to choose color pallet from ggthemes and exclude some treatments from the plot.

time-course experiment
``` R
pca_discrete (seurat_timecourse, type = "timepoint", tableu = "Tableau 10")

pca_discrete (seurat_timecourse, type = "batch", tableu = "Tableau 20", excluded ="SD7dPS")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/pca_timepoint_timecourse.png" width=40% height=40%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_batch_timecourse.png" width=40% height=40%>

*dog1-4* experiment
``` R
pca_discrete (seurat_dog1, type = "timepoint", tableu = "Tableau 10")

pca_discrete (seurat_dog1, type = "batch", tableu = "Tableau 20")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/pca_treatment_dog1.png" width=40% height=40%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_batch_dog1.png" width=40% height=40%>

Some technical parameters like number of reads, number of identified genes and fraction of background reads may affect the position of seeds on the PCA plot.
To check continuous values on PCA plots we wrote another plotting function.

time-course experiment
``` R
pca_continuous (seurat_timecourse, column = "log10_reads")

pca_continuous (seurat_timecourse, column = "n_gene")

pca_continuous (seurat_timecourse, column = "background")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/pca_reads_timecourse.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_genes_timecourse.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_background_timecourse.png" width=33% height=33%>

*dog1-4* experiment
``` R
pca_continuous (seurat_dog1, column = "log10_reads")

pca_continuous (seurat_dog1, column = "n_gene")

pca_continuous (seurat_dog1, column = "background")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/pca_reads_dog1.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_genes_dog1.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/pca_background_dog1.png" width=33% height=33%>

## Differential gene expression

To compare gene expression changes between two sets of conditions, we created a wrapper function to Seurat FindMarkers.

For the time-course experiment, we compared sequential time points.
``` R
deg_timecourse <- deg_list (seurat_timecourse, vector1 = c ("SD1h","SD1d","SD3d","SD5d","SD7d"), 
                           vector2 = c ("SD1d","SD3d","SD5d","SD7d","SD7d24h"), 
                           column = "timepoint", padj = 0.05, log2FC_threshold = 1)
```                           

For the *dog1-4* experiment, we compared mutant and wild type in two time points.
``` R
deg_dog1 <- deg_list (seurat_dog1, vector1 = c ("SD_Col0_3d","SD_Col0_7d24h"), 
                           vector2 = c ("SD_dog1_3d","SD_dog1_7d24h"), 
                           column = "timepoint", padj = 0.05, log2FC_threshold = 1)
``` 

We plotted the number of differentially expressed genes (DEG)s with the division for upregulated and downregulated genes using the deg_plot function.
It takes the parameters of the plot as input.

``` R
deg_plot (deg_timecourse, direction = TRUE, limits = c(-400,800), by = 200)
``` 
<img src="https://github.com/mk1859/single_seed/blob/main/images/number_degs.png" width=50% height=50%>

We analyzed GO terms enriched among affected genes and plotted them with a set of interdepended functions:
go_res, multiple_category, multiple_genelists, go_heatmap, go_bubble and go_pca_map.

We noticed that in the time-course experiment, three comparisons are significantly enriched in genes involved in translation.

``` R
translation_list <- list()
translation_list$SD1d_up <- rownames (deg_timecourse$SD1h_SD1d [deg_timecourse$SD1h_SD1d$avg_log2FC > 0,])
translation_list$SD3d_down <- rownames (deg_timecourse$SD1d_SD3d [deg_timecourse$SD1d_SD3d$avg_log2FC < 0,])
translation_list$SD7d24h_up <- rownames (deg_timecourse$SD7d_SD7d24h [deg_timecourse$SD7d_SD7d24h$avg_log2FC > 0,])

# function to create a data frame with enriched go terms for lists of genes 
# using rrvgo package it also generalizes results by finding parent GO terms

translation_go <- multiple_genelists (translation_list, background = rownames(filtered_timecourse), 
                                      p_value = 0.05, rrvgo_threshold=0.99)
                                      
# plot heatmap of p-values for enriched GO terms, with GO id instead of GO term names 
# with term ontology and parent GO terms indicated as neighbouring plots                                     

go_heatmap (translation_go, term_name =FALSE, term_category = TRUE, parent_term = TRUE)
``` 
<img src="https://github.com/mk1859/single_seed/blob/main/images/translation_heatmap.png" width=50% height=50%>

We decided to plot the expression of genes belonging to the translation GO term on the PCA map.
The function allows selecting column with treatments and excluding one of them. Values of expression for treatments are plotted as violin plots in the indicated order.

``` R
go_pca_map (seurat_timecourse, go_id = "GO:0006412", excluded = NULL, column = "timepoint", order =  timepoints)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/translation_pca.png" width=50% height=50%>

We also observed an interesting set of enriched GO terms for genes downregulated upon transition from 1h to 1d time point.

``` R
dry_list <- list()
dry_list$SD1d_down <- rownames (deg_timecourse$SD1h_SD1d [deg_timecourse$SD1h_SD1d$avg_log2FC < 0,])

dry_go <- multiple_genelists (dry_list, background = rownames(filtered_timecourse), 
                                   p_value = 0.05, rrvgo_threshold=0.99)

go_heatmap (dry_go, term_name =TRUE, term_category = FALSE, parent_term = FALSE)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/dry_heatmap.png" width=33% height=33%>

## Gene expression patterns

We created a function to plot normalized expression of a gene on PCA plot with violin plot inset to show its expression in treatments selected by the variable called column.
``` R
gene_exp (seurat_timecourse, gene = "AT2G36270", order = timepoints, column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/abi5_plot.png" width=50% height=50%>

We identified 500 most variably expressed genes in the time-course experiment and looked for GO terms enriched among them in BP ontology.
``` R
hvg_timecourse <- hvg (seurat_timecourse, top = 500)

hvg_go <- go_res (hvg_timecourse$gene, background = rownames(filtered_timecourse), 
                  p_value = 0.05, category = "BP", rrvgo_threshold=0.8)

go_bubble (hvg_go)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/hvg_go.png" width=50% height=50%>

We also identified co-expressed gene groups using the coexpressed function. It calculates pairwise gene expression correlations and filters them. Next, it creates a graph object, looks for highly connected groups in it and outputs groups with gene numbers above the set threshold.

time-course experiment
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
<img src="https://github.com/mk1859/single_seed/blob/main/images/cluster1_go.png" width=50% height=50%>

*dog1-4* experiment
``` R
clusters_dog1 <- coexpressed (seurat_dog1, threshold = 0.6, n_genes = 10)

lengths (clusters_dog1)
```
```
cluster_1 cluster_2 cluster_3 cluster_4 cluster_5 
      385       154        39        27        21
```

Identified co-expressed gene groups were used to create signatures that were plotted on PCA maps.

time-course experiment
``` R
seurat_timecourse <- AddModuleScore(seurat_timecourse, features = clusters_timepoint, name = "cluster_")

signature_map (seurat_timecourse, signature = "cluster_1", order = timepoints, column = "timepoint")

signature_map (seurat_timecourse, signature = "cluster_2", order = timepoints, column = "timepoint")

signature_map (seurat_timecourse, signature = "cluster_3", order = timepoints, column = "timepoint")

signature_map (seurat_timecourse, signature = "cluster_4", order = timepoints, column = "timepoint")

signature_map (seurat_timecourse, signature = "cluster_5",  rder = timepoints, column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/cluster1_pca.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/cluster2_pca.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/cluster3_pca.png" width=33% height=33%>
<img src="https://github.com/mk1859/single_seed/blob/main/images/cluster4_pca.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/cluster5_pca.png" width=33% height=33%>

*dog1-4* experiment
``` R
seurat_dog1 <- AddModuleScore(seurat_dog1, features = clusters_timepoint, name = "old_")

seurat_dog1 <- AddModuleScore(seurat_dog1, features = clusters_dog1, name = "new_")

signature_map (seurat_dog1, signature = "new_1", order = treatments, column = "timepoint")

signature_map (seurat_dog1, signature = "new_2", order = treatments, column = "timepoint")

signature_map (seurat_dog1, signature = "old_2", order = treatments, column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/new1_pca.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/new2_pca.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/old2_pca.png" width=33% height=33%>

We looked for gene overlaps between identified groups in two single seed experiments using Venn diagrams.
``` R
# plotting Venn diagrams
plot <- list(clusters_timepoint [[1]], clusters_dog1 [[1]])
names(plot) <- c("time course cluster 1","dog1-4 cluster 1")
plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF"))

plot <- list(clusters_timepoint [[2]], clusters_dog1 [[2]])
names(plot) <- c("time course cluster 2","dog1-4 cluster 2")
plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF"))
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/cluster1_venn.png" width=33% height=33%> <img src="https://github.com/mk1859/single_seed/blob/main/images/cluster2_venn.png" width=33% height=33%>

We finally plotted levels of two signatures using the sig_vs_sig function.

time-course experiment
``` R
sig_vs_sig (seurat_timecourse, "cluster_1", "cluster_2", exclude= "SD7dPS", order = timepoints, column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/timecourse_sign.png" width=50% height=50%>

*dog1-4* experiment
``` R
sig_vs_sig (seurat_dog1, "new_1", "new_2", order = treatments, column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/dog1_sign.png" width=50% height=50%>

## Gene expression variability

We looked for gene expression variability in each of the time points of the time-course experiment.

First, we created a Seurat object for each of the time points. To avoid including lowly expressed genes, we filtered them for each time point separately.
``` R
selected_genes <- select_genes(filtered_timecourse, treatments = timepoints, avg_reads = 1))
timepoint_seurats <- list_seurat (selected_genes, background = background_timecourse)
```

To estimate gene expression variability in each time point we calculated the expression variance of the top 200 variable genes.
``` R
hvg_timepoints <- lapply (timepoint_seurats, function (x) hvg (x, top = 200)) 

hvg_timepoints <- rbindlist(hvg_timepoints, idcol = "timepoint")

# set order of time point
hvg_timepoints$timepoint <- factor(hvg_timepoints$timepoint, levels = timepoints)

# plot variance
ggplot(hvg_timepoints, aes(x=timepoint, y= log10(var), color = timepoint)) +
  geom_boxplot() +
  theme_classic() + 
  scale_color_tableau() +
  ylab ("log10(variance)") +
  theme(legend.position = "none")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/hvg_boxplot.png" width=50% height=50%>

The variance of genes expression does not say anything if gene expression variability is random or create some patterns among seeds.
To find how much seeds differ in each time point, we divided them into sub-pools and performed differential gene expression analysis between them.
``` R
# Seurat function FindNeighbors
timepoint_seurats <- lapply (timepoint_seurats, function (x) {
                              FindNeighbors(object = x, dims = 1:10, verbose = FALSE)})

# resolutions to obtain exactly two sub-pools were found earlier manually
resolution <- c ("SD1h" = 0.7, "SD1d" = 0.6, "SD3d" = 0.8, "SD5d" = 0.8, "SD7d" = 0.85, 
                 "SD7d24h" = 0.3, "SD7dPS" = 0.85)

# finding seed clusters with Seurat function
timepoint_seurats <- mapply (function (x, y) {
                            FindClusters(object = x, verbose = FALSE, resolution = y)},
                            timepoint_seurats, resolution)

# finding DEGs with Seurat function
degs_subpools_timepoints <- lapply (timepoint_seurats, function (x) {
                                  FindMarkers(object = x, ident.1 = 0, ident.2 = 1, 
                                              logfc.threshold = log2(1.5), 
                                              test.use = "wilcox",only.pos =FALSE, 
                                              assay = "SCT", slot ="data", verbose = FALSE) %>%
                                  filter(.,p_val_adj < 0.05)}) 
                           
names(degs_subpools_timepoints) <- timepoints

# plot for number of DEGs
deg_plot (degs_subpools_timepoints, direction = FALSE, limits = c(0,500), by = 100)

# plot for seed sub-pools (export data from Seurat objects, combine them and add time point column)
plot <- lapply (timepoint_seurats, function (x) { 
                cbind (as.data.frame (Embeddings(object = x, reduction = "pca")) [,1:2], 
                x@meta.data [,-8])}) %>% 
        rbindlist (.)  %>% 
        mutate (., timepoint = factor(gsub('.{0,3}$', '', batch), levels = timepoints))

ggplot(plot, aes(x=PC_1, y= PC_2, color = seurat_clusters, group = timepoint)) +
  geom_point (size = 3) + 
  scale_color_tableau() +
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(vars(timepoint), scales = "free")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/clusters_each.png" width=50% height=50%>
<img src="https://github.com/mk1859/single_seed/blob/main/images/degs_each_timepoint.png" width=50% height=50%>

We identified cluster_1 and cluster_2 gene groups with largely antagonistic expression during the time course. We wanted to check if their expression pattern can distinguish seeds in each time point. AddModuleScore from Seurat does not allow creating complex gene expression signatures with some genes having positive and some negative input. To create such signature we the used Vision package (REF).

``` R
# create gene expression signature for Vision
germ_comp <-  c(rep (1, length (clusters_timepoint$cluster_1)), rep (-1, length (clusters_timepoint$cluster_2)))
germ_comp <- setNames (germ_comp, c(clusters_timepoint$cluster_1,clusters_timepoint$cluster_2))
germ_comp <- createGeneSignature(name = "germ_comp", sigData = germ_comp)

# make Vision objects
vis <- lapply (timepoint_seurats, function(x) {
                  Vision(x, signatures = list(germ_comp), 
                         meta = x@meta.data, assay = "SCT")  
        })

# analyze Vision objects
vis <- lapply (vis, function(x) {
                    analyze(x)  
        })

# make plot data frame, Vision signatures were scaled to better illustrate gradients between seeds
plot <- mapply (function (x,y) {
                    cbind (as.data.frame (Embeddings(object = x, reduction = "pca"))  [,1:2], 
                           x@meta.data[,c(4,9)], rescale(y@SigScores, to = c(-1, 1)))
      }, timepoint_seurats, vis, SIMPLIFY = FALSE)

plot <- rbindlist(plot)
plot$treatment <- factor(gsub('.{0,3}$', '', plot$batch), levels = timepoints)

# make plot of PCA maps
ggplot(plot, aes(x=PC_1, y= PC_2, color = germ_comp, group = treatment)) +
  geom_point (size = 3) + 
  scale_color_viridis() +
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(vars(treatment), nrow = 2, scales = "free")

# make boxplots
ggplot(plot, aes(x=seurat_clusters, y= germ_comp, color = seurat_clusters)) +
  geom_boxplot () + 
  scale_color_tableau() +
  theme_classic ()  +
  theme (legend.position = "none") +
  facet_wrap(vars(treatment), nrow = 1)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/sign_pca.png" width=50% height=50%>
<img src="https://github.com/mk1859/single_seed/blob/main/images/sign_boxplots.png" width=50% height=50%>

Finally, identification of groups of co-expressed genes in time point may point to presence of coherent gene expression patterns among seeds.
``` R
timepoints_clusters <- lapply (timepoint_seurats, function (x) coexpressed (x, threshold =0.5, n_genes = 10))

# number of genes in co-expressed groups
number_coexp <- lapply (timepoints_clusters, lengths )

# set 0 for empty lists
number_coexp [lengths(number_coexp) == 0] <- 0                                  

# convert each list to data frame
number_coexp <- lapply (number_coexp, function (x) {                            
                                        as.data.frame(x) %>%
                                        mutate (., cluster = rownames(.))
                                        })

# join data frames
number_coexp <- rbindlist (number_coexp, idcol = "timepoint")                   
number_coexp$timepoint <- factor(number_coexp$timepoint, levels = timepoints) 

ggplot (number_coexp, aes (x = timepoint, y = x, fill= cluster)) +
                geom_bar(stat='identity') +
                theme_classic() +
                scale_fill_tableau("Miller Stone") +
                ylab ("genes")
```

## *dog1-4* vs Col-0 dry seed pool DEGs

We sequenced mRNAs isolated from pools of dry seeds to eluciadate role of DOG1 in establishment of gene expression patterns.
We identified DEGsusing DESeq2 (REF).
``` R
# metadata
col_data <- data.frame (library = colnames(data_dry_dog1),
                        genotype = as.factor(substr(colnames(data_dry_dog1), 1, 4)),
                        replica = substr(colnames(data_dry_dog1), 7, 7))

# finding DEGs with DESeq2
dds <- DESeqDataSetFromMatrix(countData = data_dry_dog1, colData = col_data, design = ~ genotype) %>%
       DESeq(.)

deg_dry_dog1 <- as.data.frame(results(dds, alpha = 0.05, contrast= c("genotype","dog1","Col0")))

# creating Volcano plot
ggplot(deg_dry_dog1 , aes(y=-log10(padj), x= log2FoldChange, color = padj < 0.05 , alpha = padj < 0.05)) +
  geom_point(size = 1) + 
  scale_color_tableau() +
  theme_classic() +
  scale_alpha_ordinal(range = c(0.1, 1))
```

We looked for ovrlaps of identified DEGS:
``` R
affected_dog1 <- rownames(deg_dry_dog1[which(deg_dry_dog1$padj< 0.05),])

# with main co-expressed gene groups of time-course experiment
plot <- list(clusters_timepoint [[1]], clusters_timepoint [[2]], affected_dog1)
names(plot) <- c("time course cluster 1","time course cluster 2", "dry seeds affected")
plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF", "#868686FF"))

# with single seed DEGs of dog1 experiment
plot <- list(rownames(deg_dog1$SD_Col0_3d_SD_dog1_3d), 
             rownames(deg_dog1$SD_Col0_7d24h_SD_dog1_7d24h), affected_dog1)
names(plot) <- c("dog1 3d","dog1 7d24h", "dry seeds affected")
plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF", "#868686FF"))
```

We identified GO terms enriched among genes with down- and upregulated expression.
``` R
# find enriched GO terms
up_genes <- rownames (deg_dry_dog1[which(deg_dry_dog1$padj< 0.05 & deg_dry_dog1$log2FoldChange > 0),])
dw_genes <- rownames (deg_dry_dog1[which(deg_dry_dog1$padj< 0.05 & deg_dry_dog1$log2FoldChange < 0),])
background <- rownames (deg_dry_dog1[which(deg_dry_dog1$baseMean > 1),])

up_genes <- gost(query = up_genes, organism = "athaliana", 
            custom_bg = background, user_threshold = p_value,
            sources = "GO")$result

dw_genes <- gost(query = dw_genes, organism = "athaliana", 
            custom_bg = background, user_threshold = p_value,
            sources = "GO")$result

# create data frame
affected <- list(up= up_genes, dw= dw_genes)
affected <- rbindlist(affected, idcol = "change")

# sort data frame
affected <- affected [order(affected$p_value, decreasing = TRUE),]
lev_order <- as.factor(affected$term_name [!duplicated(affected$term_name)])
affected$term_name <- factor(affected$term_name,levels = lev_order)

# plot for enrichment
g1 <- ggplot(affected, aes(1, term_name)) + 
         geom_tile(aes(fill = -log10(p_value))) + 
         scale_fill_gradientn(colors =c("white","darkred"), limits= c( 0, -log10(min(affected$p_value))))+
         theme_classic() +
         facet_grid(vars(change), scales = "free", space = "free_y") +
         theme(axis.line.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.x=element_blank()) 

# plot for type of ontologies
g2 <- ggplot(affected, aes(1, term_name)) + 
         geom_tile(aes(fill = source)) + 
         theme_classic() +
         facet_grid(vars(change), scales = "free", space = "free_y") +
         theme(axis.line=element_blank(),
               axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank()) 
 
# combine plots
plot_grid(g1, g2, nrow = 1, align = "h", 
               rel_widths = c(2, 0.5))
```

Finally we crated signature of *dog1-4* affected genes using Vision and overlaid it on time-course experiment PCA map.

``` R
# create signature 
deg_dry_dog1 <- as.data.frame(deg_dry_dog1) %>% filter (., padj < 0.05) 
dog1_sign <-  sign(deg_dry_dog1$log2FoldChange)
dog1_sign <- setNames (dog1_sign, rownames(deg_dry_dog1))
dog1_sign <- createGeneSignature (name = "dog1mut_sign", sigData = dog1_sign)

# create Vision object
vis <- Vision(seurat_timecourse, signatures = list(dog1_sign), meta = seurat_timecourse@meta.data, assay = "SCT")
vis <- analyze(vis)

# create plot
seurat_timecourse@meta.data$dog1_4 <- vis@SigScores
signature_map (seurat_timecourse, signature = "dog1_4", order = timepoints, column = "timepoint")
```
