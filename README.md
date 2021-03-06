Here we provide code for the preparation of most of the figures for "Single seeds exhibit transcriptional heterogeneity during secondary dormancy induction".

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
library (UpSetR)
```

## Import the data

On GEO we deposited results of three different experiments. 
1) single seed RNA-seq during the time course of secondary dormancy induction. There were 6 time points (plus control treatment for one of them) each consisting of 3 libraries. We provide a matrix of read counts for each library in which rows are genes and columns are one of 32 seeds used for library preparation.

2) single seed RNA-seq for *dog1-4* mutant and Col-0 in two time points of secondary dormancy induction. Each of the four treatments (genotype + time point) consist of 3 libraries and count matrices have the same structure as above.

3) 3'RNA-seq results for the analysis of *dog1-4* and Col-0 dry seed pools (3 biological replicas). Here we provide read counts as separate files for each of the replicas. 

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
1) removing non-protein-coding genes
2) removing genes encoded in organelles
3) removing summary lines at last rows of the count matrix
4) filtering out genes with a low count number
5) filtering seeds with not enough reads

To do that we created the function prefilter_matrix and applied it to our single seed matrices. By default, it uses the Araport data frame with columns described above.
We require the mean expression of a gene to be at least 1 read on average per seed for a gene to remain and at least 5,000 reads per seed for a seed to remain.

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

Then, we made plots for the number of target reads and the fraction of background reads with the background_plot function.

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
A closer examination of reads' tracks in the browser showed that the distribution of background reads is not random and they tend to create hot spots laying both between genes and partially overlapping with them. In addition, the strength of genic peaks is negatively correlated with the number of background reads.
Based on these observations, we decided to remove from our analysis genes whose read count is strongly positively correlated with the number of background reads. As gene expression patterns are different between treatments, we calculated these correlations for each of them separately as well as for all seeds combined. To do that we wrote the function called correlation_table.

``` R
correlation_timecourse <- correlation_table (filtered_timecourse, background_timecourse)

correlation_dog1 <- correlation_table (filtered_dog1, background_dog1)
```

We created histograms to show the fraction of genes with correlation to background reads above a threshold of 0.3. 

time-course experiment
``` R
corr_hist (correlation_timecourse, threshold = 0.3, order = c(timepoints, "all"))
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/histograms_timepoint.png" width=50% height=50%>

*dog1-4* experiment
``` R
corr_hist (correlation_dog1, threshold = 0.3, order = c(treatments, "all"))
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

After we obtained filtered matrices of counts for both single seed experiments, we created Seurat objects with sctransform normalization. To do this, we prepared a wrapper function that takes the count matrix and extracts information about seeds from their names. Optionally fraction of background reads may be used as a variable to regress.

Time-course experiment. Here, due to high background content in some libraries, we regressed the fraction of background reads.
``` R
seurat_timecourse <- seurat_object (filtered_timecourse, background = background_timecourse)
```

*dog1-4* experiment. Here, we did not regress the fraction of background reads.
``` R
seurat_dog1 <- seurat_object (filtered_dog1, background = background_dog1, include_background = FALSE)
```

We calculated PCA during the preparation of Seurat objects. Now, we plotted it to show treatments and batches (libraries) with the pca_discrete function.
This function exports dimension reduction and metadata from the Seurat object. It is possible to choose a color pallet from ggthemes and exclude some treatments from the plot.

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

To compare gene expression changes between two sets of conditions, we created a wrapper function for Seurat FindMarkers.

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
# using rrvgo package we generalize results by finding parent GO terms

translation_go <- multiple_genelists (translation_list, background = rownames(filtered_timecourse), 
                                      p_value = 0.05, rrvgo_threshold=0.99)
                                      
# plot heatmap of p-values for enriched GO terms, with GO id instead of GO term names 
# with term ontology and parent GO terms indicated as neighbouring plots                                     

go_heatmap (translation_go, term_name =FALSE, term_category = TRUE, parent_term = TRUE)
``` 
<img src="https://github.com/mk1859/single_seed/blob/main/images/translation_heatmap.png" width=50% height=50%>

We decided to plot the expression of genes belonging to the translation GO term on the PCA map.
The function allows selecting column with treatments and excluding some of them. Values of expression for treatments are plotted as violin plots in the indicated order.

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

Germination assays after secondary dormancy induction treatment show a gradual increase in dormancy levels. We wanted to identify a subset of genes whose expression changes are correlated with that trend.

``` R
# first we identify genes whose expression changes during the treatment
deg <- deg_list (seurat_timecourse, 
                 vector1 = c ("SD1d","SD3d","SD5d"), 
                 vector2 = c ("SD7d","SD7d","SD7d"), 
                 column = "timepoint", padj = 0.05, log2FC_threshold = log2(2))

deg <- lapply (deg, function(x) mutate(x, gene = rownames(x)))

genes <- (rbindlist(deg))$gene
genes <- genes [-which(duplicated(genes))]

# we export normalized gene expression and average it for seeds representing time points
norm_reads <- as.matrix(seurat_timecourse@assays$SCT@data)

norm_genes<- data.frame(SD_1d = rowMeans(norm_reads [,grep("SD1d", colnames(norm_reads))]),
                        SD_3d = rowMeans(norm_reads [,grep("SD3d", colnames(norm_reads))]),
                        SD_5d = rowMeans(norm_reads [,grep("SD5d", colnames(norm_reads))]),
                        SD_7d = rowMeans(norm_reads [,grep("SD7d.", colnames(norm_reads))])) 

# we filter genes affected during treatment and scale their expression
norm_genes <- norm_genes [which(rownames(norm_genes)%in% genes),]

norm_genes <- as.data.frame(t(scale(t(norm_genes))))

# we cluster gene expression into 12 groups
gene_clusters <- norm_genes %>% 
  dist(.) %>%
  hclust(., method = "complete") %>%
  cutree(., k = 12) %>%
  enframe(., name = "gene", value = "cluster")
  
 # finally we prepare data frame for plotting
 norm_genes <- norm_genes %>%
  mutate (.,gene = rownames(norm_genes)) %>%
  pivot_longer(., cols = SD_1d:SD_7d, names_to = "SD", values_to = "exp") %>%
  merge (.,gene_clusters, by= "gene") 

ggplot(norm_genes, aes(SD, exp, color = as.factor(cluster))) +
  geom_line(aes(group = gene), alpha = 0.3) +
  facet_wrap(~ cluster, ncol = 2)+ 
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_color_tableau("Tableau 20")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/clustered_plots.png" width=33% height=33%>

## Gene expression patterns

We created a function to plot normalized expression of a gene on PCA plot with violin plot inset to show its expression in treatments selected by the variable called column.
``` R
gene_exp (seurat_timecourse, gene = "AT2G36270", order = timepoints, column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/abi5_plot.png" width=50% height=50%>

We identified 500 most variably expressed genes in the time-course experiment and looked for GO terms enriched among them in BP ontology.
``` R
hvg_timecourse <- FindVariableFeatures(subset(seurat_timecourse, idents = c("SD7dPS"), invert = TRUE), nfeatures = 500)@assays$SCT@var.features

hvg_go <- go_res (hvg_timecourse$gene, background = rownames(filtered_timecourse), 
                  p_value = 0.05, category = "BP", rrvgo_threshold=0.8)

go_bubble (hvg_go)
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/hvg_go_all.png" width=50% height=50%>

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

signature_map (seurat_timecourse, signature = "cluster_5",  order = timepoints, column = "timepoint")
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
sig_vs_sig (seurat_timecourse, "cluster_1", "cluster_2", order = timepoints, column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/sig_vs_sig_timecourse.png" width=50% height=50%>

*dog1-4* experiment
``` R
sig_vs_sig (seurat_dog1, "new_1", "new_2", order = treatments, column = "timepoint")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/dog1_sign.png" width=50% height=50%>

## Gene expression variability

We looked for gene expression variability in each of the time points of the time-course experiment.

First, we created a Seurat object for each of the time points. To avoid including lowly expressed genes, we filtered them for each time point separately.
``` R
selected_genes <- select_genes(filtered_timecourse, treatments = timepoints, avg_reads = 1)
timepoint_seurats <- list_seurat (selected_genes, background = background_timecourse)
```

To estimate gene expression variability in each time point we checked values of residual variance after sctransform procedure.
``` R

hvg_timepoints <- lapply (timepoint_seurats, HVFInfo)
hvg_timepoints <- rbindlist(hvg_timepoints, idcol = "timepoint")

# set order of time point
hvg_timepoints$timepoint <- factor(hvg_timepoints$timepoint, levels = timepoints)

# plot residual variance
ggplot(plot, aes(x=log10(gmean), y= residual_variance, color = residual_variance > 4)) +
  geom_point (size = 2) + 
  scale_color_tableau() +
  theme_classic () +
  facet_wrap(vars(timepoint), nrow = 1)+
  theme(strip.background = element_blank())

# select only highly variable genes
hvg_timepoints <- hvg_timepoints [which(hvg_timepoints$residual_variance > 4),]

# plot them
ggplot(hvg_timepoints, aes(x=timepoint, y= log10(residual_variance), color = timepoint)) +
  geom_jitter(size = 2) +
  theme_classic() + 
  scale_color_tableau() +
  ylab ("log10 residual_variance") 
  
# find overlaps between highly variable genes from different time points

overlap <- aggregate(gene ~ timepoint, hvg_timepoints, c)$gene
names (overlap) <- levels(as.factor(hvg_timepoints$timepoint))
upset(fromList(overlap), order.by = "freq", nsets = 7, nintersects = NA, group.by = "degree")
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/hvg_variance.png" width=90% height=90%>
<img src="https://github.com/mk1859/single_seed/blob/main/images/jitter_variance.png" width=30% height=30%> <img src="https://github.com/mk1859/single_seed/blob/main/images/upsetR_overlap.png" width=30% height=30%>

Calculate percentage of variance explained by PC1-50 in each time point.
``` R
var_exp <- lapply (timepoint_seurats, function (x){
                    (x@reductions$pca@stdev)^2/sum((x@reductions$pca@stdev)^2)*100})
var_exp <- data.frame(timepoint = rep(names(var_exp), each = 50), PC = 1:50, var = unlist(var_exp))

var_exp$timepoint <- factor(var_exp$timepoint, 
                      levels = levels(as.factor(var_exp$timepoint)) [c(2,1,3:7)])

ggplot(var_exp, aes(x=PC, y= var, color = timepoint)) +
  geom_point (size = 2) + 
  scale_color_tableau() +
  theme_classic () +
  facet_wrap(vars(timepoint), nrow = 1)+
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/PC_explained.png" width=90% height=90%>

The variance of genes expression does not say if gene expression variability is random or creates some patterns among seeds.
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
<img src="https://github.com/mk1859/single_seed/blob/main/images/degs_each_tiempoint.png" width=50% height=50%>

We identified cluster_1 and cluster_2 gene groups with largely antagonistic expression during the time course. We wanted to check if their expression pattern can distinguish seeds at each time point. AddModuleScore from Seurat does not allow creating complex gene expression signatures with some genes having positive and some negative input. To create such a signature we the used Vision package.

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

Finally, the identification of co-expressed genes groups in time point may point to the presence of coherent gene expression patterns among seeds.
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
                ylab ("genes") +
                coord_flip()
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/coexpressed_timecourse.png" width=50% height=50%>

## *dog1-4* vs Col-0 seed pool DEGs

We sequenced mRNAs isolated from pools of dry seeds to elucidate the role of DOG1 in the establishment of gene expression patterns.
We identified DEGs using DESeq2.
``` R
# metadata
col_data <- data.frame (library = colnames(data_dry_dog1),
                        genotype = as.factor(substr(colnames(data_dry_dog1), 1, 4)),
                        replica = substr(colnames(data_dry_dog1), 7, 7))

# finding DEGs with DESeq2
dds <- DESeqDataSetFromMatrix(countData = data_dry_dog1, colData = col_data, design = ~ genotype) %>%
       DESeq(.)

deg_dry_dog1 <- as.data.frame(results(dds, alpha = 0.05, contrast= c("genotype","dog1","Col0")))

# Volcano plot
ggplot(deg_dry_dog1 , aes(y=-log10(padj), x= log2FoldChange, color = padj < 0.05 , alpha = padj < 0.05)) +
  geom_point(size = 1) + 
  scale_color_tableau() +
  theme_classic() +
  scale_alpha_ordinal(range = c(0.1, 1))
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/volcano.png" width=33% height=33%>

We looked for overlaps of identified DEGs:
``` R
affected_dog1 <- rownames(deg_dry_dog1[which(deg_dry_dog1$padj< 0.05),])

# with main co-expressed gene groups of the time-course experiment
plot <- list(clusters_timepoint [[1]], clusters_timepoint [[2]], affected_dog1)
names(plot) <- c("time course cluster 1","time course cluster 2", "dry seeds affected")
plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF", "#868686FF"))

# with single seed DEGs of the dog1-4 experiment
plot <- list(rownames(deg_dog1$SD_Col0_3d_SD_dog1_3d), 
             rownames(deg_dog1$SD_Col0_7d24h_SD_dog1_7d24h), affected_dog1)
names(plot) <- c("dog1 3d","dog1 7d24h", "dry seeds affected")
plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF", "#868686FF"))
```
<img src="https://github.com/mk1859/single_seed/blob/main/images/clusters_dog1.png" width=20% height=20%> <img src="https://github.com/mk1859/single_seed/blob/main/images/ss_dog1_dog1.png" width=20% height=20%>

We identified GO terms enriched among genes with down- and upregulated expression.
``` R
# find enriched GO terms
up_genes <- rownames (deg_dry_dog1[which(deg_dry_dog1$padj< 0.05 & deg_dry_dog1$log2FoldChange > 0),])
dw_genes <- rownames (deg_dry_dog1[which(deg_dry_dog1$padj< 0.05 & deg_dry_dog1$log2FoldChange < 0),])
background <- rownames (deg_dry_dog1[which(deg_dry_dog1$baseMean > 1),])

up_genes <- gost(query = up_genes, organism = "athaliana", 
            custom_bg = background, user_threshold = 0.05,
            sources = "GO")$result

dw_genes <- gost(query = dw_genes, organism = "athaliana", 
            custom_bg = background, user_threshold = 0.05,
            sources = "GO")$result

# create data frame
affected <- list(up= up_genes, dw= dw_genes)
affected <- rbindlist(affected, idcol = "change")

# sort data frame
affected <- affected [order(affected$p_value, decreasing = TRUE),]
lev_order <- as.factor(affected$term_name [!duplicated(affected$term_name)])
affected$term_name <- factor(affected$term_name,levels = lev_order)

# plot for enrichments
g1 <- ggplot(affected, aes(1, term_name)) + 
         geom_tile(aes(fill = -log10(p_value))) + 
         scale_fill_gradientn(colors =c("white","darkred"), limits= c( 0, -log10(min(affected$p_value))))+
         theme_classic() +
         facet_grid(vars(change), scales = "free", space = "free_y") +
         theme(axis.line.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.x=element_blank()) 

# plot for types of ontologies
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
<img src="https://github.com/mk1859/single_seed/blob/main/images/degs_dog1_go.png" width=50% height=50%>

Finally, we created the signature of *dog1-4* affected genes using Vision and overlaid it on the time-course experiment PCA map.

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
<img src="https://github.com/mk1859/single_seed/blob/main/images/dog1_4_sign.png" width=50% height=50%>

```
session_info()

R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.1.1                GO.db_3.12.1                eulerr_6.1.1                RBGL_1.66.0                
 [5] scater_1.18.6               scran_1.18.7                SingleCellExperiment_1.12.0 patchwork_1.1.1            
 [9] cowplot_1.1.1               rlist_0.4.6.2               sctransform_0.3.2.9006      viridis_0.6.2              
[13] viridisLite_0.4.0           ggbeeswarm_0.6.0            UpSetR_1.4.0                data.table_1.14.2          
[17] org.At.tair.db_3.12.0       AnnotationDbi_1.52.0        rrvgo_1.2.0                 gprofiler2_0.2.1           
[21] ggthemes_4.2.4              forcats_0.5.1               stringr_1.4.0               dplyr_1.0.7                
[25] purrr_0.3.4                 readr_2.1.1                 tidyr_1.1.4                 tibble_3.1.6               
[29] ggplot2_3.3.5               tidyverse_1.3.1             SeuratObject_4.0.4          Seurat_4.0.6               
[33] outliers_0.14               onewaytests_2.6             car_3.0-12                  carData_3.0-5              
[37] Biostrings_2.58.0           XVector_0.30.0              DESeq2_1.30.1               SummarizedExperiment_1.20.0
[41] Biobase_2.50.0              MatrixGenerics_1.2.1        matrixStats_0.61.0          graph_1.68.0               
[45] biomaRt_2.46.3              VISION_2.1.0                GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
[49] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1        

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3            scattermore_0.7           bit64_4.0.5               knitr_1.37               
  [5] irlba_2.3.5               DelayedArray_0.16.3       rpart_4.1-15              RCurl_1.98-1.5           
  [9] generics_0.1.1            RSQLite_2.2.9             RANN_2.6.1                future_1.23.0            
 [13] wordspace_0.2-6           bit_4.0.4                 tzdb_0.2.0                spatstat.data_2.1-2      
 [17] xml2_1.3.3                lubridate_1.8.0           httpuv_1.6.5              assertthat_0.2.1         
 [21] xfun_0.29                 hms_1.1.1                 evaluate_0.14             promises_1.2.0.1         
 [25] fansi_1.0.0               progress_1.2.2            dbplyr_2.1.1              readxl_1.3.1             
 [29] igraph_1.2.11             DBI_1.1.2                 geneplotter_1.68.0        htmlwidgets_1.5.4        
 [33] sparsesvd_0.2             spatstat.geom_2.3-1       ellipsis_0.3.2            backports_1.4.1          
 [37] permute_0.9-5             annotate_1.68.0           gridBase_0.4-7            sparseMatrixStats_1.2.1  
 [41] moments_0.14              deldir_1.0-6              vctrs_0.3.8               ROCR_1.0-11              
 [45] abind_1.4-5               cachem_1.0.6              withr_2.4.3               treemap_2.4-3            
 [49] vegan_2.5-7               prettyunits_1.1.1         mclust_5.4.9              goftest_1.2-3            
 [53] cluster_2.1.0             lazyeval_0.2.2            crayon_1.4.2              webutils_1.1             
 [57] genefilter_1.72.1         edgeR_3.32.1              pkgconfig_2.0.3           slam_0.1-50              
 [61] labeling_0.4.2            vipor_0.4.5               nlme_3.1-148              wordcloud_2.6            
 [65] rlang_0.4.12              globals_0.14.0            lifecycle_1.0.1           miniUI_0.1.1.1           
 [69] loe_1.1                   BiocFileCache_1.14.0      modelr_0.1.8              rsvd_1.0.5               
 [73] cellranger_1.1.0          polyclip_1.10-0           lmtest_0.9-39             Matrix_1.4-0             
 [77] zoo_1.8-9                 beeswarm_0.4.0            reprex_2.0.1              ggridges_0.5.3           
 [81] pheatmap_1.0.12           png_0.1-7                 bitops_1.0-7              KernSmooth_2.23-17       
 [85] DelayedMatrixStats_1.12.3 blob_1.2.2                parallelly_1.30.0         iotools_0.3-2            
 [89] beachmat_2.6.4            memoise_2.0.1             magrittr_2.0.1            plyr_1.8.6               
 [93] ica_1.0-2                 zlibbioc_1.36.0           compiler_4.0.2            dqrng_0.3.0              
 [97] RColorBrewer_1.1-2        fitdistrplus_1.1-6        cli_3.1.0                 listenv_0.8.0            
[101] pbapply_1.5-0             MASS_7.3-51.6             mgcv_1.8-31               tidyselect_1.1.1         
[105] stringi_1.7.6             swagger_3.33.1            yaml_2.2.1                GOSemSim_2.16.1          
[109] BiocSingular_1.6.0        askpass_1.1               locfit_1.5-9.4            ggrepel_0.9.1            
[113] pbmcapply_1.5.0           grid_4.0.2                tools_4.0.2               future.apply_1.8.1       
[117] rstudioapi_0.13           bluster_1.0.0             logging_0.10-108          gridExtra_2.3            
[121] farver_2.1.0              Rtsne_0.15                digest_0.6.29             shiny_1.7.1              
[125] nortest_1.0-4             Rcpp_1.0.7                broom_0.7.11              scuttle_1.0.4            
[129] later_1.3.0               RcppAnnoy_0.0.19          httr_1.4.2                colorspace_2.0-2         
[133] rvest_1.0.2               XML_3.99-0.8              fs_1.5.2                  tensor_1.5               
[137] reticulate_1.22           splines_4.0.2             statmod_1.4.36            uwot_0.1.11              
[141] spatstat.utils_2.3-0      plotly_4.10.0             xtable_1.8-4              jsonlite_1.7.2           
[145] NLP_0.2-1                 plumber_1.1.0             R6_2.5.1                  tm_0.7-8                 
[149] pillar_1.6.4              htmltools_0.5.2           mime_0.12                 glue_1.6.0               
[153] fastmap_1.1.0             BiocParallel_1.24.1       BiocNeighbors_1.8.2       codetools_0.2-16         
[157] utf8_1.2.2                lattice_0.20-41           spatstat.sparse_2.1-0     curl_4.3.2               
[161] leiden_0.3.9              openssl_1.4.6             limma_3.46.0              survival_3.1-12          
[165] rmarkdown_2.11            fastICA_1.2-3             munsell_0.5.0             GenomeInfoDbData_1.2.4   
[169] haven_2.4.3               reshape2_1.4.4            gtable_0.3.0              spatstat.core_2.3-2  
```
