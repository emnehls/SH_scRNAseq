---
title: "SH_scRNAseq"
author: "Michael Nehls"
date: "08/19/2020"
output: 
  html_document:
    keep_md: true
---



## Importing the data and integrating the two Naive Data Sets

Using Seurat, data is imported from each 10X run. 
Then anchor genes are identified between the KO and WT conditions, and the data sets are integrated into 1. 
From here, PCA can be performed and UMAP plots can be generated. 
I imported the "filtered" matrices from the CellRanger output provided by the genome institute. 
It seems like they have already undergone some quality control. 
We also have the "raw" matrices and the initial FASTQ files. 
The packaged Mock treated files are not unpacking nicely (they are corrupted) on the cluster, 
so I still need to verify whether they are alright on our RAID system. 


```r
# loading Naive condition datasets for BMDMs from WT and KO genotypes
# setting them as Seurat Objects
wt_naive.data <- Read10X(data.dir = "WT_Naive/filtered_feature_bc_matrix/")
wt_naive <- CreateSeuratObject(
                               counts = wt_naive.data,
                               project = "WT_Naive",
                               min.cells = 3,
                               min.features = 200
)
wt_naive@meta.data[, "condition"] <- "wt_naive"
wt_naive
```

```
## An object of class Seurat 
## 16140 features across 6405 samples within 1 assay 
## Active assay: RNA (16140 features, 0 variable features)
```

```r
ko_naive.data <- Read10X(data.dir = "KO_Naive/filtered_feature_bc_matrix/")
ko_naive <- CreateSeuratObject(
                               counts = ko_naive.data,
                               project = "WT_Naive",
                               min.cells = 3,
                               min.features = 200
)

ko_naive@meta.data[, "condition"] <- "ko_naive"
ko_naive
```

```
## An object of class Seurat 
## 16730 features across 15116 samples within 1 assay 
## Active assay: RNA (16730 features, 0 variable features)
```

```r
#adding both naive objects into a list
naive.list <- c(wt_naive, ko_naive)

# before finding anchors between datasets, 
# do standard preprocessing and identify variable features
for (i in seq_along(naive.list)) {
  naive.list[[i]] <- NormalizeData(naive.list[[i]], verbose = FALSE)
  naive.list[[i]] <- FindVariableFeatures(
                                    naive.list[[i]],
                                    selection.method = "vst",
                                    nfeatures = 2000,
                                    verbose = FALSE
  )
}  # vst = variance stabilizing transformation

# identify anchors using FindIntegrationAnchors function,
# takes list of Seurat objects as input

naive.anchors <- FindIntegrationAnchors(object.list = naive.list, dims = 1:30)
```

```
## Warning in CheckDuplicateCellNames(object.list = object.list): Some cell names
## are duplicated across objects provided. Renaming to enforce unique cell names.
```

```r
naive.integrated <- IntegrateData(anchorset = naive.anchors, dims = 1:30)
```

```
## Warning: Adding a command log without an assay associated with it
```

```r
# switch to integrated assay. The variable features of this assay are 
# automatically set during IntegrateData

DefaultAssay(naive.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
naive.integrated <- ScaleData(naive.integrated, verbose = FALSE)
naive.integrated <- RunPCA(naive.integrated, npcs = 30, verbose = FALSE)
naive.integrated <- RunUMAP(naive.integrated, reduction = "pca", dims = 1:30)
```

```
## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## This message will be shown once per session
```

```r
p1 <- DimPlot(
              naive.integrated,
              reduction = "umap",
              group.by = "condition",
              label = TRUE,
              repel = TRUE
)
```

```
## Warning: Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
## Please use `as_label()` or `as_name()` instead.
## This warning is displayed once per session.
```

```r
p1
```

![](SH_scRNAseq_files/figure-html/Import-1.png)<!-- -->

## We can also identify clusters via Seurat!

### Pooled Cluster Plot


```r
naive.integrated <- FindNeighbors(
                                  naive.integrated,
                                  reduction = "pca",
                                  dims = 1:30
)
naive.integrated <- FindClusters(naive.integrated, resolution = 0.5)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 21521
## Number of edges: 787280
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9017
## Number of communities: 15
## Elapsed time: 5 seconds
```

```r
p2 <- DimPlot(naive.integrated, reduction = "umap", label = TRUE)
cowplot::plot_grid(p1, p2)
```

![](SH_scRNAseq_files/figure-html/Clustering_merged-1.png)<!-- -->

```r
p2
```

![](SH_scRNAseq_files/figure-html/Clustering_merged-2.png)<!-- -->

### Cluster by condition

```r
DimPlot(naive.integrated, reduction = "umap", split.by = "condition")
```

![](SH_scRNAseq_files/figure-html/Clustering_by_condition-1.png)<!-- -->

## Now we can plot the data and identify the gene expression in each of the clusters!

There are two primary ways to look at gene expression that I have found. 

### Looking per cluster at the highest expressed genes

The first is to ask what genes are most highly expressed in each cluster. 
This takes a while but gives a ton of data about which genes are determining cluster identity. 

I played around with setting cut-offs to sort the data similar to the Artyomov's shiny app. 

The last table generated by this code called "naive.markers" seems to give a summary of what genes
are most responsible for determining the identity of each cluster. 
The command that generates this can be fed parameters for only looking genes that are 
positively correlated with a cluster identity, defining a minimum percentage expression within the cluster,
and defining a minimum log-fold-change over other clusters to give more or less restrictive results.
My list was generated to only identify genes that are positively correlated with a cluster identity,
are expressed in a minimum of 25 percent of the cells in the cluster, and are at least expressed 0.25 
log-fold-change higher within the cluster. 


```r
# Try to identify the cell type identities of the clusters
DefaultAssay(naive.integrated) <- "RNA"

# A function to get a dataframe of markers if one doesn't exist, or renew it if it does
# Also saves dataframe as both an rds and csv
renew_markers <- function(cluster_number) {
  markers <- FindConservedMarkers(
    naive.integrated,
    ident.1 = cluster_number,
    grouping.var = "condition",
    verbose = FALSE
  )
  write_rds(markers, str_c("./files/markers_rds/c", as.character(cluster_number), "markers.rds")) 
  write.csv(markers, 
    file = str_c("./files/csv_xlsx_files/c", as.character(cluster_number), "markers.csv"), 
    append = FALSE, row.names = TRUE)
  return(markers)
}

# Checks for an rds or renews the markers frame - set force = TRUE to renew all markers frames
get_markers <- function(cluster_number, force = FALSE) {
  markers_name = str_c("c", as.character(cluster_number), "markers")
  markers_file = str_c("./files/markers_rds/", markers_name, ".rds")
  if (!force & file.exists(markers_file)){
    markers <- read_rds(markers_file)
    return(markers)
  } else {
    markers <- renew_markers(cluster_number)
    return(markers)
  }
}


for (cluster in seq_along(levels(naive.integrated))) {
  markers <- get_markers(cluster - 1) 
  assign(str_c("c", as.character(cluster - 1), "markers"), markers, envir = .GlobalEnv)
}

c0markers_filtered <- c0markers[c0markers$wt_naive_pct.1 > 0.7 & c0markers$wt_naive_pct.2 < 0.7,]
head(c0markers_filtered)
```

```
##        wt_naive_p_val wt_naive_avg_logFC wt_naive_pct.1 wt_naive_pct.2
## Clec4d  1.720232e-213          0.9368841          0.870          0.519
## Lst1    1.995378e-151          0.7619648          0.851          0.645
## Dock10   1.920961e-29          0.4435840          0.746          0.680
## Cd68    1.148800e-112          0.6311869          0.866          0.674
## Card19  1.069034e-119          0.5875151          0.842          0.693
## F10      8.689163e-44          0.6572991          0.710          0.581
##        wt_naive_p_val_adj ko_naive_p_val ko_naive_avg_logFC ko_naive_pct.1
## Clec4d      2.957767e-209   0.000000e+00          0.6969278          0.892
## Lst1        3.430853e-147  1.115581e-277          0.5586259          0.825
## Dock10       3.302901e-25  5.843079e-257          0.7888004          0.717
## Cd68        1.975247e-108  3.767497e-248          0.5916928          0.803
## Card19      1.838096e-115  2.418859e-239          0.5727743          0.763
## F10          1.494015e-39  8.780598e-235          0.7369743          0.761
##        ko_naive_pct.2 ko_naive_p_val_adj      max_pval minimump_p_val
## Clec4d          0.733       0.000000e+00 1.720232e-213   0.000000e+00
## Lst1            0.715      1.918130e-273 1.995378e-151  2.231162e-277
## Dock10          0.595      1.004659e-252  1.920961e-29  1.168616e-256
## Cd68            0.667      6.477834e-244 1.148800e-112  7.534994e-248
## Card19          0.692      4.158986e-235 1.069034e-119  4.837717e-239
## F10             0.658      1.509736e-230  8.689163e-44  1.756120e-234
```

```r
c1markers_filtered <- c1markers[c1markers$wt_naive_pct.1 > 0.7 & c1markers$wt_naive_pct.2 < 0.6, ]
head(c1markers_filtered)
```

```
##         wt_naive_p_val wt_naive_avg_logFC wt_naive_pct.1 wt_naive_pct.2
## Csf3r     0.000000e+00           1.855398          0.824          0.166
## Cxcr2     0.000000e+00           1.632500          0.747          0.082
## Mmp9      0.000000e+00           1.605736          0.872          0.183
## Sell      0.000000e+00           1.472886          0.788          0.143
## Pglyrp1  2.611227e-208           1.118199          0.786          0.247
## Slfn1    2.557595e-204           1.296735          0.729          0.234
##         wt_naive_p_val_adj ko_naive_p_val ko_naive_avg_logFC ko_naive_pct.1
## Csf3r         0.000000e+00              0          1.8354566          0.890
## Cxcr2         0.000000e+00              0          1.6189756          0.688
## Mmp9          0.000000e+00              0          1.2719350          0.665
## Sell          0.000000e+00              0          1.3040014          0.615
## Pglyrp1      4.489745e-204              0          0.6234252          0.670
## Slfn1        4.397528e-200              0          1.2815761          0.677
##         ko_naive_pct.2 ko_naive_p_val_adj      max_pval minimump_p_val
## Csf3r            0.343                  0  0.000000e+00              0
## Cxcr2            0.132                  0  0.000000e+00              0
## Mmp9             0.249                  0  0.000000e+00              0
## Sell             0.186                  0  0.000000e+00              0
## Pglyrp1          0.295                  0 2.611227e-208              0
## Slfn1            0.287                  0 2.557595e-204              0
```

```r
c2markers_filtered <- c2markers[c2markers$wt_naive_pct.1 > 0.7 & c2markers$wt_naive_pct.2 < 0.6, ]
head(c2markers_filtered)
```

```
##        wt_naive_p_val wt_naive_avg_logFC wt_naive_pct.1 wt_naive_pct.2
## Cd36     0.000000e+00          1.5476255          0.906          0.232
## Mmp19    0.000000e+00          1.4303439          0.936          0.258
## Anpep    0.000000e+00          0.8724305          0.917          0.238
## Itgax   9.074890e-281          1.2192506          0.981          0.435
## Lpl     1.244749e-278          2.1184529          1.000          0.516
## Rnf128  1.158122e-275          0.8550737          0.757          0.146
##        wt_naive_p_val_adj ko_naive_p_val ko_naive_avg_logFC ko_naive_pct.1
## Cd36         0.000000e+00              0          1.6991462          0.915
## Mmp19        0.000000e+00              0          1.6348371          0.952
## Anpep        0.000000e+00              0          0.8124444          0.905
## Itgax       1.560337e-276              0          1.1716830          0.959
## Lpl         2.140222e-274              0          1.8982026          0.998
## Rnf128      1.991275e-271              0          0.7885583          0.828
##        ko_naive_pct.2 ko_naive_p_val_adj      max_pval minimump_p_val
## Cd36            0.172                  0  0.000000e+00              0
## Mmp19           0.260                  0  0.000000e+00              0
## Anpep           0.187                  0  0.000000e+00              0
## Itgax           0.266                  0 9.074890e-281              0
## Lpl             0.512                  0 1.244749e-278              0
## Rnf128          0.180                  0 1.158122e-275              0
```

```r
## Compare clusters to identify the most significant genes in determining their identity
renew_naive_markers = FALSE
if (!file.exists("./files/markers_rds/naive_markers.rds") | renew_naive_markers == TRUE) {
  naive.markers <- FindAllMarkers(naive.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write_rds(naive.markers, "./files/markers_rds/naive_markers.rds")
  write.csv(naive.markers, "./files/csv_xlsx_files/naive_markers.csv", row.names = TRUE, append = FALSE)
} else {
  naive.markers <- read_rds("./files/markers_rds/naive_markers.rds")
}
naive.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC) # top 4 cluster defining genes for each cluster
```

```
## # A tibble: 60 x 7
## # Groups:   cluster [15]
##        p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene   
##        <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
##  1 0.            0.950 0.799 0.619 0.        0       Thbs1  
##  2 0.            0.808 0.882 0.788 0.        0       Creg1  
##  3 2.58e-261     1.12  0.273 0.095 4.44e-257 0       Ear6   
##  4 1.04e-256     0.821 0.636 0.486 1.79e-252 0       Slc7a11
##  5 0.            1.89  0.878 0.287 0.        1       Csf3r  
##  6 0.            1.65  0.699 0.116 0.        1       Cxcr2  
##  7 0.            1.61  0.836 0.365 0.        1       Il1r2  
##  8 1.35e-298     1.45  0.413 0.167 2.32e-294 1       Ly6a   
##  9 0.            2.38  1     0.79  0.        2       Mmp12  
## 10 0.            1.97  0.998 0.513 0.        2       Lpl    
## # … with 50 more rows
```

## Comparing relative expression between closely related clusters

A similar method can be used to compare two clusters directly.
This can be helpful for determining identities of closely related clusters.


```r
# a function to generate the comparison table and save it to a csv 
comparison_markers <- function(cluster1, cluster2) {
  markers <- FindMarkers(naive.integrated, ident.1 = cluster1, ident.2 = cluster2, min.pct = 0.25) 
  write_rds(markers, str_c("./files/markers_rds/c", as.character(cluster1), "v", as.character(cluster2), "markers.rds")) 
  write.csv(markers, 
    file = str_c("./files/csv_xlsx_files/c", as.character(cluster1), "v", as.character(cluster2), "markers.csv"), 
    append = FALSE, row.names = TRUE)
  return(markers)
}

# Checks for an rds or renews the markers frame - set force = TRUE to renew all markers frames
get_comparison_markers <- function(cluster1, cluster2, force = FALSE) {
  markers_name = str_c("c", as.character(cluster1), "v", as.character(cluster2), "markers")
  markers_file = str_c("./files/markers_rds/", markers_name, ".rds")
  if (!force & file.exists(markers_file)){
    markers <- read_rds(markers_file)
    return(markers)
  } else {
    markers <- comparison_markers(cluster1, cluster2)
    return(markers)
  }
}

c2v3markers <- get_comparison_markers(2, 3)
c0v1markers <- get_comparison_markers(0, 1)
c5v8markers <- get_comparison_markers(5, 8)
c1v6markers <- get_comparison_markers(1, 6)
```


### Looking at expression of predetermined genes

The second is by asking what the expression profile is of makers that we know are associated with different cell types.


```r
## Plot gene expression by feature - way faster!
p_macs <- FeaturePlot(naive.integrated, features = c("Cd14", "Cd68", "Cd62l", "Cd31", "Cx3cr1", "Ace", "Ccr2", "Lyz2", "Lpl"), min.cutoff = "q9")
```

```
## Warning in FetchData(object = object, vars = c(dims, "ident", features), : The
## following requested variables were not found: Cd62l, Cd31
```

```r
p_macs
```

![](SH_scRNAseq_files/figure-html/featurePlot-1.png)<!-- -->

```r
p_PMN <- FeaturePlot(naive.integrated, features = c("Clec4d", "S100a8", "S100a9", "Acod1"), min.cutoff = "q9")
p_PMN
```

![](SH_scRNAseq_files/figure-html/featurePlot-2.png)<!-- -->

```r
p_DC <- FeaturePlot(naive.integrated, features = c("Ctsh", "Ptms", "Napsa", "Ccr3", "Itgae","Flt3","Cxcl16", "Syngr2"), min.cutoff = "q9")
p_DC
```

![](SH_scRNAseq_files/figure-html/featurePlot-3.png)<!-- -->

## Differential expression analysis

Now it is time for differential expression analysis! I am thinking about approaching it in two ways: 

* A differential expression gene table

* A volcano style-plot

### Generating the differential expression tables


```r
naive.integrated$celltype.condition <- paste(Idents(naive.integrated), naive.integrated$condition, sep = "_")
naive.integrated$celltype <- Idents(naive.integrated)
Idents(naive.integrated) <- "celltype.condition" # changes levels in naive integrated to one per cluster per condition

# A function to generate the differential expression table, assign it to a global variable, and save a csv
cluster_differential <- function(cluster) {
  differential_name = str_c("c", as.character(cluster), "_differential")
  differential_file_name = str_c("./files/csv_xlsx_files/", differential_name, ".csv")
  differential <- FindMarkers(naive.integrated, ident.1 = str_c(as.character(cluster), "_wt_naive"), 
                               ident.2 = str_c(as.character(cluster), "_ko_naive"), verbose = FALSE)
  write.csv(differential, file = differential_file_name, append = FALSE, row.names = TRUE)
  assign(differential_name, differential, envir = .GlobalEnv)
}

for (cluster_plus_one in seq_along(levels(naive.integrated$seurat_clusters))) {
  cluster_differential(cluster_plus_one - 1)
}
```

```
## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored

## Warning in write.csv(differential, file = differential_file_name, append =
## FALSE, : attempt to set 'append' ignored
```

### Generating volcano-style plots


```r
# A function to plot the differential expression in a volcano style plot with genes labeled
plot_cluster <- function(cluster, genes = c("Bhlhe40")) {
    Idents(naive.integrated) <- "celltype"
    cluster_subset <- subset(naive.integrated, idents = as.character(cluster))
    Idents(cluster_subset) <- "condition"
    avg_exp <- log1p(AverageExpression(cluster_subset, verbose = FALSE)$RNA)
    avg_exp$gene <- rownames(avg_exp)

    cluster_plot <- ggplot(avg_exp, aes(wt_naive, ko_naive)) +
        geom_point() +
        ggtitle(paste("cluster", as.character(cluster), sep = " "))
    cluster_plot <- LabelPoints(
                                cluster_plot,
                                points = genes,
                                repel = TRUE
    )
    cluster_plot
}
```

### Cluster 0


```r
# Look at the beginning of the differential expression file to choose genes to label
head(c0_differential, 10)
```

```
##                 p_val avg_logFC pct.1 pct.2     p_val_adj
## Bhlhe40 6.954264e-160 0.6015306 0.733 0.251 1.195716e-155
## Rps3a1  4.709849e-137 0.4873808 0.984 0.941 8.098114e-133
## Rps7    1.049153e-132 0.6306760 0.961 0.829 1.803914e-128
## Rpl9    8.982399e-132 0.5531451 0.954 0.848 1.544434e-127
## Rpl5    1.290399e-124 0.5568293 0.955 0.810 2.218712e-120
## Rpl23   2.399503e-119 0.4958967 0.985 0.947 4.125706e-115
## Rps24   4.191520e-119 0.4693191 0.986 0.938 7.206900e-115
## Rpl17   4.120666e-118 0.5263964 0.978 0.908 7.085073e-114
## Rps3    2.246752e-117 0.5172118 0.961 0.869 3.863065e-113
## Abcd2   4.891208e-116 0.2592643 0.308 0.051 8.409943e-112
```

```r
genes_to_label <- c("Bhlhe40", "Rps3a1", "Rps7", "Rp19")

p_c0 <- plot_cluster(0, genes_to_label)
plot_grid(p2, p_c0) # p2 is included to give reference to which cluster we are looking into
```

![](SH_scRNAseq_files/figure-html/cluster0_diffex_plot-1.png)<!-- -->


### Cluster 1


```r
head(c1_differential, 10)
```

```
##               p_val  avg_logFC pct.1 pct.2    p_val_adj
## Apoe   7.749334e-86  1.6736564 0.266 0.028 1.332421e-81
## Mmp12  1.216946e-78 -1.1213501 0.311 0.661 2.092417e-74
## Dhrs3  1.745614e-59  0.3195280 0.216 0.028 3.001408e-55
## Rpl9   3.778507e-59  0.5033640 0.943 0.764 6.496765e-55
## Scd2   4.687399e-51  0.5432857 0.430 0.140 8.059515e-47
## Adgre5 6.037115e-50  0.3346634 0.198 0.029 1.038022e-45
## Rps11  5.778174e-48  0.4824353 0.974 0.816 9.934992e-44
## Trf    2.671131e-46  0.6212355 0.518 0.208 4.592743e-42
## Rps27a 8.177961e-46  0.4159935 0.985 0.909 1.406119e-41
## Fau    1.492318e-45  0.3764314 1.000 0.987 2.565891e-41
```

```r
genes_to_label <- c("Apoe", "Mmp12", "Dhrs3", "Rp19")
p_c1 <- plot_cluster(1, genes_to_label)
plot_grid(p2, p_c1)
```

![](SH_scRNAseq_files/figure-html/cluster1_diffexpression-1.png)<!-- -->


### Cluster 2


```r
head(c2_differential, n = 10)
```

```
##                  p_val  avg_logFC pct.1 pct.2     p_val_adj
## Bhlhe40  3.378246e-183  0.7844712 0.934 0.493 5.808555e-179
## AA467197 3.913614e-169  1.5145675 0.992 0.770 6.729068e-165
## Ccl2     2.648943e-140 -1.6007202 0.551 0.955 4.554593e-136
## Vim      1.961268e-132 -0.6054507 0.996 1.000 3.372203e-128
## S100a9   2.043107e-115 -0.9325852 0.705 0.917 3.512918e-111
## Klk1b1   3.976240e-112  0.5130597 0.514 0.099 6.836747e-108
## Klk1b11  1.413795e-111  0.8233073 0.580 0.138 2.430880e-107
## S100a6   1.465303e-110 -0.6657772 0.983 0.997 2.519442e-106
## Ccl7     1.149144e-108 -1.4666799 0.135 0.683 1.975838e-104
## Rgl1      2.755901e-98 -0.6673743 0.590 0.881  4.738495e-94
```

```r
genes_to_label <- c("Bhlhe40", "AA467197", "Ccl2", "Vim", "S100a9")

p_c2 <- plot_cluster(2, genes_to_label)
plot_grid(p2, p_c2)
```

![](SH_scRNAseq_files/figure-html/cluster2_diffexpression-1.png)<!-- -->


### Cluster 3


```r
head(c3_differential, n = 10)
```

```
##                p_val  avg_logFC pct.1 pct.2     p_val_adj
## Ccl2   8.571257e-229 -1.9495986 0.216 0.874 1.473742e-224
## Rps9   6.921446e-165  0.5782080 0.999 0.999 1.190073e-160
## Rps27a 2.995029e-160  0.4296387 1.000 0.999 5.149653e-156
## Rbpms  8.698002e-152  0.4277515 0.712 0.178 1.495534e-147
## Rpl5   1.703785e-150  0.5204401 0.995 0.985 2.929488e-146
## Rps3a1 4.921082e-150  0.4268346 1.000 0.999 8.461309e-146
## Rps14  2.115263e-149  0.5358750 1.000 0.997 3.636984e-145
## Axl    2.787599e-149  0.7754919 0.904 0.477 4.792998e-145
## Rpl17  1.160402e-147  0.5329156 1.000 0.995 1.995196e-143
## Emp1   1.018207e-145 -0.7138134 0.317 0.770 1.750706e-141
```

```r
genes_to_label <- c("Ccl2", "Rps9", "Rps27a", "Rbpms", "Rpl5")

p_c3 <- plot_cluster(3, genes_to_label)
plot_grid(p2, p_c3)
```

![](SH_scRNAseq_files/figure-html/cluster3_diffexpression-1.png)<!-- -->


### Cluster 4


```r
head(c4_differential, n = 10)
```

```
##                 p_val  avg_logFC pct.1 pct.2    p_val_adj
## Ccl2     7.483327e-62 -1.5915345 0.115 0.521 1.286683e-57
## AA467197 1.983385e-52  1.1866553 0.654 0.292 3.410232e-48
## Mmp12    5.017082e-47 -1.3076659 0.574 0.842 8.626370e-43
## Rpl32    8.921177e-47  0.6296226 0.930 0.808 1.533907e-42
## Rps7     8.265880e-45  0.8129538 0.833 0.578 1.421235e-40
## Rps14    1.107081e-44  0.6775904 0.900 0.731 1.903515e-40
## Tmsb4x   1.132870e-42  0.3941733 0.998 0.995 1.947857e-38
## Rps3a1   1.786508e-41  0.5688486 0.917 0.807 3.071722e-37
## Rpl23    9.184313e-41  0.5207116 0.954 0.858 1.579151e-36
## Rpl17    1.159263e-40  0.6681364 0.879 0.704 1.993236e-36
```

```r
genes_to_label <- c("Ccl2", "AA467197", "Mmp12", "Rpl32", "Rps7")

p_c4 <- plot_cluster(4, genes_to_label)
plot_grid(p2, p_c4)
```

![](SH_scRNAseq_files/figure-html/cluster4_diffexpression-1.png)<!-- -->

### Cluster 5


```r
head(c5_differential, n = 10)
```

```
##                 p_val  avg_logFC pct.1 pct.2     p_val_adj
## Hpgds   1.527173e-164 -0.9766181 0.215 0.887 2.625822e-160
## Lmna    2.764033e-136 -0.9535482 0.460 0.897 4.752478e-132
## Bhlhe40 3.315825e-127  0.7098046 0.938 0.534 5.701229e-123
## Angptl2 1.042759e-117 -0.6919939 0.122 0.688 1.792920e-113
## Crip1   4.611598e-102  1.0784219 0.955 0.706  7.929181e-98
## Klk1b11  2.297853e-95  1.3769576 0.884 0.407  3.950928e-91
## Cdkn1a   1.078654e-94 -0.9262530 0.810 0.977  1.854638e-90
## Mrfap1   2.589254e-92  0.6783343 0.977 0.918  4.451963e-88
## Sdc4     5.597997e-90 -0.7301269 0.162 0.655  9.625196e-86
## Trf      1.978515e-89  0.6874814 0.935 0.650  3.401860e-85
```

```r
genes_to_label <- c("Hpgds", "Lmna", "Bhlhe40", "Angptl2", "Crip1")

p_c5 <- plot_cluster(5)
plot_grid(p2, p_c5)
```

![](SH_scRNAseq_files/figure-html/cluster5_diffexpression-1.png)<!-- -->

### Cluster 6


```r
head(c6_differential, n = 10)
```

```
##                 p_val  avg_logFC pct.1 pct.2    p_val_adj
## Mmp12    9.548411e-43 -0.9049781 0.335 0.727 1.641754e-38
## Apoe     2.711194e-31  1.1543665 0.199 0.019 4.661627e-27
## Lyz2     5.814871e-31  0.5781136 1.000 0.987 9.998089e-27
## Prdx5    8.847032e-29  0.3836727 0.993 0.964 1.521159e-24
## Nov      4.758120e-27 -0.7854010 0.382 0.719 8.181111e-23
## Arhgap24 4.961666e-22  0.2992745 0.257 0.058 8.531089e-18
## Klhl24   1.287542e-21  0.3388515 0.522 0.206 2.213800e-17
## Ypel3    2.552217e-20  0.4257260 0.868 0.605 4.388282e-16
## Txn1     4.984508e-20  0.4211647 0.963 0.870 8.570362e-16
## Rpl18a   6.087111e-20  0.4071151 0.982 0.858 1.046618e-15
```

```r
genes_to_label <- c("Mmp12", "Apoe", "Lyz2", "Prdx5")

p_c6 <- plot_cluster(6, genes_to_label)
plot_grid(p2, p_c6)
```

![](SH_scRNAseq_files/figure-html/cluster6_diffexpression-1.png)<!-- -->


### Cluster 7


```r
head(c7_differential, n = 10)
```

```
##                  p_val  avg_logFC pct.1 pct.2     p_val_adj
## Bhlhe40  4.543844e-114  0.7879753 0.950 0.488 7.812686e-110
## S100a9    1.175195e-83 -1.7917023 0.728 0.936  2.020630e-79
## Ccl2      1.690272e-80 -1.4191841 0.380 0.872  2.906254e-76
## Clec4d    9.277244e-75 -0.7337378 0.845 0.952  1.595129e-70
## Hpgds     6.369496e-71 -0.5197185 0.358 0.760  1.095171e-66
## AA467197  1.423884e-70  1.2012469 0.945 0.699  2.448226e-66
## Axl       5.339681e-66  0.6182239 0.915 0.557  9.181048e-62
## Rbpms     1.418615e-61  0.2594284 0.682 0.190  2.439166e-57
## Ccl7      5.544445e-55 -1.0259297 0.065 0.514  9.533119e-51
## Tmsb4x    7.100147e-52  0.2967838 1.000 1.000  1.220799e-47
```

```r
genes_to_label <- c("Bhlhe40", "S100a9", "Ccl2", "Clec4d", "Hpgds")

p_c7 <- plot_cluster(7, genes_to_label)
plot_grid(p2, p_c7)
```

![](SH_scRNAseq_files/figure-html/cluster7_diffexpression-1.png)<!-- -->

### Cluster 8


```r
head(c8_differential, n = 10)
```

```
##                  p_val  avg_logFC pct.1 pct.2    p_val_adj
## Sphk1    1.486203e-101 -1.0685672 0.113 0.775 2.555378e-97
## AA467197  3.177189e-72  1.2672922 0.961 0.717 5.462859e-68
## Bhlhe40   6.090276e-69  0.7472052 0.944 0.772 1.047162e-64
## Mmp12     3.736333e-68 -1.3320456 0.559 0.912 6.424252e-64
## Satb1     1.525196e-65 -0.9451787 0.092 0.596 2.622422e-61
## Ninj1     1.153600e-59 -0.7323671 0.619 0.863 1.983500e-55
## Ndufa4    4.931661e-58 -0.6822237 0.711 0.866 8.479497e-54
## Psap      1.232995e-50 -0.8764742 0.965 0.988 2.120012e-46
## Mrfap1    1.539631e-48  0.7260068 0.966 0.891 2.647242e-44
## Malt1     1.156163e-47 -0.8713482 0.713 0.897 1.987906e-43
```

```r
genes_to_label <- c("Sphk1", "AA467197", "Bhlhe40", "Mmp12")

p_c8 <- plot_cluster(8, genes_to_label)
plot_grid(p2, p_c8)
```

![](SH_scRNAseq_files/figure-html/cluster8_diffexpression-1.png)<!-- -->


### Cluster 9


```r
head(c9_differential, n = 10)
```

```
##                 p_val  avg_logFC pct.1 pct.2    p_val_adj
## Rbpms    1.039109e-41  0.4116352 0.654 0.122 1.786645e-37
## Nov      3.980935e-39 -0.9786411 0.605 0.939 6.844820e-35
## S100a6   8.740415e-39 -0.7709357 0.962 1.000 1.502827e-34
## Heatr1   1.219463e-37 -0.6219832 0.411 0.819 2.096745e-33
## Trf      6.932017e-37  0.9257301 0.919 0.475 1.191891e-32
## Emilin2  4.106133e-30 -0.5362387 0.924 0.991 7.060085e-26
## Bhlhe40  1.070350e-29  0.5675734 0.778 0.383 1.840360e-25
## Ccl2     9.492188e-29 -1.3249827 0.108 0.564 1.632087e-24
## AW112010 8.641940e-28  0.5953417 0.676 0.229 1.485895e-23
## Rps7     3.783464e-27  0.3761633 1.000 1.000 6.505287e-23
```

```r
genes_to_label <- c("Rbpms", "Nov", "S100a6", "Heatr1", "Trf")

p_c9 <- plot_cluster(9, genes_to_label)
plot_grid(p2, p_c9)
```

![](SH_scRNAseq_files/figure-html/cluster9_diffexpression-1.png)<!-- -->


### Cluster 10


```r
head(c10_differential, n = 10)
```

```
##                p_val  avg_logFC pct.1 pct.2    p_val_adj
## Fth1    3.912476e-29 -0.6144991 1.000 1.000 6.727110e-25
## Ftl1    1.041726e-23 -0.6501416 0.973 0.990 1.791143e-19
## Lgals3  7.679980e-22 -0.7447641 0.858 0.977 1.320496e-17
## Bhlhe40 4.253901e-14  0.5250242 0.236 0.029 7.314157e-10
## Vim     6.097447e-14 -0.7692216 0.568 0.757 1.048395e-09
## S100a6  1.430993e-13 -0.6993532 0.520 0.721 2.460450e-09
## Ccl6    3.958499e-13 -0.3397800 0.818 0.958 6.806244e-09
## Prdx1   1.554256e-12 -0.7302395 0.608 0.799 2.672388e-08
## Wfdc17  3.565152e-12 -0.5736814 0.750 0.888 6.129922e-08
## Mmp12   3.723497e-12 -0.5782741 0.764 0.924 6.402181e-08
```

```r
genes_to_label <- c("Fth1", "Ftl1", "Lgals3", "Bhlhe40", "Vim")

p_c10 <- plot_cluster(10, genes_to_label)
plot_grid(p2, p_c10)
```

![](SH_scRNAseq_files/figure-html/cluster10_diffexpression-1.png)<!-- -->


### Cluster 11


```r
head(c11_differential, n = 10)
```

```
##                 p_val  avg_logFC pct.1 pct.2    p_val_adj
## Hpgds    1.337129e-49 -0.9384844 0.309 0.857 2.299060e-45
## Bhlhe40  3.449034e-47  0.7529492 0.964 0.578 5.930269e-43
## Lmna     1.372641e-38 -0.8159697 0.558 0.912 2.360120e-34
## AA467197 2.748202e-32  1.0436442 0.855 0.395 4.725259e-28
## Psap     1.249782e-30 -0.6743963 1.000 1.000 2.148875e-26
## Cldnd1   2.989231e-30  0.8645028 0.950 0.701 5.139683e-26
## Sdc4     4.179671e-30 -0.3521851 0.053 0.476 7.186527e-26
## Angptl2  5.609954e-29 -0.5241424 0.338 0.776 9.645755e-25
## Crip1    2.318548e-27  1.0071183 0.944 0.823 3.986511e-23
## Scp2     5.537610e-27 -0.4148359 0.985 1.000 9.521367e-23
```

```r
genes_to_label <- c("Hpgds", "Bhlhe40", "Lmna", "AA467197", "Psap")

p_c11 <- plot_cluster(11, genes_to_label)
plot_grid(p2, p_c11)
```

![](SH_scRNAseq_files/figure-html/cluster11_diffexpression-1.png)<!-- -->


### Cluster 12


```r
head(c12_differential, n = 10)
```

```
##                p_val  avg_logFC pct.1 pct.2    p_val_adj
## Mmp12   5.046758e-20 -0.7383272 0.614 0.910 8.677396e-16
## Gpnmb   6.247396e-09 -0.4467025 0.357 0.680 1.074177e-04
## Klhl24  4.128074e-08  0.2605828 0.586 0.242 7.097811e-04
## Adgre5  1.772135e-07  0.2750302 0.629 0.288 3.047010e-03
## Nov     2.289196e-07 -0.4576389 0.386 0.649 3.936044e-03
## Rpl41   7.166714e-07 -0.3334319 0.986 1.000 1.232245e-02
## Bhlhe40 1.259625e-06  0.2581269 0.557 0.271 2.165798e-02
## Atp5k   5.077538e-06 -0.2699531 0.900 0.939 8.730319e-02
## Pim1    6.375843e-06 -0.3138310 0.871 0.939 1.096262e-01
## Ifitm1  8.082751e-06 -0.3815488 0.829 0.949 1.389748e-01
```

```r
genes_to_label <- c("Mmp12", "Gpnmb", "Klhl24", "Adgre5", "Nov")

p_c12 <- plot_cluster(12, genes_to_label)
plot_grid(p2, p_c12)
```

![](SH_scRNAseq_files/figure-html/cluster12_diffexpression-1.png)<!-- -->


### Cluster 13


```r
head(c13_differential, n = 10)
```

```
##                p_val  avg_logFC pct.1 pct.2    p_val_adj
## Bhlhe40 2.855784e-16  1.9134470  0.88 0.165 4.910236e-12
## Elp5    1.054790e-08  0.5640104  0.40 0.031 1.813605e-04
## Rab10   2.130151e-07  0.6557639  0.80 0.268 3.662582e-03
## Mocos   4.313562e-07  0.2922728  0.28 0.016 7.416738e-03
## Mmp12   1.925518e-06 -1.1021425  0.44 0.819 3.310736e-02
## Ilkap   2.638766e-06  0.3046454  0.52 0.102 4.537093e-02
## Azi2    2.739080e-06  0.3561392  0.44 0.079 4.709573e-02
## Tnks2   3.057715e-06  0.4700137  0.56 0.134 5.257435e-02
## Glg1    3.149959e-06  0.4953702  0.52 0.118 5.416040e-02
## Snx30   4.827033e-06  0.5225079  0.36 0.055 8.299601e-02
```

```r
genes_to_label <- c("Bhlhe40", "Elp5", "Rab10", "Mocos", "Mmp12")

p_c13 <- plot_cluster(13, genes_to_label)
plot_grid(p2, p_c13)
```

![](SH_scRNAseq_files/figure-html/cluster13_diffexpression-1.png)<!-- -->


### Cluster 14


```r
head(c14_differential, n = 10)
```

```
##                p_val  avg_logFC pct.1 pct.2 p_val_adj
## Rundc1  1.255170e-05  0.4474208 0.357 0.016 0.2158140
## Cebpb   2.540184e-05 -1.0818905 0.429 0.906 0.4367592
## Lsm11   3.102335e-05  0.4186176 0.429 0.047 0.5334155
## Kdelr1  1.178642e-04  1.0207980 0.786 0.297 1.0000000
## Lgals1  1.552543e-04 -1.0264985 0.857 0.984 1.0000000
## Parp8   1.900117e-04  0.3745440 0.214 0.000 1.0000000
## Cdyl2   2.195729e-04  0.8734852 0.286 0.016 1.0000000
## Tgfb1   2.929054e-04 -0.8902553 0.571 0.891 1.0000000
## S100a6  3.655606e-04 -1.3189652 0.571 0.953 1.0000000
## Abhd17b 5.205270e-04  0.2588214 0.571 0.141 1.0000000
```

```r
genes_to_label <- c("Rundc1", "Cebpb", "Lsm11", "Kdelr1", "Lgals1")

p_14 <- plot_cluster(14, genes_to_label)
plot_grid(p2, p_14)
```

![](SH_scRNAseq_files/figure-html/cluster14_diffexpression-1.png)<!-- -->

## An attemp at pseudotime analysis

This requires the SeuratWrappers library and Monocle3.
It is pretty easy to generate!
I need to read/experiment more to know how to modify the way that Monocle3 clusters and partitions.
For example, I want to try to partition by what we thing the cells are.
For now, I have used the cluster with maximal c-kit expression as the root (time 0) of the analysis.


```r
cds <- as.cell_data_set(naive.integrated)
cds <- cluster_cells(cds)
p_monocle_clusters <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("Monocle3's Clusters")
```

```
## Warning: The `add` argument of `group_by()` is deprecated as of dplyr 1.0.0.
## Please use the `.add` argument instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_warnings()` to see where this warning was generated.
```

```r
p_monocle_partitions <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + ggtitle("Monocle3's Partitions")
wrap_plots(p_monocle_clusters, p_monocle_partitions) 
```

![](SH_scRNAseq_files/figure-html/pseudotime_experimentation-1.png)<!-- -->

```r
integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
```

```
##   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
```

```r
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
```

```
## Warning: `select_()` is deprecated as of dplyr 0.7.0.
## Please use `select()` instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_warnings()` to see where this warning was generated.
```

![](SH_scRNAseq_files/figure-html/pseudotime_experimentation-2.png)<!-- -->

```r
# Set max Kit expression as the root of our pseudotime tree
max.kit <- which.max(unlist(FetchData(integrated.sub, "Kit")))
max.kit <- colnames(integrated.sub)[max.kit]
cds <- order_cells(cds, root_cells = max.kit)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = TRUE) +
            ggtitle("Pseudotime map with max Kit expression as root")
```

![](SH_scRNAseq_files/figure-html/pseudotime_experimentation-3.png)<!-- -->

```r
# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = "integrated")
p_pseudotime <- FeaturePlot(integrated.sub, "monocle3_pseudotime")
wrap_plots(p2, p_pseudotime)
```

![](SH_scRNAseq_files/figure-html/pseudotime_experimentation-4.png)<!-- -->

Look at that! This is pretty exciting! I don't have any breakthrough interpretations, but it looks like clusters 2 and 5 have the highest expression of c-kit. It looks like cluster 6 and 13 are very mature, which fits with our interpretation of them as mature PMNs. That is all I really can say, given how early on I am in understanding this analysis. But just for fun, I will overinterpret it! It looks like there may be two distinct paths forward in pseudotime. The one to the left would fit with the interpretation that those clusters are DC like. The one to the right is fuzzier, expecially because given cell-types don't cluster as nicely. It looks like cluster 12 may be a less mature PMN cluster. And maybe cluster 9 is a mature macrophage like cluster (expresses Ccr2 highly)? It is a bit weird that clusters 0 and 1 have multiple pseudotimes attributed to them. I don't know what that means. I also need to better understand the significance of the branch points. Let me know if you have any questions! That will help guide what I learn about!  
