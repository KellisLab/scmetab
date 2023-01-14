## ------------------------------------
##
## Script name: fancy_tsne
##
## Purpose of script: Implementation of suggestions offered in Nature Communication "The art of using t-SNE for single-cell transcriptomics". The script covers data normalization, selected highly variable gene removal, PCA, tSNE, UMAP, density- and hierachical-based clustering
##
## Author: Jackie Yang
##
## Date Created: 2021-01-13
##
## Email: jkyang@mit.edu
##
##--------------------------------------

### data is the data object
### rm_mt_hsp is a logical argument indicating if reads mapped to mitochondrial genome or genes encoding heatshock proteins are removed with a default of false
### species is among "mouse" and "human"
### norm_mode is among "regular" (median normalization and log transformation) or "sctransform"
### cluster is a logical argument indicating if neighborhood identification and hierachical- and density-based clustering are performed with a default of false
### nthreads is the number of threads used for this function
### dims is the number of PCs used for running the UMAP, neighbor and clustering functions
### res is the resolution used for the clustering function

## library setup
library(Seurat)
library(dbscan)
fftRtsne_dir <- "~/data-lab5/programs/FIt-SNE/"
source(paste0(fftRtsne_dir, "fast_tsne.R"))

fancy_tsne <- function(data, very_large, rm_mt_hsp, species, norm_mode, cluster, nthreads, dims, res) {
  if (DefaultAssay(data) != "integrated") {
    if (norm_mode == "regular") {
        data <- NormalizeData(data, scale.factor = median(data$nCount_RNA))
        data <- ScaleData(data, verbose = FALSE)
    } else if (norm_mode == "sctransform") {
        data <- SCTransform(data, verbose = FALSE)
    }
    data <- FindVariableFeatures(data)
  } else {
    data <- ScaleData(data, verbose = FALSE)
  }
  if (rm_mt_hsp) {
      # exclude MT/heat_shock/hypoxia genes from the highly variable genes used for calculating PCs and downstream analyses
      mtGenes <- grep("^mt-", rownames(data), value = T, ignore.case = T)
      hspGenes <- as.character(read.csv(paste0("~/data-lab5/resource/heat_shock_protein_gene_list_", species, ".csv"), header = F)$V1)
      if (norm_mode == "regular") {
          data@assays$RNA@var.features <- setdiff(VariableFeatures(data), c(mtGenes, hspGenes))
      } else if (norm_mode == "sctransform") {
          data@assays$SCT@var.features <- setdiff(VariableFeatures(data), c(mtGenes, hspGenes))
      }
  }
  data <- RunPCA(data, features = VariableFeatures(data)) # default 50 pcs
  PCA.init <- data@reductions$pca@cell.embeddings[, 1:2] / sd(data@reductions$pca@cell.embeddings[, 1]) * 0.0001
  if (very_large) {
      temp <- fftRtsne(data@reductions$pca@cell.embeddings, perplexity = 30, initialization = PCA.init, learning_rate = ncol(data)/12, max_iter = 1000, nthreads = nthreads, exaggeration_factor = 4)
  } else {
      temp <- fftRtsne(data@reductions$pca@cell.embeddings, perplexity_list = c(30, floor(ncol(data)/100)), perplexity = 0, initialization = PCA.init, learning_rate = ncol(data)/12, max_iter = 1000, nthreads = nthreads)
  }
  colnames(temp) <- paste0("tSNE_", 1:2)
  rownames(temp) <- colnames(data)
  data[["tsne"]] <- CreateDimReducObject(embeddings = temp, key = "tSNE_", assay = DefaultAssay(data))
  data <- RunUMAP(data, reduction = "pca", dims = 1:dims)
  if (cluster) {
      data <- FindNeighbors(data, reduction = "pca", dims = 1:dims)
      data <- FindClusters(data, res = res)
      data$umap_dbscan <- dbscan(Embeddings(object = data[["umap"]]), eps = res)$cluster
  }
  return(data)
}