## ------------------------------------
##
## Script name: add metadata
##
## Purpose of script: Integrate samples by merging/concatenating or rPCA or CCA
##
## Author: Jiekun (Jackie) Yang
##
## Date Created: 2022-12-14
##
## Email: jkyang@mit.edu
##
##--------------------------------------

### object_path is the path to the seurat object that needed the metadata
### object is the seurat object name if already loaded in the environment
### metadata is a data frame with the metadata to be added

## library setup
library(Seurat)
options(future.globals.maxSize = 128000 * 1024^2)
library(future)
library(future.apply)
library(qs)
plan("multiprocess", workers = 16)

## function
add_metadata <- function(object_path = NULL, object = NULL, metadata = NULL) {
  if (!is.null(object_path)) {
    temp.data <- qread(object_path)
  }
  if (!is.null(object)) {
    temp.data <- get(object)
  }
  temp.data@meta.data <- cbind(temp.data@meta.data, pheno_df[match(temp.data$orig.ident, rownames(pheno_df)),])
  if (!is.null(object_path)) {
    qsave(temp.data, object_path)
  }
  return(temp.data)
}
