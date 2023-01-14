## ------------------------------------
##
## Script name: feature plot
##
## Purpose of script: Map continuous and categorical features to the lower-dimensional embedding space of any seurat project
##
## Author: Jiekun (Jackie) Yang
##
## Date Created: 2022-12-14
##
## Email: jkyang@mit.edu
##
##--------------------------------------

### object_path is the path to the seurat object that needed the metadata
### features is the feature name that you want to map
### var_class is the class of the feature, being either categorical or continuous
### reduction is the embedding space, being either tsne or umap
### feature_colors is a list with colors specified for the feature
### target_folder is where you want to save the plots
### file_format is the format of the generated plot being either png or pdf
### label is a logical argument specifying if you want to have labels on the plot with a default of FALSE
### legend is a categorical argument specifying if you want to include legend or not or have the legend generate separately

## library setup
library(Seurat)
options(future.globals.maxSize = 128000 * 1024^2)
library(future)
library(future.apply)
library(qs)
plan("multiprocess", workers = 16)

## function
feature_plot <- function(object_path = NULL, object = NULL, feature = NULL, var_class = c("categorical", "continuous"), reduction = c("tsne", "umap"), feature_colors = NULL, target_folder = NULL, file_format = c("png", "pdf"), label = FALSE, legend = c("none", "include", "separate")) {
  if (!is.null(object_path)) {
    temp.data <- qread(object_path)
    temp.name <- tools::file_path_sans_ext(basename(object_path))
  }
  if (!is.null(object)) {
    temp.data <- get(object)
  }
  if (ncol(temp.data) > 400000) {
    temp.point.size <- 0.001
    temp.label.size <- 18
  } else if (ncol(temp.data) > 100000 & ncol(temp.data) < 400000) {
    temp.point.size <- 0.1
    temp.label.size <- 20
  } else {
    temp.point.size <- 1
    temp.label.size <- 22
  }
  if (file_format == "png") {
    if (legend != "include") {
      temp.width <- 1500
    } else {
      temp.width <- 1600
    }
    png(paste0(target_folder, temp.name, "_", feature, "_", reduction, ".", file_format), width = temp.width, height = 1500)
  } else {
    if (legend != "include") {
      temp.width <- 8
    } else {
      temp.width <- 9
    }
    pdf(paste0(target_folder, temp.name, "_", feature, "_", reduction, ".", file_format), width = temp.width, height = 8)
  }
  if (var_class == "categorical") {
    if (!is.null(feature_colors)) {
      if (legend == "include") {
        print(DimPlot(temp.data, pt.size = temp.point.size, reduction = tolower(reduction), group.by = feature, cols = feature_colors, label = label, label.size = temp.label.size, raster = FALSE) + NoAxes())
      } else {
        print(DimPlot(temp.data, pt.size = temp.point.size, reduction = tolower(reduction), group.by = feature, cols = feature_colors, label = label, label.size = temp.label.size, raster = FALSE) + NoAxes() + NoLegend())
      }
    } else {
      if (legend == "include") {
        print(DimPlot(temp.data, pt.size = temp.point.size, reduction = tolower(reduction), group.by = feature, label = label, label.size = temp.label.size, raster = FALSE) + NoAxes())
      } else {
        print(DimPlot(temp.data, pt.size = temp.point.size, reduction = tolower(reduction), group.by = feature, label = label, label.size = temp.label.size, raster = FALSE) + NoAxes() + NoLegend())
      }
    }
  } else {
    if (!is.null(feature_colors)) {
      if (legend == "include") {
        print(FeaturePlot(temp.data, pt.size = temp.point.size, reduction = tolower(reduction), features = feature, cols = feature_colors, label = label, label.size = temp.label.size, min.cutoff = "q1", max.cutoff = "q9", raster = FALSE) + NoAxes())
      } else {
        print(FeaturePlot(temp.data, pt.size = temp.point.size, reduction = tolower(reduction), features = feature, cols = feature_colors, label = label, label.size = temp.label.size, min.cutoff = "q1", max.cutoff = "q9", raster = FALSE) + NoAxes() + NoLegend())
      }
    } else {
      if (legend == "include") {
        print(FeaturePlot(temp.data, pt.size = temp.point.size, reduction = tolower(reduction), features = feature, label = label, label.size = temp.label.size, min.cutoff = "q1", max.cutoff = "q9", raster = FALSE) + NoAxes())
      } else {
        print(FeaturePlot(temp.data, pt.size = temp.point.size, reduction = tolower(reduction), features = feature, label = label, label.size = temp.label.size, min.cutoff = "q1", max.cutoff = "q9", raster = FALSE) + NoAxes() + NoLegend())
      }
    }
  }
  dev.off()
  if(legend == "separate") {
    temp.legend <- data.frame("groups" = unique(temp.data@meta.data[, feature]))
    pdf(paste0(target_folder, feature, "_legend.pdf"), width = 3, height = 3, useDingbats = FALSE)
    if (!is.null(feature_colors)) {
      print(ggplot(temp.legend, aes("1", groups)) +
              geom_point(aes(color = groups)) +
              scale_color_manual(values = feature_colors) +
              theme_classic() +
              theme(legend.position = "None"))
    } else {
      print(ggplot(temp.legend, aes("1", groups)) +
              geom_point(aes(color = groups)) +
              theme_classic() +
              theme(legend.position = "None"))
    }
    dev.off()
  }
}
