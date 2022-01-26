# library
library(DESeq2)
library(IHW)
library(rgl)
library(pheatmap)
library(RColorBrewer)
hm_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
annotation_colors <- list(diet = c(HFD = "#E41A1C", Chow = "#377EB8"), training = c(Training = "#4DAF4A", Sedentary = "#984EA3"), tissue = c(SkM = "#1B9E77", scWAT = "#D95F02", vWAT = "#7570B3"))

# data import
bulk_data <- read.table("bulk_sample_level_data.txt", header = T, sep = "\t", stringsAsFactors = F, check.names = F)
# data filter and metadata setup
cts <- bulk_data[bulk_data$`Gene type` == "protein_coding", grepl("intCt|MGI symbol", colnames(bulk_data))]
cts <- cts[!duplicated(cts$`MGI symbol`) & !is.na(cts$`MGI symbol`),]
rownames(cts) <- cts$`MGI symbol`
cts <- cts[, -ncol(cts)]
coldata <- as.data.frame(colnames(cts))
colnames(coldata) <- "sample"
coldata$tissue <- sapply(strsplit(coldata$sample, "\\."), `[`, 2)
coldata$pheno_class <- substr(sapply(strsplit(coldata$sample, "\\."), `[`, 1), 1, 2)
coldata$rep <- substr(sapply(strsplit(coldata$sample, "\\."), `[`, 1), 3, length(sapply(strsplit(coldata$sample, "\\."), `[`, 1)))
coldata <- as.data.frame(coldata %>% group_by(tissue, pheno_class) %>% mutate(rep_recode = row_number()))
coldata$diet <- ifelse(coldata$pheno_class %in% c("SC", "TC"), "Chow", "HFD")
coldata$training <- ifelse(coldata$pheno_class %in% c("SC", "SH"), "Sedentary", "Training")
coldata$pheno_class <- factor(coldata$pheno_class, levels = c("SC", "TC", "SH", "TH"))
rownames(coldata) <- coldata$sample
# DESeq2
for (temp.tissue in c("scWAT", "vWAT", "SkM")) {
    dds <- DESeqDataSetFromMatrix(countData = cts[, grepl(temp.tissue, colnames(cts))], colData = coldata[coldata$tissue == temp.tissue,], design = ~ pheno_class)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    ## Data quality assessment
    vsd <- vst(dds, blind = FALSE)
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- vsd$sample
    colnames(sampleDistMatrix) <- vsd$sample
    pdf(paste0(temp.tissue, "_dist_cluster_protein_coding.pdf"), width = 7, height = 5)
    print(pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, annotation_col = coldata[coldata$tissue == temp.tissue, rev(c("tissue", "training", "diet"))], show_rownames = FALSE, show_colnames = FALSE, clustering_method = "ward.D2", border_color = NA, color = hm_colors, annotation_colors = annotation_colors))
    ### diet effect dominates training effect based on unsupervised clustering
    dev.off()
    rv <- rowVars(assay(vsd))
    select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
    pca <- prcomp(t(assay(vsd)[select,]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    pca_anno <- merge(pca$x, coldata[coldata$tissue == temp.tissue,], by.x = "row.names", by.y = "row.names")
    pca_anno$diet <- factor(pca_anno$diet, levels = names(annotation_colors$diet))
    pca_anno$diet_color <- annotation_colors$diet[as.numeric(pca_anno$diet)]
    pca_anno$training <- factor(pca_anno$training, levels = names(annotation_colors$training))
    pca_anno$training_color <- annotation_colors$training[as.numeric(pca_anno$training)]
    print(ggplot(pca_anno, aes(PC1, PC2)) +
        geom_point(aes(color = pheno_class), size = 2) +
        xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
        ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
        theme_bw() +
        theme(panel.grid = element_blank()))
    print(plot3d(x = pca_anno$PC1, y = pca_anno$PC2, z = pca_anno$PC3, col = pca_anno$training_color, type = "s", radius = .5, xlab = paste0("PC1: ", round(percentVar[1] * 100), "% variance"), ylab = paste0("PC2: ", round(percentVar[2] * 100), "% variance"), zlab = paste0("PC3: ", round(percentVar[3] * 100), "% variance")))
    rgl.postscript(paste0(temp.tissue, "_vsd_pca_pc1-3_protein_coding_top500_training.pdf"),"pdf")
    ## DE analysis
    dds <- DESeq(dds)
    for (temp.comp in resultsNames(dds)[2:3]) {
        res <- results(dds, name = temp.comp, filterFun = ihw)
        resLFC <- lfcShrink(dds, coef = temp.comp, type = "ashr")
        summary(resLFC, alpha = 0.05) # adjusted p < 0.05
        print(plotMA(resLFC, ylim=c(-2,2), alpha = 0.05))
        sig <- as.data.frame(resLFC[resLFC$padj < 0.05 & !is.na(resLFC$padj),])
        sig <- merge(sig, bulk_data[, c("Row.names", "MGI symbol")], by.x = "row.names", by.y = "Row.names")
        write.table(sig, paste0(temp.tissue, "_sig_", temp.comp, "_0.05.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    dds$pheno_class <- relevel(dds$pheno_class, ref = "SH")
    dds <- nbinomWaldTest(dds)
    resultsNames(dds)
    res_TH_SH <- results(dds, name = "pheno_class_TH_vs_SH", filterFun = ihw)
    resLFC_TH_SH <- lfcShrink(dds, coef = "pheno_class_TH_vs_SH", type = "ashr")
    summary(resLFC_TH_SH, alpha = 0.05)
    print(plotMA(resLFC_TH_SH, ylim=c(-2,2), alpha = 0.05))
    sig_TH_SH <- as.data.frame(resLFC_TH_SH[resLFC_TH_SH$padj < 0.05 & !is.na(resLFC_TH_SH$padj),])
    sig_TH_SH <- merge(sig_TH_SH, bulk_data[, c("Row.names", "MGI symbol")], by.x = "row.names", by.y = "Row.names")
    write.table(sig_TH_SH, paste0(temp.tissue, "_sig_pheno_class_TH_vs_SH_0.05.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    saveRDS(dds, paste0("dds_", temp.tissue, ".rds"))
}
