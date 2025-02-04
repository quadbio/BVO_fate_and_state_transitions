library(Seurat)
library(simspec)
library(dplyr)
library(presto)
source("~/Tools/scripts/feature_plots.r")

# endothelial cell heterogeneity
seurat <- readRDS("H9_AV_merged_filtered_css_EC_only.rds")
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
blacklist <- unique(c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1]))
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
seurat <- ScaleData(seurat, vars.to.regress = c("S.Score","G2M.Score")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
FeaturePlot(seurat, features=c("CXCR4","GJA4","IGFBP3","ARL4A","CPE","NR2F2","COL15A1","TM4SF18","PLVAP"), reduction="umap", order=T, cols=beach_colscheme(30)) & NoAxes() & NoLegend()

seurat <- cluster_sim_spectrum(seurat, label_tag = "group", cluster_resolution = 1) %>%
  run_PCA(reduction = "css", npcs = 20, reduction.name="css_pca", reduction.key="CSSPCA_")
cor(Embeddings(seurat, "css_pca"), seurat@meta.data[,c("S.Score","G2M.Score")])
cor(Embeddings(seurat, "css_pca"), seurat@meta.data[,c("S.Score","G2M.Score")], method="spearman")

seurat[['css_pca_nocc']] <- CreateDimReducObject(apply(Embeddings(seurat, "css_pca"), 2, function(x) residuals(lm(x ~ seurat$G2M.Score+seurat$S.Score))), key="CSSPCANOCC_", assay=DefaultAssay(seurat))

seurat <- RunUMAP(seurat, reduction = "css_pca_nocc", dims = c(1:20), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
FeaturePlot(seurat, features=c("CXCR4","GJA4","IGFBP3","UNC5B","ARL4A","CPE","NR2F2","DUSP23","COL15A1","TM4SF18","PLVAP","NKX2-3","CLDN5","DCN","HAND1","MKI67"), reduction="umap_css", order=T, cols=beach_colscheme(30)) & NoAxes() & NoLegend()

seurat <- FindNeighbors(seurat, reduction="css_pca_nocc", dims = 1:20)
seurat[['RNA_css_nn']] <- seurat[['RNA_nn']]
seurat[['RNA_css_snn']] <- seurat[['RNA_snn']]
seurat <- FindClusters(seurat, graph.name = "RNA_css_snn", resolution = 0.1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 0.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 2)
DimPlot(seurat, group.by="RNA_css_snn_res.0.1", reduction="umap_css", label=T) & NoLegend() & NoAxes()
DimPlot(seurat, group.by="RNA_css_snn_res.0.5", reduction="umap_css", label=T) & NoLegend() & NoAxes()

saveRDS(seurat, file="H9_AV_EC_processed.seurat.rds")

DE_cl_res0.5 <- wilcoxauc(seurat, "RNA_css_snn_res.0.5")
DE_cl_res1.5 <- wilcoxauc(seurat, "RNA_css_snn_res.1.5")
markers_cl_res0.5 <- lapply(split(DE_cl_res0.5, DE_cl_res0.5$group), function(x){
  x <- x[which(x$padj < 0.01 & x$logFC > log(1.2) & x$auc > 0.6 & x$pct_out < 30 & x$pct_in - x$pct_out > 20),]
  x <- x[order(x$pct_in - x$pct_out, decreasing=T),]
  return(x)
})
markers_cl_res1.5 <- lapply(split(DE_cl_res1.5, DE_cl_res1.5$group), function(x){
  x <- x[which(x$padj < 0.01 & x$logFC > log(1.2) & x$auc > 0.6 & x$pct_in > 20 & x$pct_out < 20 & x$pct_in - x$pct_out > 20),]
  x <- x[order(x$pct_out),]
  return(x)
})
FeaturePlot(seurat, features=c("SLC1A3","ANKRD1","PDGFRB","COL6A3"), reduction="umap_css", order=T, cols=beach_colscheme(30)) & NoAxes() & NoLegend()
FeaturePlot(seurat, features=c("PDGFRB","COL1A2","ACTA2","CLDN5"), reduction="umap_css", order=T, cols=beach_colscheme(30)) & NoAxes() & NoLegend()





# further subset to clean up the PDGFRB+ population
seurat <- subset(seurat, subset = RNA_css_snn_res.0.1 != "2")
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
seurat <- ScaleData(seurat, vars.to.regress = c("S.Score","G2M.Score")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
seurat <- cluster_sim_spectrum(seurat, label_tag = "group", cluster_resolution = 1) %>%
  run_PCA(reduction = "css", npcs = 20, reduction.name="css_pca", reduction.key="CSSPCA_")
seurat[['css_nocc']] <- CreateDimReducObject(apply(Embeddings(seurat, "css"), 2, function(x) residuals(lm(x ~ seurat$G2M.Score+seurat$S.Score))), key="CSSNOCC_", assay=DefaultAssay(seurat))
seurat <- run_PCA(seurat, reduction = "css_nocc", npcs = 20, reduction.name="css_nocc_pca", reduction.key="CSSNOCCPCA_")
seurat[['css_pca_nocc']] <- CreateDimReducObject(apply(Embeddings(seurat, "css_pca"), 2, function(x) residuals(lm(x ~ seurat$G2M.Score+seurat$S.Score))), key="CSSPCANOCC_", assay=DefaultAssay(seurat))

seurat <- FindNeighbors(seurat, reduction="css_nocc_pca", dims = 1:10)
seurat[['RNA_css_nn']] <- seurat[['RNA_nn']]
seurat[['RNA_css_snn']] <- seurat[['RNA_snn']]
seurat <- FindClusters(seurat, graph.name = "RNA_css_snn", resolution = 0.1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 0.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 2)
DimPlot(seurat, group.by="RNA_css_snn_res.0.1", reduction="umap_css", label=T) & NoLegend() & NoAxes()

saveRDS(seurat, file="H9_AV_EC_PDGFRB_cells_excl.seurat.rds")

freq <- table(seurat$RNA_css_snn_res.0.5, seurat$condition)[hcl_cl$order,c(1,3,2,4,6,5)]
barplot(apply(freq,1,function(x)x/sum(x)), col = prettyrainbow_colscheme(6))

freq <- table(seurat$Phase, seurat$condition)[c(1,3,2),c(1,3,2,4,6,5)]
layout(matrix(1:2,nrow=1))
barplot(freq, col=c(G1="#cdcdcd",S="#909090",G2M="#303030"))
barplot(apply(freq,2,function(x)x/sum(x)), col=c(G1="#cdcdcd",S="#909090",G2M="#303030"))




# signaling target gene expression in conditions
genes <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", # Notch receptors
           "JAG1","JAG2", "DLL1","DLL3","DLL4", # Notch ligands
           "RBPJ", "SNW1", "MAML1", # nuclear effectors
           "HES1", "HES2", "HES4", "HES5", "HEY1", "HEY2", # target genes
           "DLL4", "NRP1", "KDR", "FLT4", # tip cells
           "CXCR4", "GJA4", "NR2F2") 
avg_expr_cond <- sapply(levels(seurat$condition), function(cond) rowMeans(seurat@assays$RNA@data[genes,which(seurat$condition == cond)]))
mat_norm <- apply(avg_expr_cond, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)))
heatmap.2(t(mat_norm), Rowv=NA, Colv=NA, dendrogram="none", trace="none", key=F, keysize=0.5, col=bluered_colscheme(30), scale="row", margins=c(7,9), cexRow=1.3, cexCol=1)

