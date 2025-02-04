library(Seurat)
library(dplyr)
library(presto)
library(simspec)

seurat <- readRDS("../data/timecourse.seurat.rds")
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 5000) %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt")) %>% RunPCA(npcs = 50, verbose = F)

cc_genes <- read.table("../data/regev_lab_cell_cycle_genes.txt", stringsAsFactors=F)[,1]
s_features <- cc_genes[1:43]
g2m_features <- cc_genes[44:97]
g2m_features_more <- read.table("~/Work/databases/GeneLists/G2M_genes.txt", stringsAsFactors=F)[,1]
seurat@meta.data$group <- paste0(seurat@meta.data$line, "_", seurat@meta.data$sample)
seurat <- CellCycleScoring(seurat, s.features = s_features, g2m.features = g2m_features, set.ident = TRUE)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 3000)
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), unique(c(cc_genes, g2m_features_more)))
seurat <- seurat %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt","S.Score", "G2M.Score")) %>%
  RunPCA(npcs = 50, verbose = F)
seurat <- cluster_sim_spectrum(seurat, dims_use = 1:20, label_tag = "group", cluster_resolution = 0.5, merge_spectrums = F, merge_height_prop = 0.02)
seurat <- RunUMAP(seurat, reduction = "css", dims = 1:ncol(seurat@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
DimPlot(seurat, reduction = "umap_css", group.by="sample")
DimPlot(seurat, reduction = "umap_css", group.by="line")
seurat <- seurat %>% FindNeighbors(reduction = "css", dims = 1:ncol(seurat@reductions$css@cell.embeddings)) %>% FindClusters(resolution = 0.6)
seurat <- FindClusters(resolution = 1)
DimPlot(seurat, reduction = "umap_css", label = T)
markers <- wilcoxauc(seurat)
saveRDS(seurat, file="../data/timecourse.seurat.rds")


# RNA velocity
counts_spliced_unspliced <- readRDS("../data/dropEst/spliced_unspliced_count_matrices.rds")
seurat[["spliced"]] <- CreateAssayObject(counts_spliced_unspliced$spliced)
seurat[["unspliced"]] <- CreateAssayObject(counts_spliced_unspliced$unspliced)

library(velocyto.R)
library(SeuratWrappers)
seurat <- RunVelocity(object = seurat, reduction = "css", graph = "RNA_snn", deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = 50)
saveRDS(seurat, file="../data/timecourse.seurat.rds")

vel_embed <- show.velocity.on.embedding.cor(emb = Embeddings(object = seurat, reduction = "umap_css"),
                               vel = Tool(object = seurat, slot = "RunVelocity"),
                               n = 1000, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = TRUE, cell.border.alpha = 0.1)

# Generate loom file for scVelo
library(loomR)
loom <- loomR::create("../data/timecourse_velocity.loom",
                      data = seurat@assays$spliced@data,
                      layers = list(spliced = seurat@assays$spliced@data, unspliced = seurat@assays$unspliced@data),
                      cell.attrs = c(as.list(seurat@meta.data), list(umap_css_embedding = seurat@reductions$umap_css@cell.embeddings, umap_embedding = seurat@reductions$umap@cell.embeddings, pca_embedding = seurat@reductions$pca@cell.embeddings, css_embedding = seurat@reductions$css@cell.embeddings)))
loom$close_all()


## after scVelo, save the result into the seurat object
library(Matrix)
library(hdf5r)
hdf5_scvelo <- h5file("../data/timecourse_velocity.h5ad", mode = "r")
hdf5_scvelo$close_all()
trans_mat <- readMM("../data/timecourse_scvelo_transition_mat.mtx")
rownames(trans_mat) <- rownames(seurat@meta.data)
colnames(trans_mat) <- rownames(seurat@meta.data)

cl_avg_tm_to <- sapply(levels(seurat@meta.data$RNA_snn_res.1), function(cl) colMeans(trans_mat[which(seurat@meta.data$RNA_snn_res.1 == cl),]))
cl_avg_tm_from <- sapply(levels(seurat@meta.data$RNA_snn_res.1), function(cl) rowMeans(apply(trans_mat[,which(seurat@meta.data$RNA_snn_res.1 == cl)],2,function(x)x/sum(x)), na.rm=T))



# markers and velocity genes
library(presto)

## markers
markers <- wilcoxauc(seurat, group_by = "RNA_snn_res.1")
markers_cl <- setNames(lapply(unique(markers$group), function(g){
  DE <- markers[markers$group == g,]
  idx <- which(DE$padj < 0.01 & DE$pct_in > 20 & DE$logFC > log(1.3))
  m_cl <- DE$feature[idx[order(DE$auc[idx], decreasing=T)]]
}), unique(markers$group))
cl_patterns <- sapply(markers_cl, function(x) colSums(seurat@assays$RNA@data[x,]))

### endothelial cells end-state: comparing C4, C16, C21 and C6
DE_endo_ends <-wilcoxauc(X = seurat@assays$RNA@data[,seurat@meta.data$RNA_snn_res.1 %in% c("4","16","21","6")],
                         y = as.character(seurat@meta.data[seurat@meta.data$RNA_snn_res.1 %in% c("4","16","21","6"), "RNA_snn_res.1"]))
endo_ends_markers <- setNames(lapply(unique(DE_endo_ends$group), function(x){
  idx <- which(DE_endo_ends$group == x & DE_endo_ends$padj < 0.01 & DE_endo_ends$auc > 0.7 & DE_endo_ends$logFC > log(1.5) & DE_endo_ends$pct_in > 20)
  return(DE_endo_ends$feature[idx[order(DE_endo_ends$auc[idx], decreasing=T)]])
}), unique(DE_endo_ends$group))
endo_ends_markers <- lapply(endo_ends_markers, function(x) setdiff(x, unlist(markers_cl[setdiff(names(markers_cl),as.character(c(4,16,21,6,13,10,20)))])))

### pericytes end-state: comparing C0, C3 and C14
DE_peri_ends <-wilcoxauc(X = seurat@assays$RNA@data[,seurat@meta.data$RNA_snn_res.1 %in% c("0","3","14")],
                         y = as.character(seurat@meta.data[seurat@meta.data$RNA_snn_res.1 %in% c("0","3","14"), "RNA_snn_res.1"]))
peri_ends_markers <- setNames(lapply(unique(DE_peri_ends$group), function(x){
  idx <- which(DE_peri_ends$group == x & DE_peri_ends$padj < 0.01 & DE_peri_ends$auc > 0.7 & DE_peri_ends$logFC > log(1.5) & DE_peri_ends$pct_in > 20)
  return(DE_peri_ends$feature[idx[order(DE_peri_ends$auc[idx], decreasing=T)]])
}), unique(DE_peri_ends$group))



## pan-endothelial markers and pan-pericyte markers
library(presto)
idx_endo <- which(seurat@meta.data$RNA_snn_res.1 %in% c("13","16","4","6","21"))
idx_peri <- which(seurat@meta.data$RNA_snn_res.1 %in% c("14","3","5","12","0","17"))
idx_meso <- which(seurat@meta.data$RNA_snn_res.1 %in% c("7","15","8"))
idx_psc <- which(seurat@meta.data$RNA_snn_res.1 %in% c("1","2","9"))
DE_major <- wilcoxauc(X = seurat@assays$RNA@data[,c(idx_psc, idx_meso, idx_peri, idx_endo)],
                      y = rep(c("PSC","Meso","Peri","Endo"), c(length(idx_psc), length(idx_meso), length(idx_peri), length(idx_endo))))
major_markers <- setNames(lapply(unique(DE_major$group), function(g){
  DE <- DE_major[DE_major$group == g,]
  pct_in_other <- tapply(DE_major$pct_in[DE_major$group!=g], DE_major$feature[DE_major$group!=g], max)
  idx <- which(DE$padj < 0.01 & DE$logFC > log(1.3) & DE$pct_in - pct_in_other[DE$feature] > 30)
  DE$feature[idx[order(DE$logFC[idx], decreasing=T)]]
}), unique(DE_major$group))
major_TFs <- lapply(major_markers, function(x) intersect(TFs,x))
major_receptors <- lapply(major_markers, function(x) intersect(receptors,x))

