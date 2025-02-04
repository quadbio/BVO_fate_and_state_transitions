library(Seurat)
library(simspec)
library(Matrix)
library(dplyr)
library(presto)
source("~/Tools/scripts/feature_plots.r")

seurat <- readRDS("H9_AV_merged_filtered_css.rds")
seurat_endo <- readRDS("H9_AV_EC_PDGFRB_cells_excl.seurat.rds")

freq <- table(seurat$RNA_snn_res.0.8 %in% c(1,3,13,10), seurat$condition)[,c(1,3,2,4,6,5)]
data.frame(t(sapply(2:ncol(freq), function(i){
  test <- fisher.test(cbind(freq[,1],freq[,i]))
  res <- setNames(c(test$estimate, test$p.value), c("oddsratio","pval"))
  return(res)
})), row.names=colnames(freq)[-1])
layout(matrix(1:2,nrow=1))
barplot(freq, col=c("#f1c40f","#e74c3c"))
barplot(apply(freq,2,function(x)x/sum(x)), col=c("#f1c40f","#e74c3c"))

freq <- table(seurat$Phase, seurat$condition)[c(1,3,2),c(1,3,2,4,6,5)]
layout(matrix(1:2,nrow=1))
barplot(freq, col=c(G1="#cdcdcd",S="#909090",G2M="#303030"))
barplot(apply(freq,2,function(x)x/sum(x)), col=c(G1="#cdcdcd",S="#909090",G2M="#303030"))

freq <- table(seurat_endo$Phase, seurat_endo$condition)[c(1,3,2),c(1,3,2,4,6,5)]
layout(matrix(1:2,nrow=1))
barplot(freq, col=c(G1="#cdcdcd",S="#909090",G2M="#303030"))
barplot(apply(freq,2,function(x)x/sum(x)), col=c(G1="#cdcdcd",S="#909090",G2M="#303030"))

saveRDS(seurat, file="H9_AV_merged_filtered_processed.rds")



# MC subset
seurat_mc <- subset(seurat, cells = colnames(seurat)[which(! seurat$RNA_snn_res.0.8 %in% c(1,3,13,10))]) %>%
  FindVariableFeatures(nfeatures = 3000)
blacklist <- unique(c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat_mc), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1]))
VariableFeatures(seurat_mc) <- setdiff(VariableFeatures(seurat_mc), blacklist)
seurat_mc <- ScaleData(seurat_mc, vars.to.regress = c("S.Score","G2M.Score")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
seurat_mc <- cluster_sim_spectrum(seurat_mc, label_tag = "condition", cluster_resolution = 1) %>%
  run_PCA(reduction = "css", npcs = 20, reduction.name="css_pca", reduction.key="CSSPCA_")
seurat_mc[['css_pca_nocc']] <- CreateDimReducObject(apply(Embeddings(seurat_mc, "css_pca"), 2, function(x) residuals(lm(x ~ seurat_mc$G2M.Score+seurat_mc$S.Score))), key="CSSPCANOCC_", assay=DefaultAssay(seurat_mc))

seurat_mc <- RunUMAP(seurat_mc, reduction = "css_pca", dims = c(1:20), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
seurat_mc <- RunUMAP(seurat_mc, reduction = "css_pca_nocc", dims = c(1:20), reduction.name = "umap_css_nocc", reduction.key = "UMAPCSSNOCC_")
DimPlot(seurat_mc, group.by="condition", reduction="umap_css_nocc") & NoAxes()
DimPlot(seurat_mc, group.by="Phase", reduction="umap_css_nocc") & NoAxes()
DimPlot(seurat_mc, group.by="RNA_snn_res.0.8", reduction="umap_css_nocc", label=T) & NoAxes()

seurat_mc <- FindNeighbors(seurat_mc, reduction="css_pca_nocc", dims = 1:20) %>%
  FindClusters(resolution = 0.2) %>%
  FindClusters(resolution = 0.3)

saveRDS(seurat_mc, file="H9_AV_MC_processed.seurat.rds")

freq <- table(seurat_mc$Phase, seurat_mc$condition)[c(1,3,2),c(1,3,2,4,6,5)]
layout(matrix(1:2,nrow=1))
barplot(freq, col=c(G1="#cdcdcd",S="#909090",G2M="#303030"))
barplot(apply(freq,2,function(x)x/sum(x)), col=c(G1="#cdcdcd",S="#909090",G2M="#303030"))


## cluster composition
freq <- table(seurat_mc$RNA_snn_res.0.3, seurat_mc$condition)[names(col_cl_mc),c(1,3,2,4,6,5)]
layout(matrix(1:2,nrow=1)); par(mar=c(5,5,1,1))
barplot(apply(freq,1,function(x)x/sum(x)), col=col_cond)
barplot(apply(freq,2,function(x)x/sum(x)), col=col_cl_mc)


