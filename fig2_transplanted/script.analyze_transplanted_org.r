library(Seurat)
library(dplyr)

# read data, create Seurat object and QC
dat <- readRDS("../../data/data.counts_meta_transplanted_human.rds")
seurat <- CreateSeuratObject(dat$count, project="transplanted", meta.data=dat$meta)
seurat@meta.data$orig.ident <- "transplanted"
seurat@active.ident <- as.factor(setNames(seurat$orig.ident, colnames(seurat)))
seurat[['percent.mt']] <- PercentageFeatureSet(seurat, pattern = "^MT\\.")
seurat <- subset(seurat, subset = nFeature_RNA > 1500 & nFeature_RNA < 6000 & percent.mt < 15)

# routine preprocessing
seurat <- NormalizeData(seurat) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress = c("nCount_RNA","percent.mt")) %>% RunPCA(npcs = 50) %>% RunUMAP(dims = 1:20)
seurat$group <- paste0(seurat$sample, "_", seurat$line)
saveRDS(seurat, file="../../data/data.seurat_transplanted_human.rds")

## feature plots
source("~/Tools/scripts/feature_plots.r")
layout(matrix(1:4, nrow=2, byrow=T))
par(mar=c(1,1,3,1))
plotFeature(Embeddings(seurat, "umap"), seurat$sample, cex=0.5, axis=F, main="Sample", cex.main=1.5, colorPal=div_colscheme)
legend("topright", bty="n", border=NA, fill=div_colscheme(length(unique(seurat$sample))), legend=levels(as.factor(seurat$sample)))
plotFeature(Embeddings(seurat, "umap"), seurat$line, cex=0.5, axis=F, main="Line", cex.main=1.5, colorPal=div_colscheme)
legend("topright", bty="n", border=NA, fill=div_colscheme(length(unique(seurat$line))), legend=levels(as.factor(seurat$line)))
for(gene in c("CLDN5","COL1A2")) plotFeature(Embeddings(seurat,"umap"), as.numeric(seurat@assays$RNA@data[gene,]), cex=0.5, axis=F, main = gene, cex.main=1.5, colorPal=beach_colscheme)

# CSS integration of reference samples, defined as sample-line with more than 500 cells
source("~/Tools/scripts/css_integration.r")
source("~/Tools/scripts/dimension_reduction.r")
ref_groups <- names(table(seurat$group))[table(seurat$group)>500]
seurat_ref <- subset(seurat, subset = group %in% ref_groups)
ref_css_model <- cluster_sim_spectrum(seurat_ref, label_tag="group", cluster_resolution=0.6, return_seuratObj=F)
seurat_ref[['CSS']] <- CreateDimReducObject(ref_css_model$sim2profiles, key="CSS_")
ref_umapcss_model <- get_umap_model(seurat_ref, reduction="CSS", dims = 1:ncol(Embeddings(seurat_ref, "CSS")))
seurat_ref[["umap_CSS"]] <- CreateDimReducObject(ref_umapcss_model$embedding, key="UMAPCSS_")
rownames(seurat_ref@reductions$umap_CSS@cell.embeddings) <- colnames(seurat_ref)
saveRDS(seurat_ref, file="../../data/data.seurat_transplanted_human_ref.rds")

## save CSS integration models, which will be used for non-reference cell projection for a unified embedding and cell representation
save_uwot(ref_umapcss_model, file="~/Work/vascular_organoids/analysis/transplanted/res.umap_model_ref-CSS.uwot")
saveRDS(ref_css_model, file="~/Work/vascular_organoids/analysis/transplanted/res.css_model_ref.rds")

## feature plots
layout(matrix(1:4, nrow=2, byrow=T))
plotFeature(Embeddings(seurat_ref, "umap_CSS"), seurat_ref$sample, cex=0.5, axis=F, main="Sample", cex.main=1.5, colorPal=div_colscheme)
plotFeature(Embeddings(seurat_ref, "umap_CSS"), seurat_ref$line, cex=0.5, axis=F, main="Line", cex.main=1.5, colorPal=div_colscheme)
for(gene in c("CLDN5","COL1A2")) plotFeature(Embeddings(seurat_ref,"umap_CSS"), as.numeric(seurat_ref@assays$RNA@data[gene,]), cex=0.5, axis=F, main = gene, cex.main=1.5, colorPal=beach_colscheme)

# CSS representation of cells in non-reference samples and the projected UMAP
query_groups <- setdiff(unique(seurat$group), ref_groups)
seurat_query <- subset(seurat, subset = group %in% query_groups)
rss2ref <- do.call(cbind, lapply(ref_css_model$model$profiles, function(ref) RefSimSpec::calculateRSS(as.matrix(seurat_query@assays$RNA@data), ref)))
umapcss_query <- umap_transform(rss2ref, ref_umapcss_model)
rownames(umapcss_query) <- rownames(rss2ref)

seurat[['CSS']] <- CreateDimReducObject(rbind(Embeddings(seurat_ref, "CSS"), rss2ref)[colnames(seurat),], key="CSS")
seurat[['umap_CSSproj']] <- CreateDimReducObject(rbind(Embeddings(seurat_ref, "umap_CSS"), umapcss_query)[colnames(seurat),], key="UMAPCSSPROJ_")
saveRDS(seurat, file="../../data/data.seurat_transplanted_human.rds")

## feature plots
layout(matrix(1:4, nrow=2, byrow=T))
par(mar=c(1,1,3,1))
plotFeature(Embeddings(seurat, "umap_CSSproj"), seurat$sample, cex=0.5, axis=F, main="Sample", cex.main=1.5, colorPal=div_colscheme)
legend("topleft", bty="n", border=NA, fill=div_colscheme(length(unique(seurat$sample))), legend=levels(as.factor(seurat$sample)))
plotFeature(Embeddings(seurat, "umap_CSSproj"), seurat$line, cex=0.5, axis=F, main="Line", cex.main=1.5, colorPal=div_colscheme)
legend("topleft", bty="n", border=NA, fill=div_colscheme(length(unique(seurat$line))), legend=levels(as.factor(seurat$line)))
for(gene in c("CLDN5","COL1A2")) plotFeature(Embeddings(seurat,"umap_CSSproj"), as.numeric(seurat@assays$RNA@data[gene,]), cex=0.5, axis=F, main = gene, cex.main=1.5, colorPal=beach_colscheme)

