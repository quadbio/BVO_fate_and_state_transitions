library(Matrix)
library(Seurat)
library(simspec)
library(dplyr)
library(presto)

# data integration and clustering
seurat <- readRDS('../data/transplanted.seurat.rds')
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 3000)
blacklist <- c(read.table("../data/RPgenes.txt")[,1], grep("^MT-", rownames(seurat),value=T), unlist(cc.genes.updated.2019), read.table("../data/G2M_genes.txt")[,1])
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
seurat <- CellCycleScoring(seurat, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes)
seurat <- ScaleData(seurat, vars.to.regress = c("G2M.Score", "S.Score")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
seurat <- cluster_sim_spectrum(seurat, label_tag="group", cluster_resolution = 1) %>%
  run_PCA(reduction = "css", npcs = 20, reduction.name="css_pca", reduction.key="CSSPCA_") %>%
  RunUMAP(reduction = "css_pca", dims = 1:20, reduction.name="umap_css", reduction.key="UMAPCSS_")
seurat <- FindNeighbors(seurat, reduction="css_pca", dims = 1:20)
seurat[['RNA_css_nn']] <- seurat[['RNA_nn']]
seurat[['RNA_css_snn']] <- seurat[['RNA_snn']]
seurat <- FindClusters(seurat, graph.name="RNA_css_snn", resolution = 0.1) %>%
  FindClusters(graph.name="RNA_css_snn", resolution = 0.5) %>%
  FindClusters(graph.name="RNA_css_snn", resolution = 0.8) %>%
  FindClusters(graph.name="RNA_css_snn", resolution = 1)

FeaturePlot(seurat, c("HAND1","DCN","TCF21","PDGFRB","ACTA2","COL1A2","VEGFD","SIX1","CLDN5","CXCR4","APLNR","PROX1","EPCAM","SOX10","MKI67"), reduction="umap_css", cols=beach_colscheme(30), pt.size=0.1, order=T) & NoLegend() & NoAxes()

col_cl_2 <- setNames(colorRampPalette(c("#A4D371", "#A64A97", "#CA6778", "#DCCB7C", "#5E4FA2", "#3C7AB6", "#138D75", "#80B1D3", "#4DAA99", "#52BE80","#756BB1", "#322A84", "#F8C471", "#7BCAA4", "#9E0142"))(length(levels(seurat$RNA_css_snn_res.0.8))),
                     levels(seurat$RNA_css_snn_res.0.8)[hclust(as.dist(1-cor(avg_expr_cl)), method="ward.D2")$order])

saveRDS(seurat, file="../data/transplanted.seurat.rds")


## compare to fetal atlas (Yu et al.)
avg_expr_cl_ref <- readRDS("~/Work/intestine_organoids/dat.fetal_atlas_avg_expr_cl_highvar.rds")
cl_ref_cols <- readRDS("~/Work/intestine_organoids/dat.fetal_atlas_cl_cols.rds")
avg_expr_cl <- sapply(levels(seurat$RNA_css_snn_res.0.8), function(cl) rowMeans(seurat@assays$RNA@data[,which(seurat$RNA_css_snn_res.0.8 == cl)]))
sim2ref <- ref_sim_spectrum(avg_expr_cl[VariableFeatures(seurat),], avg_expr_cl_ref, scale=F)
rownames(sim2ref) <- colnames(avg_expr_cl)
colnames(sim2ref) <- colnames(avg_expr_cl_ref)
idx_ref_cl <- c(grep("Mesenchymal",colnames(sim2ref))[c(2,9,6,7,4,5,1,3,8)], NA,
                grep("Endothelial",colnames(sim2ref)), NA,
                grep("Epithelial",colnames(sim2ref))[c(7,3,2,5,4,1,6,8)], NA,
                grep("Neuronal",colnames(sim2ref))[c(2,1)], NA,
                grep("Immune",colnames(sim2ref))[c(1,5,3,4,2)])
idx_cl <- c(c(6,18,9),c(3,7,12,4,10,1,13,15,14,0,17,2,8,11,5,16))
idx_cl <- names(col_cl_2)
pdf("plot.transplanted_all_sim2ref_guttube.pdf")
gplots::heatmap.2(sim2ref[as.character(idx_cl),idx_ref_cl],
                  Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", key=F, keysize=0.2, col=bluered_colscheme(30),
                  #RowSideColors=col_cl[as.character(idx_cl)],
                  #RowSideColors=setNames(colorRampPalette(RColorBrewer::brewer.pal(8,"Set1"))(length(levels(seurat$RNA_css_snn_res.0.8))), levels(seurat$RNA_css_snn_res.0.8))[as.character(idx_cl)],
                  RowSideColors=col_cl_2[as.character(idx_cl)],
                  ColSideColors = cl_ref_cols[colnames(sim2ref)[idx_ref_cl]])
dev.off()





# focus on endothelial cells
seurat_endo <- subset(seurat, subset = RNA_css_snn_res.0.8 %in% c(6,9)) # cluster 18 is likely doublet
seurat_endo <- FindVariableFeatures(seurat_endo, nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(dims = 1:10)
seurat_endo <- cluster_sim_spectrum(seurat_endo, label_tag="group", cluster_resolution = 1, use_fast_rank = T, dims_use = 1:20)
seurat_endo <- RunUMAP(seurat_endo, reduction="css", dims = 1:ncol(Embeddings(seurat_endo, "css")), reduction.name="umap_css", reduction.key="UMAPCSS_")
seurat_endo_mnn <- RunFastMNN(SplitObject(seurat_endo, "group"))
seurat_endo[['mnn']] <- CreateDimReducObject(Embeddings(seurat_endo_mnn,"mnn"), key="MNN_")
seurat_endo <- RunUMAP(seurat_endo, reduction="mnn", dims = 1:10, reduction.name="umap_mnn", reduction.key="UMAPMNN_")

seurat_endo <- FindNeighbors(seurat_endo, reduction="css", dims = 1:ncol(Embeddings(seurat_endo, "css"))) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.2) %>%
  FindClusters(resolution = 0.3) %>% 
  FindClusters(resolution = 0.5)

saveRDS(seurat_endo, file="../data/data.seuratlanted_endothelial.rds")

FeaturePlot(seurat_endo, c("CXCR4","EFNB2","DLL4","GJA5","APLNR","NR2F2","EPHB4","CPE","MKI67"), cols=beach_colscheme(30), reduction="umap_css", order=T) & NoLegend() & NoAxes()
png("plot.transplanted_endo_umap_featureplots.png", height=26, width=24, unit="cm", res=500)
plotMultiFeatures(cbind(-Embeddings(seurat_endo,"umap_css")[,1],Embeddings(seurat_endo,"umap_css")[,2]),
                  seurat_endo@assays$RNA@data[c("CXCR4","EFNB2","DLL4","GJA5","APLNR","NR2F2","CLDN5","CPE","MKI67"),],
                  colorPal = blue_colscheme, cex=1.2, random_order=F, sort_by_value=T, cex.main = 2, mask_col = "#cdcdcd")
dev.off()

par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat,"umap_css"), seurat$RNA_css_snn_res.0.8, colorPal=col_cl, cex=0.6, pt_border=T, lwd_border=0.15, emphasize=which(seurat$RNA_css_snn_res.0.8 %in% c(6,9)), mask_col="#dedede")
plotFeature(Embeddings(seurat_endo,"umap_css"), seurat_endo$RNA_css_snn_res.0.8, colorPal=col_cl, cex=1, pt_border=T, lwd_border=0.3)
col_endo_cl <- setNames(c("#2980B9","#5DADE2","#EC7063","#C0392B","#1ABC9C"),c(3,2,1,0,4))
coord_endo_cl <- t(sapply(levels(seurat_endo$RNA_snn_res.0.2), function(x) colMeans(Embeddings(seurat_endo, "umap_css")[which(seurat_endo$RNA_snn_res.0.2 == x),])))
png("plot.transplanted_endo_umap_clusters.png", height=8, width=8, unit="cm", res=500); par(mar=c(1,1,1,1), cex=0.5)
plotFeature(Embeddings(seurat_endo,"umap_css"), seurat_endo$RNA_snn_res.0.2, colorPal=col_endo_cl, cex=1.5, pt_border=T, lwd_border=0.2, do_label = F, label_round = T, cex.label = 2.5)
dev.off()


# cluster markers
DE_cl <- wilcoxauc(seurat_endo, group_by = "RNA_snn_res.0.2")
DEG_cl <- lapply(split(DE_cl, DE_cl$group), function(x){
  x <- x[which(x$padj < 0.01 & x$auc > 0.65 & x$logFC > log(1.2) & x$pct_in - x$pct_out > 20 & x$pct_out < 20),]
  return(x[order(x$pct_in - x$pct_out, decreasing=T),])
})
features2plot <- unique(unlist(lapply(DEG_cl, function(x) x$feature[1:min(c(nrow(x),5))])))
mat <- sapply(levels(seurat_endo$RNA_snn_res.0.2), function(cl) rowMeans(seurat_endo@assays$RNA@data[features2plot,which(seurat_endo$RNA_snn_res.0.2 == cl)]))
pdf("plot.transplanted_endo_heatmap_cl_markers.pdf")
heatmap.2(apply(mat,1,function(x) (x-min(x))/(max(x)-min(x))), Rowv=NA, Colv=NA, dendrogram="none", trace="none", key=F, keysize=0.5, col=greyscale_colscheme(30), scale="none", margins=c(7,9), cexRow=1.3, cexCol=1)
dev.off()




# compare to fetal endothelial cells (Yu et al.)
avg_expr_endo <- sapply(as.character(c(4,3,2,1,0)), function(cl) rowMeans(seurat_endo@assays$RNA@data[,which(seurat_endo$RNA_snn_res.0.2 == cl)]))
seurat_invitro_endo <- readRDS("../data/data.seurat_untransplanted_noWTC_endoBranch_cleanedup.rds")
avg_expr_invitro_endo <- sapply(as.character(c(18,11,5,6,20)), function(cl) rowMeans(seurat_invitro_endo@assays$RNA@data[,which(seurat_invitro_endo$RNA_css_snn_res.0.8 == cl)]))

seurat_fetal_endo <- readRDS("~/Work/intestine_organoids/fetal_atlas_endothelial.seurat.rds")
seurat_fetal_endo <- subset(seurat_fetal_endo, cells = colnames(seurat_fetal_endo)[which(!seurat_fetal_endo$celltype %in% c("COL1A2+","Lymphatic EC"))]) %>% FindVariableFeatures(nfeatures = 3000)
avg_expr_fetal_endo <- sapply(levels(droplevels(seurat_fetal_endo$celltype)), function(ct) rowMeans(seurat_fetal_endo@assays$RNA@data[,which(seurat_fetal_endo$celltype == ct)]))
ref_fetal_endo_ct <- avg_expr_fetal_endo[setdiff(VariableFeatures(seurat_fetal_endo),blacklist),]

sim2ref_endo <- ref_sim_spectrum(seurat_endo@assays$RNA@data, ref = ref_fetal_endo_ct, scale=F)
dimnames(sim2ref_endo) <- list(colnames(seurat_endo), colnames(ref_fetal_endo_ct))
avg_corr_ref_endo <- t(sapply(as.character(c(4,3,1,2,0)), function(cl) colMeans(sim2ref_endo[which(seurat_endo$RNA_snn_res.0.2 == cl),])))
corr_ref_endo_cl <- ref_sim_spectrum(avg_expr_endo, ref = ref_fetal_endo_ct, scale=F)
dimnames(corr_ref_endo_cl) <- list(colnames(avg_expr_endo), colnames(ref_fetal_endo_ct))

sim2ref_endo_invitro <- ref_sim_spectrum(seurat_invitro_endo@assays$RNA@data, ref = ref_fetal_endo_ct, scale=F)
dimnames(sim2ref_endo_invitro) <- list(colnames(seurat_invitro_endo), colnames(ref_fetal_endo_ct))
corr_ref_invitro_endo_cl <- ref_sim_spectrum(avg_expr_invitro_endo, ref = ref_fetal_endo_ct, scale=F)
dimnames(corr_ref_invitro_endo_cl) <- list(colnames(avg_expr_invitro_endo), colnames(ref_fetal_endo_ct))

corr <- rbind(corr_ref_invitro_endo_cl,NA,corr_ref_endo_cl)
corr_rowNorm <- t(apply(corr, 1, function(x) (x-min(x))/(max(x)-min(x))))
heatmap.2(corr[,c(1,2,4,3,5:8)], Rowv=NA, Colv=NA, dendrogram="none", trace="none", key=F, keysize=0.5, col=bluered_colscheme(30), scale="none", margins=c(10,5), RowSideColors = rep(c("#bdbdbd",NA,"#696969"),c(5,1,5)), cexRow=1.3, cexCol=1.3)
heatmap.2(corr_rowNorm[,c(1,2,4,3,5:8)], Rowv=NA, Colv=NA, dendrogram="none", trace="none", key=F, keysize=0.5, col=bluered_colscheme(30), scale="none", margins=c(10,5))

seurat_fetal_endo_sub <- subset(seurat_fetal_endo, subset = celltype %in% c("Arterial EC","Early aEC","Venous EC","Early vEC"))
#DE_fetal_endo <- wilcoxauc(seurat_fetal_endo, "celltype")
DE_fetal_endo <- wilcoxauc(seurat_fetal_endo_sub, "celltype")
markers_fetal_endo <- lapply(split(DE_fetal_endo,DE_fetal_endo$group), function(x){
  res <- x[which(x$logFC > log(1.2) & x$auc > 0.7 & x$padj < 0.01 & x$pct_in - x$pct_out > 20),]
  return(setdiff(res$feature[order(res$logFC,decreasing=T)], blacklist))
})
genes <- intersect(sort(unique(unlist(lapply(markers_fetal_endo, function(x) x[1:min(c(length(x),10))])))), rownames(seurat_endo))

#mat <- cbind(avg_expr_invitro_endo[genes,], NA, avg_expr_endo[genes,], NA, ref_fetal_endo_ct[genes,c(1,2,4,3,5:8)])
#mat <- mat[order(apply(ref_fetal_endo_ct[genes,c(1,2,4,3,5:8)],1,which.max)),]
mat <- cbind(avg_expr_invitro_endo[genes,], NA, avg_expr_endo[genes,], NA, avg_expr_fetal_endo[genes,c(1,2,4,3)])
mat <- mat[order(apply(avg_expr_fetal_endo[genes,c(1,2,4,3)],1,which.max)),]
mat_norm <- t(apply(mat, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
#heatmap.2(t(mat_norm), Rowv=NA, Colv=NA, dendrogram="none", trace="none", key=F, keysize=0.5, col=bluered_colscheme(30), scale="none", margins=c(7,9), RowSideColors = rep(c("#bdbdbd",NA,"#696969",NA,"#303030"),c(5,1,5,1,8)), cexRow=1.3, cexCol=1)
heatmap.2(t(mat_norm), Rowv=NA, Colv=NA, dendrogram="none", trace="none", key=F, keysize=0.5, col=bluered_colscheme(30), scale="none", margins=c(7,9), RowSideColors = rep(c("#bdbdbd",NA,"#696969",NA,"#303030"),c(5,1,5,1,4)), cexRow=1.3, cexCol=1)


# compare the most venous cluster in transplanted and in vitro organoids
idx_vEC_transp <- which(seurat_endo$RNA_snn_res.0.2 == 3)
idx_vEC_invitro <- which(seurat_invitro_endo$RNA_css_snn_res.0.8 == 18)
mat <- cbind(seurat_endo@assays$RNA@data[,idx_vEC_transp], seurat_invitro_endo@assays$RNA@data[,idx_vEC_invitro])
lab <- rep(c("transplanted","in_vitro"), c(length(idx_vEC_transp), length(idx_vEC_invitro)))
DE_vEC <- wilcoxauc(mat, lab)

DEG_vEC_invitro_transp <- DE_vEC[which(DE_vEC$padj<0.01 & DE_vEC$logFC > log(1.2) & DE_vEC$auc > 0.7 & DE_vEC$pct_in > 30 & DE_vEC$pct_in - DE_vEC$pct_out > 20),]
DEG_vEC_invitro_transp <- DEG_vEC_invitro_transp[which(! DEG_vEC_invitro_transp$feature %in% blacklist),]
DEG_vEC_invitro_transp <- data.frame(do.call(rbind, lapply(split(DEG_vEC_invitro_transp, DEG_vEC_invitro_transp$group), function(x) x[order(x$logFC, decreasing=T),])), row.names=NULL)
write.table(DEG_vEC_invitro_transp, file="tab.DE_vEC_invitro_transp.tsv", quote=F, sep="\t")





# artery-vein trajectory of endothelial cells
seurat_endo <- readRDS("../data/data.seuratlanted_endothelial.rds")
seurat_fetal_endo <- readRDS("~/Work/intestine_organoids/fetal_atlas_endothelial.seurat.rds")

## AV score
library(presto)
DE_av <- wilcoxauc(seurat_fetal_endo, "celltype", groups_use = c("Arterial EC","Venous EC"))
DEG_av <- lapply(split(DE_av, DE_av$group), function(x){
  x <- x[which(x$logFC > log(1.2) & x$auc > 0.7 & x$padj < 0.01 & x$pct_in - x$pct_out > 20 & x$pct_out < 20),]
  return(x[order(x$pct_in - x$pct_out,decreasing=T),])
})
seurat_endo <- AddModuleScore(seurat_endo, features = lapply(DEG_av, function(x) x$feature))
seurat_endo$AV <- seurat_endo$Cluster1 - seurat_endo$Cluster2
saveRDS(seurat_endo, file="../data/data.seuratlanted_endothelial.rds")

FeaturePlot(seurat_endo, "AV", reduction="umap_css", cols=bluewhitered_colscheme(30)) & NoAxes() & NoLegend()
png("plot.transplanted_endo_umap_AV.png", height=8, width=8, unit="cm", res=500); par(mar=c(1,1,1,1), cex=0.6)
plotFeature(Embeddings(seurat_endo,"umap_css"), seurat_endo$AV, colorPal=bluered_colscheme, pt_border = T, lwd_border = 0.1, cex=1.2)
dev.off()

avg_expr_endo_av_bins <- sapply(10:1, function(i) rowMeans(seurat_endo@assays$RNA@data[,ceiling(rank(seurat_endo$AV)/ncol(seurat_endo)*10)==i]))
genes <- intersect(unique(unlist(lapply(DEG_av, function(x) x$feature[1:10]))), rownames(seurat_endo))
mat <- cbind(avg_expr_invitro_endo[genes,], NA, avg_expr_endo_av_bins[genes,], NA, avg_expr_fetal_endo[genes,c(1,3)])
mat <- mat[order(apply(avg_expr_fetal_endo[genes,c(1,2,4,3)],1,which.max)),]
mat_norm <- t(apply(mat, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
pdf("plot.heatmap_av_markers_invitro_transplanted_primary.pdf", height=5)
heatmap.2(t(mat_norm), Rowv=NA, Colv=NA, dendrogram="none", trace="none", key=F, keysize=0.2, col=bluered_colscheme(30), scale="none", margins=c(7,9), RowSideColors = rep(c("#bdbdbd",NA,"#696969",NA,"#303030"),c(5,1,10,1,2)), cexRow=1.3, cexCol=1)
dev.off()


## unbiased analysis (diffusion map)
library(destiny)
css_pca_resi <- apply(Embeddings(seurat_endo,"css_pca"), 2, function(x) residuals(lm(x ~ seurat_endo$G2M.Score + seurat_endo$S.Score)))
seurat_endo[['css_pca_ccreg']] <- CreateDimReducObject(css_pca_resi, key = "CSSPCACCREG_")
dm_endo <- DiffusionMap(Embeddings(seurat_endo,"css_pca_ccreg"), n_pcs = NA)
dpt_endo <- DPT(dm_endo)
dpt_endo_final <- branch_divide(dpt_endo, integer(0L))[which(dpt_endo@branch[,1L]==min(dpt_endo@branch[,1L],na.rm=T) & dpt_endo@tips[,1L]),]
layout(matrix(1:2,nrow=1)); par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat_endo, "umap_css"), -(dpt_endo_final), colorPal = bluewhitered_colscheme, cex=1)
plotFeature(Embeddings(seurat_endo, "umap_css"), seurat_endo$AV, colorPal = bluewhitered_colscheme, cex=1)

seurat_endo_noprolif <- subset(seurat_endo, subset = RNA_snn_res.0.2 != "4")
dm_endo_noprolif <- DiffusionMap(Embeddings(seurat_endo_noprolif, "css_pca")[,1:20], n_pcs = NA)
dpt_endo_noprolif <- DPT(dm_endo_noprolif)
dpt_endo_noprolif_final <- branch_divide(dpt_endo_noprolif, integer(0L))[which(dpt_endo_noprolif@branch[,1L]==min(dpt_endo_noprolif@branch[,1L],na.rm=T) & dpt_endo_noprolif@tips[,1L]),]
layout(matrix(1:2,nrow=1)); par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat_endo_noprolif, "umap_css"), rank(dpt_endo_noprolif_final), colorPal = bluewhitered_colscheme, cex=1)
plotFeature(Embeddings(seurat_endo_noprolif, "umap_css"), seurat_endo_noprolif$AV, colorPal = bluewhitered_colscheme, cex=1, balance_col = F)

seurat_endo_noprolif$pt <- dpt_endo_noprolif_final
seurat_endo$pt <- max(dpt_endo_final) - dpt_endo_final
seurat_endo$pt_nonrepl <- seurat_endo_noprolif$pt[colnames(seurat_endo)]
saveRDS(seurat_endo, file="../data/data.seuratlanted_endothelial.rds")




# tip-stalk cells
tip_markers <- intersect(c("NRP1","DLL4","SOX17","KDR","CXCL12","CXCR4","APLN","ESM1","ANG2","PDGFB","KCNE3"), rownames(seurat_endo))
stalk_markers <- intersect(c("NR2F2","APLNR","JAG1","ETS1"), rownames(seurat_endo))




# the two arterious EC clusters
DE_c0_c1 <- wilcoxauc(seurat_endo, group_by = "RNA_snn_res.0.2", groups_use = c("0","1"))
DEG_c0_c1 <- lapply(split(DE_c0_c1,DE_c0_c1$group), function(x){
  x <- x[which(x$padj < 0.01 & x$logFC > log(1.2) & x$auc > 0.65 & x$pct_in - x$pct_out > 20 & x$pct_out < 20),]
  return(x[order(x$pct_out),])
})
scores <- sapply(DEG_c0_c1, function(x) colMeans(seurat_fetal_endo@assays$RNA@data[intersect(x$feature, rownames(seurat_fetal_endo)),]))
plotMultiFeatures(Embeddings(seurat_fetal_endo, "umap_css"), t(scores), colorPal=bluewhitered_colscheme, cex=1.2, random_order=F, sort_by_value=T)





# signaling (after running the regulome_signaling script)
load("res.signaling_pathway_activities_hBVO_invitro_transplanted.rdata")
seurat_endo[['signlaing_pc1']] <- CreateAssayObject(t(pc_signaling[colnames(seurat_endo),]))
seurat_invitro_endo[['signlaing_pc1']] <- CreateAssayObject(t(pc_signaling[colnames(seurat_invitro_endo),]))

avg_sig_invitro_endo <- sapply(as.character(c(18,11,5,6,20)), function(cl) rowMeans(seurat_invitro_endo@assays$signlaing_pc1@data[,which(seurat_invitro_endo$RNA_css_snn_res.0.8 == cl)]))
avg_sig_transpl_endo <- sapply(10:1, function(i) rowMeans(seurat_endo@assays$signlaing_pc1@data[,ceiling(rank(seurat_endo$AV)/ncol(seurat_endo)*10)==i]))
avg_sig_fetal_endo <- sapply(c("Arterial EC","Venous EC"), function(ct) rowMeans(seurat_fetal_endo@assays$signlaing_pc1@data[,which(seurat_fetal_endo$celltype == ct)]))
mat <- cbind(avg_sig_invitro_endo, NA, avg_sig_transpl_endo, NA, avg_sig_fetal_endo)
#mat <- mat[order(apply(mat,1,which.max)),]
mat_norm <- t(apply(mat, 1, function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))))
pdf("plot.heatmap_signaling_invitro_transplanted_primary.pdf")
heatmap.2(t(mat_norm), Rowv=NA, dendrogram="none", trace="none", key=F, keysize=0.2, col=bluered_colscheme(30), scale="none", margins=c(7,9), RowSideColors = rep(c("#bdbdbd",NA,"#696969",NA,"#303030"),c(5,1,10,1,2)), cexRow=1.3, cexCol=1, labCol = gsub("_", " ", names(m_signaling)))
dev.off()














# look at mesenchymal cells
seurat <- readRDS(file="../data/data.seuratlanted.rds")
seurat_mural <- subset(seurat, subset = RNA_css_snn_res.0.8 %in% c(5,11,8,2,0,17,14,15,1,4,13,12,7,3,10,16))
seurat_mural <- FindVariableFeatures(seurat_mural, nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(dims = 1:20) %>%
  cluster_sim_spectrum(label_tag="group", cluster_resolution = 1, use_fast_rank = T, dims_use = 1:20) %>%
  run_PCA(reduction="css", npcs = 20, reduction.name = "css_pca", reduction.key="CSSPCA_") %>%
  RunUMAP(reduction="css_pca", dims = 1:10, reduction.name="umap_css", reduction.key="UMAPCSS_")
seurat_mural_mnn <- RunFastMNN(SplitObject(seurat_mural, "group"))
seurat_mural[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mural_mnn,"mnn"), key="MNN_")
seurat_mural <- RunUMAP(seurat_mural, reduction="mnn", dims = 1:10, reduction.name="umap_mnn", reduction.key="UMAPMNN_")

seurat_mural <- FindNeighbors(seurat_mural, reduction = "css_pca", dims = 1:10) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1)

seurat_mural@active.ident <- seurat_mural$RNA_snn_res.0.5 # use the re-clustering result
(DimPlot(seurat_mural, reduction="umap_css", label=T)|DimPlot(seurat_mural, reduction="umap_mnn", label=T)) & NoAxes() & NoLegend()

saveRDS(seurat_mural, file="../data/data.seuratlanted_mural.rds")

seurat_mural@active.ident <- droplevels(seurat_mural$RNA_css_snn_res.0.8) # use the whole-dataset clustering result

DE_mural_cl <- wilcoxauc(seurat_mural)
DEG_mural_cl <- lapply(split(DE_mural_cl, DE_mural_cl$group), function(x){
  x <- x[which(x$padj < 0.01 & x$auc > 0.7 & x$logFC > log(1.2) & x$pct_in - x$pct_out > 20 & x$pct_out < 20),]
  return(x[order(x$pct_in - x$pct_out),])
})
markers_mural <- sort(unique(unlist(lapply(DEG_mural_cl, function(x) x$feature[1:min(c(5,nrow(x)))]))))

avg_expr_mural_cl <- sapply(levels(seurat_mural@active.ident), function(cl) rowMeans(seurat_mural@assays$RNA@data[,which(seurat_mural@active.ident == cl)]))
hcl_cl <- hclust(as.dist(1-cor(avg_expr_mural_cl[VariableFeatures(seurat_mural),])), method="ward.D2")
hcl_markers <- hclust(as.dist(1-cor(t(avg_expr_mural_cl[markers_mural,]))), method="ward.D2")
mat <- apply(avg_expr_mural_cl[markers_mural,],1,function(x) (x-min(x))/(max(x)-min(x)))
#mat <- mat[,order(apply(mat[rev(hcl_cl$order),],2,which.max))]
#gplots::heatmap.2(mat, Rowv=as.dendrogram(hcl_cl), Colv=NA, dendrogram="row", scale="none", trace="none", col = bluewhitered_colscheme(30), key=F, keysize=0.2, margins = c(8,3), cexCol = 1.2, cexRow = 1.5) # 1300x500

mat <- mat[intersect(names(col_cl_2),rownames(mat)),]
mat <- mat[,order(apply(mat,2,which.max))]
pdf("plot.transplanted_mesen_heatmap_cl_markers.pdf", width=10, height=5)
gplots::heatmap.2(mat, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", col = greyscale_colscheme(30), key=F, keysize=0.2, margins = c(8,3), cexCol = 1, cexRow = 1.5) # 1300x500
dev.off()

layout(matrix(1:2,nrow=1)); par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat_mural, "umap_css"), seurat_mural@active.ident, cex=0.5, pt_border = T, lwd_border = 0.1, do_label = T, label_round = T, cex.label = 1.2, colorPal = setNames(prettyrainbow_colscheme(length(levels(seurat_mural@active.ident))), levels(seurat_mural@active.ident)[hcl_cl$order]))
plotFeature(Embeddings(seurat_mural, "umap_css"), seurat_mural$line, cex=0.5, pt_border = T, lwd_border = 0.1, do_legend = T, legend_pos = "bottomleft", legend_cex = 1.2, colorPal = c("#303030","#cdcdcd"))

### compare to Cao et al. cell types
avg_expr_ct_atlas <- read.csv("~/Work/public_datasets/Cao_Science_2020_fetal_human_atlas_scRNAseq/normalized_expression/gene_expression_celltype.txt", row.names=1)
DE_ct_atlas <- read.csv("~/Work/public_datasets/Cao_Science_2020_fetal_human_atlas_scRNAseq/differential_expr_genes/DE_gene_77_main_cell_type.csv", quote = "\"'")
markers_atlas <- intersect(VariableFeatures(seurat_mural), unique(DE_ct_atlas$gene_short_name[which(DE_ct_atlas$padj < 0.01 & DE_ct_atlas$fold.change > 1.2)]))
df_genes_atlas <- readRDS("~/Work/public_datasets/Cao_Science_2020_fetal_human_atlas_scRNAseq/df_gene.RDS")
df_genes_atlas$gene_short_name <- as.character(df_genes_atlas$gene_short_name)
rownames(avg_expr_ct_atlas) <- make.unique(df_genes_atlas$gene_short_name)

sim_mural2atlas <- ref_sim_spectrum(avg_expr_mural_cl, as.matrix(avg_expr_ct_atlas)[markers_atlas,], scale=F)
#heatmap.2(sim_mural2atlas-min(sim_mural2atlas), Rowv=as.dendrogram(hcl_cl), scale="none", trace="none", col=bluewhitered_colscheme(30), cexCol = 0.8, cexRow = 1.5, margins = c(15,3), key=F, keysize=0.2)
#heatmap.2((sim_mural2atlas-min(sim_mural2atlas))[,which(apply(sim_mural2atlas,2,max)>0.5)], Rowv=as.dendrogram(hcl_cl), scale="none", trace="none", col=bluewhitered_colscheme(30), cexCol = 1.2, cexRow = 1.5, margins = c(20,3), key=F, keysize=0.2)
heatmap.2((sim_mural2atlas-min(sim_mural2atlas))[intersect(names(col_cl_2),rownames(mat)),which(apply(sim_mural2atlas,2,max)>0.5)],
          Rowv=NA, dendrogram="col", scale="none", trace="none", col=bluewhitered_colscheme(30), cexCol = 1.2, cexRow = 1.2, margins = c(20,3), key=F, keysize=0.2)

avg_expr_tissue_atlas <- read.csv("~/Work/public_datasets/Cao_Science_2020_fetal_human_atlas_scRNAseq/normalized_expression/gene_expression_tissue.txt", row.names=1)
markers_tissue_atlas <- read.csv("~/Work/public_datasets/Cao_Science_2020_fetal_human_atlas_scRNAseq/differential_expr_genes/DE_gene_by_organ.csv", quote = "\"'")

### compare to Yu et al. mesenchymal subtypes
seurat_gut_mesen <- readRDS("~/Work/public_datasets/Yu_Cell_2021_guttube_atlas/multiorgan_mesenchyme.seurat.rds")

avg_expr_gut_mesen <- sapply(sort(unique(seurat_gut_mesen$Cell_type)), function(ct) rowMeans(seurat_gut_mesen@assays$RNA@data[,which(seurat_gut_mesen$Cell_type == ct)]))
hcl_gut_mesen_ct <- hclust(as.dist(1-cor(avg_expr_gut_mesen[VariableFeatures(seurat_gut_mesen),])), method="ward.D2")
DE_mesen_ct <- read.table("~/Work/public_datasets/Yu_Cell_2021_guttube_atlas/markers_mesenchyme.tsv", header=T, sep="\t")
DE_mesen_ct <- lapply(split(DE_mesen_ct, DE_mesen_ct$group), function(x){
  x <- x[which(x$logFC > log(1.2) & x$auc > 0.7 & x$pct_in - x$pct_out > 20),]
  return(x)
})

markers_mesen_ct <- sort(unique(unlist(lapply(DE_mesen_ct, function(x) x$feature[1:min(c(20,nrow(x)))]))))
sim_mural2gut <- ref_sim_spectrum(avg_expr_mural_cl, avg_expr_gut_mesen[markers_mesen_ct,], scale=F)
heatmap.2(sim_mural2gut-min(sim_mural2gut), Rowv=as.dendrogram(hcl_cl), Colv=as.dendrogram(hcl_gut_mesen_ct), scale="none", trace="none", col=bluewhitered_colscheme(30), cexCol = 1.2, cexRow = 1.5, margins = c(15,3), key=F, keysize=0.5)
pdf("plot.transplanted_mesen_heatmap_sim2gut.pdf", width=10, height=7)
heatmap.2((sim_mural2gut-min(sim_mural2gut))[intersect(names(col_cl_2),rownames(sim_mural2gut)),],
          Rowv=NA, Colv=as.dendrogram(hcl_gut_mesen_ct), dendrogram="none",scale="none", trace="none", col=bluewhitered_colscheme(30), cexCol = 1.2, cexRow = 1.5, margins = c(15,3), key=F, keysize=0.2)
dev.off()

markers_mesen_ct <- sort(unique(unlist(lapply(DE_mesen_ct, function(x) x$feature[1:min(c(10,nrow(x)))]))))
mat <- apply(avg_expr_gut_mesen[markers_mesen_ct,],1,function(x) (x-min(x))/(max(x)-min(x)))
order_markers <- order(apply(mat[rev(hcl_gut_mesen_ct$order),],2,which.max))
heatmap.2(mat[,order_markers], Rowv=as.dendrogram(hcl_gut_mesen_ct), Colv=NA, dendrogram="row", scale="none", trace="none", col=bluewhitered_colscheme(30), cexCol = 0.7, cexRow = 1, margins = c(5,11), key=F, keysize=0.3)

markers_mesen_ct <- intersect(sort(unique(unlist(lapply(DE_mesen_ct, function(x) x$feature[1:min(c(10,nrow(x)))])))), rownames(seurat_mural))
mat <- apply(cbind(avg_expr_gut_mesen[markers_mesen_ct,],NA,avg_expr_mural_cl[markers_mesen_ct,]),1,function(x) (x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)))
order_markers <- order(apply(mat[rev(hcl_gut_mesen_ct$order),],2,which.max))
heatmap.2(mat[c(rev(hcl_gut_mesen_ct$order),c(0,rev(hcl_cl$order))+ncol(avg_expr_gut_mesen)+1),order_markers],
          Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", col=bluewhitered_colscheme(30), cexCol = 0.7, cexRow = 1, margins = c(5,11), key=F, keysize=0.3)
