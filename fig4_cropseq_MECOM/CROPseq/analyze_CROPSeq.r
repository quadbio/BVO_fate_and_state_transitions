library(Seurat)
library(simspec)
library(dplyr)
library(ggplot2)
library(Matrix)

seurat <- readRDS("seurat.rds")
seurat <- subset(seurat, subset = sample != "I_J_d14")

pdf("plot.guide_number_batch.pdf", height=4, width=4)
layout(matrix(1:2, nrow=1, byrow=T)); par(mar=c(5,4,3,1))
for(batch in unique(seurat$batch)){
  guide_num_per_cell <- colSums(seurat@assays$guide@counts[,which(seurat$batch == batch)] > 0)
  guide_num_per_cell <- sapply(0:max(c(min_max_num,guide_num_per_cell)), function(i) sum(guide_num_per_cell == i))
  barplot(guide_num_per_cell/100, col = "#efefef", xlab="# of guides", ylab="Cell number (hundred)", main=batch, cex.lab=1.2, cex.main=1.5, names.arg=0:(length(guide_num_per_cell)-1))
}
dev.off()

# data processing and integration
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 3000)
blacklist <- c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt", stringsAsFactors=F)[,1])
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
seurat <- CellCycleScoring(seurat, s.features=cc.genes.updated.2019$s.genes, g2m.features=cc.genes.updated.2019$g2m.genes, set.ident=F)
seurat <- ScaleData(seurat, vars.to.regress = c("G2M.Score","S.Score")) %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

seurat <- cluster_sim_spectrum(seurat, label_tag="sample") %>%
  run_PCA(reduction="css", npcs = 20, reduction.name = "css_pca", reduction.key = "CSSPCA_") %>%
  RunUMAP(reduction="css_pca", dims = 1:10, reduction.name = "umap_css", reduction.key = "UMAPCSS_")

seurat <- FindNeighbors(seurat, "css_pca", dims = 1:10) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1) %>%
  FindClusters(resolution = 2)
seurat$RNA_css_snn_res.0.1 <- seurat$RNA_snn_res.0.1
seurat$RNA_css_snn_res.0.5 <- seurat$RNA_snn_res.0.5
seurat$RNA_css_snn_res.1 <- seurat$RNA_snn_res.1
seurat$RNA_css_snn_res.2 <- seurat$RNA_snn_res.2


# cluster annotation
seurat$major_ct <- "mural"
seurat$major_ct[which(seurat$RNA_css_snn_res.2 %in% c(2,5,23,27))] <- "endothelial"
seurat$major_ct[which(seurat$RNA_css_snn_res.2 %in% c(31,33))] <- "blood"
seurat$major_ct[which(seurat$RNA_css_snn_res.2 %in% c(28,17))] <- "other-1"
seurat$major_ct[which(seurat$RNA_css_snn_res.2 %in% c(21))] <- "other-2"
seurat$major_ct[which(seurat$RNA_css_snn_res.2 %in% c(22,24))] <- "other-3"
seurat$major_ct <- factor(seurat$major_ct)



# mural and endothelial heterogeneity
seurat_mc <- subset(seurat, subset = major_ct == "mural")
seurat_ec <- subset(seurat, subset = major_ct == "endothelial")

### mural cells
seurat_mc <- FindVariableFeatures(seurat_mc, nfeatures = 5000)
blacklist <- c(grep("^MT-", rownames(seurat), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt", stringsAsFactors=F)[,1])
VariableFeatures(seurat_mc) <- setdiff(VariableFeatures(seurat_mc), blacklist)
seurat_mc <- ScaleData(seurat_mc) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
(DimPlot(seurat_mc, reduction="umap", group.by = "Phase") | DimPlot(seurat_mc, reduction="umap", group.by = "sample")) & NoAxes()

seurat_mc <- cluster_sim_spectrum(seurat_mc, label_tag="sample", cluster_resolution = 0.6, use_fast_rank=T) %>%
  run_PCA(reduction="css", reduction.name="css_pca", reduction.key="CSSPCA_", npcs=20)
seurat_mc <- RunUMAP(seurat_mc, reduction = "css_pca", dims = 1:10, n.neighbors = 50, reduction.name = "umap_css", reduction.key = "UMAPCSS_")

seurat_mc <- FindNeighbors(seurat_mc, reduction = "css_pca", dims = 1:10) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1)
seurat_mc$ct <- factor(setNames(c(rep("mural_mesenchyme",1),
                                  rep("mural_TCF21+",4),
                                  rep("mural_VEGFD+",3),
                                  rep("mural_pericyte",2),
                                  rep("mural_SM",3)),
                                c(c(5),c(1,11,3,0),c(7,4,10),c(8,2),c(9,6)))[as.character(seurat_mc$RNA_snn_res.0.5)])


### endothelial cells
seurat_ec <- FindVariableFeatures(seurat_ec, nfeatures = 5000)
blacklist <- c(grep("^MT-", rownames(seurat_ec), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt", stringsAsFactors=F)[,1])
VariableFeatures(seurat_ec) <- setdiff(VariableFeatures(seurat_ec), blacklist)
seurat_ec <- ScaleData(seurat_ec) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20) %>%
  cluster_sim_spectrum(label_tag="sample", cluster_resolution = 0.6, use_fast_rank=T) %>%
  run_PCA(reduction="css", reduction.name="css_pca", reduction.key="CSSPCA_", npcs=20) %>%
  RunUMAP(reduction = "css_pca", dims = 1:10, n.neighbors = 50, reduction.name = "umap_css", reduction.key = "UMAPCSS_")

seurat_ec <- FindNeighbors(seurat_ec, "css_pca", dims = 1:10) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.2) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1)

seurat$ct <- as.character(seurat$major_ct)
seurat$ct[colnames(seurat_mc)] <- as.character(seurat_mc$ct)
seurat$ct[colnames(seurat_ec)] <- paste0("endothelial-", as.numeric(seurat_ec$RNA_snn_res.0.2)+1)





# KO effect: finer resolution composition change of example genes (ETV2, MECOM)
idx1 <- which(seurat$targets=="ETV2")
idx2 <- which(seurat$targets=="MECOM")
idx3 <- which(seurat$targets=="DUMMY")
p1 <- ggplot(data.frame(Embeddings(seurat,"umap_css")[idx1,]), aes(x=UMAPCSS_1, y=UMAPCSS_2) ) +
  geom_hex(bins = 20) +
  scale_fill_gradientn(colours = greyscale_colscheme(30)) + theme_void() + theme(legend.position = "none", plot.title = element_blank())
p2 <- ggplot(data.frame(Embeddings(seurat,"umap_css")[idx2,]), aes(x=UMAPCSS_1, y=UMAPCSS_2) ) +
  geom_hex(bins = 20) +
  scale_fill_gradientn(colours = greyscale_colscheme(30)) + theme_void() + theme(legend.position = "none", plot.title = element_blank())
p3 <- ggplot(data.frame(Embeddings(seurat,"umap_css")[idx3,]), aes(x=UMAPCSS_1, y=UMAPCSS_2) ) +
  geom_hex(bins = 20) +
  scale_fill_gradientn(colours = greyscale_colscheme(30)) + theme_void() + theme(legend.position = "none", plot.title = element_blank())
p1 | p2 | p3



# KO effect: DE analysis on endothelial cells and mural cells separately (no subtype considered)
detection_rates <- rowMeans(seurat@assays$RNA@data > 0)
gene_cands <- names(which(detection_rates > 0.01))
blacklist <- c(grep("^MT-",rownames(seurat),value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1], grep("-gene$",rownames(seurat),value=T))
gene_cands <- setdiff(gene_cands, blacklist)

DE_ct_prob <- setNames(lapply(c("mural","endothelial"), function(ct){
  setNames(lapply(sort(setdiff(rownames(seurat@assays$target), "DUMMY")), function(target){
    idx_single_target <- which(colSums(seurat@assays$target@data>0)==1 & seurat$major_ct == ct)
    idx_target <- intersect(which(seurat@assays$target@data[target,] > 0), idx_single_target)
    idx_bg <- intersect(which(seurat@assays$target@data["DUMMY",] > 0), idx_single_target)
    idx <- c(idx_target, idx_bg)
    
    DE_overall_prob <- ancova_group_test(seurat@assays$RNA@data[gene_cands,idx],
                                         group = seurat$perturb_prob[idx],
                                         covar = data.frame(sample = seurat$sample[idx]),
                                         num_threads = 50)
    DE_overall_prob <- data.frame(target = target, DE_overall_prob)
    message(paste0("...done DE for ",target, " in ", ct))
    return(DE_overall_prob)
  }), sort(setdiff(rownames(seurat@assays$target), "DUMMY")))
}), c("mural","endothelial"))
saveRDS(DE_ct_prob, file="res.DE-anova_perturb-prop_target_vs_dummy_ct.rds")

df_DE_ct <- do.call(rbind, lapply(names(DE_ct_prob), function(ct)
  do.call(rbind, lapply(DE_ct_prob[[ct]], function(x)
    data.frame(celltype = ct, gene = rownames(x), x[,c(1:5,8)], padj = p.adjust(x$p_ANOVA, method="bonferroni"))))))
rownames(df_DE_ct) <- NULL

num_DEGs <- sapply(split(df_DE_ct,df_DE_ct$celltype), function(x) tapply(x$padj < 0.1, x$target, sum, na.rm=T))
pooled_DEGs <- lapply(split(df_DE_ct,df_DE_ct$target), function(x) unique(unlist(lapply(split(x,x$celltype), function(x) x$gene[which(x$padj<0.1)]))))
summary_DEGs_ct <- data.frame(num = lengths(pooled_DEGs),
                              cor = sapply(names(pooled_DEGs), function(x) cor(df_DE_ct$coef[df_DE_ct$target==x & df_DE_ct$celltype=="mural" & df_DE_ct$gene %in% pooled_DEGs[[x]]],
                                                                               df_DE_ct$coef[df_DE_ct$target==x & df_DE_ct$celltype=="endothelial" & df_DE_ct$gene %in% pooled_DEGs[[x]]], use="complete.obs")),
                              p = sapply(names(pooled_DEGs), function(x){
                                test = try(cor.test(df_DE_ct$coef[df_DE_ct$target==x & df_DE_ct$celltype=="mural" & df_DE_ct$gene %in% pooled_DEGs[[x]]], df_DE_ct$coef[df_DE_ct$target==x & df_DE_ct$celltype=="endothelial" & df_DE_ct$gene %in% pooled_DEGs[[x]]], use="complete.obs"))
                                return(ifelse(sum(class(test)=="try-error")>0, 1, test$p.value))
                              }))

layout(matrix(1:2,nrow=1)); par(mar=c(5,5,1,1))
plot(num_DEGs, bty="n", cex=1.5, pch=16, xlab="# DEGs in endothelial cells", ylab="# DEGs in mural cells", cex.lab=1.5)
abline(v=20, h=30, lty=2)
text(x=num_DEGs[which(num_DEGs[,1]>20 | num_DEGs[,2]>30),1], y=num_DEGs[which(num_DEGs[,1]>20 | num_DEGs[,2]>30),2]-4,
     labels = names(which(num_DEGs[,1]>20 | num_DEGs[,2]>30)), cex=1.2)
plot(summary_DEGs_ct[,2],-log10(p.adjust(summary_DEGs_ct[,3],method="BH")), pch=16, ylab="-Log10(p-adj)", xlab="cor(Endo. vs. Mural)", cex=1.5, cex.lab=1.5, bty="n", xlim=c(-1,1))
abline(v=0, lty=2)
text(x=summary_DEGs_ct[which(p.adjust(summary_DEGs_ct[,3],method="BH")<0.1),2],
     y=-log10(p.adjust(summary_DEGs_ct[,3],method="BH")[which(p.adjust(summary_DEGs_ct[,3],method="BH")<0.1)])-0.1,
     labels = rownames(summary_DEGs_ct)[which(p.adjust(summary_DEGs_ct[,3],method="BH")<0.1)], cex=1.2)


# DEG characterization for MECOM in EC
DE_EC_MECOM <- df_DE_ct %>% filter(celltype=="endothelial" & target == "MECOM")
plot(DE_EC_MECOM$coef, -log10(DE_EC_MECOM$p_ANOVA), pch=16, col=ifelse(DE_EC_MECOM$padj<0.1, "#303030", "#cdcdcd"), bty="n", xlab="DE in EC", ylab="-log10(P)", cex.lab=1.5)
abline(h = -log10(max(DE_EC_MECOM$p_ANOVA[DE_EC_MECOM$padj<0.1], na.rm=T)), v=0, lty=2)

