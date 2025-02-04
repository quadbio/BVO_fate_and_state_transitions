library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(Seurat)
library(tidyverse)
library(RCurl)
library(ggplot2)
library(SeuratWrappers)
library(ggrepel)
library(uwot)

setwd("/home/marinani/PhD_Projects/Vascular_Organoids/Analysis/MECOM_KO/MECOM_KO_day7_day21_together/")
seurat <- readRDS("/home/marinani/PhD_Projects/Vascular_Organoids/Analysis/MECOM_KO/MECOM_KO_day7_day21_together/Objects/MECOM_KO_css_integrated_seurat.rds")

seurat_sce <- as.SingleCellExperiment(seurat)
seurat_milo <- Milo(seurat_sce)

seurat_milo <- buildGraph(seurat_milo, k = 30, d = 30)


seurat_milo <- makeNhoods(seurat_milo, prop = 0.2, k = 30, d=30, refined = TRUE,refinement_scheme="graph")

plotNhoodSizeHist(seurat_milo)

seurat_milo <- countCells(seurat_milo, meta.data = data.frame(colData(seurat_milo)), sample="orig.ident")
head(nhoodCounts(seurat_milo))

seurat_design <- data.frame(colData(seurat_milo))[,c("orig.ident", "condition", "sample","timepoint", "replicate","KO_replicate")]

seurat_design$replicate <- as.factor(seurat_design$replicate)
seurat_design <- distinct(seurat_design)
rownames(seurat_design) <- seurat_design$orig.ident

saveRDS(seurat_milo, "/home/marinani/PhD_Projects/Vascular_Organoids/Analysis/MECOM_KO/MECOM_KO_day7_day21_together/Objects/seurat_milo_k30_prop0.2.rds")


seurat_design <- data.frame(colData(seurat_milo))[,c("orig.ident", "condition", "sample","timepoint", "replicate","KO_replicate")]

seurat_design$condition <- factor(seurat_design$condition, levels = c("WT","KO"))
da_results <- testNhoods(seurat_milo, design = ~ timepoint + condition, design.df = seurat_design, fdr.weighting="graph-overlap")

da_results %>%
  arrange(- SpatialFDR) %>%
  head()

pdf("Plots/NEW_histo_da_results_miloR_timepoint_condition_k30_prop0.2.pdf", width = 10, height = 7)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50) + theme_classic()
dev.off()

pdf("Plots/NEW_volcano_da_results_miloR_timepoint_condition_k30_prop0.2.pdf", width = 5, height = 3.5)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point(colour = "grey", alpha = 0.5, size = 0.8) +
  geom_point(data = filter(da_results, -log10(SpatialFDR) > 1), size = 1.2,shape = 21,fill = "black",colour = "black") +
  theme_classic()
  # geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
dev.off()


seurat_milo <- buildNhoodGraph(seurat_milo)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(seurat_milo, da_results, layout="UMAP_CSS",alpha=0.1)


da_results <- annotateNhoods(seurat_milo, da_results, coldata_col = "RNA_css_snn_res.0.4")
head(da_results)

# While neighbourhoods tend to be homogeneous, we can define a threshold for celltype_fraction to exclude neighbourhoods that are a mix of cell types.
ggplot(da_results, aes(RNA_css_snn_res.0.4_fraction)) + geom_histogram(bins=50)

da_results$RNA_css_snn_res.0.4 <- ifelse(da_results$RNA_css_snn_res.0.4_fraction < 0.5, "Overlapping", da_results$RNA_css_snn_res.0.4)

pdf("Plots/NEW_plotDAbeeswarm_miloR_timepoint_condition_k30_prop0.2.pdf", width = 4, height = 8)
plotDAbeeswarm(da_results, group.by = "RNA_css_snn_res.0.4") + theme_classic()
dev.off()

pdf("Plots/NEW_umap_miloR_timepoint_condition_k30_prop0.2.pdf", width = 8, height = 4.5)
nh_graph_pl + plot_layout(guides="collect")
dev.off()

pdf("Plots/NEW_umap_miloR_by_neighborhood_k30_prop0.2.pdf", width = 6, height = 4.5)
plotNhoodGraph(seurat_milo, colour_by =  "RNA_css_snn_res.0.4", layout="UMAP_CSS")
dev.off()

png("Plots/NEW_umap_miloR_timepoint_condition_k30_prop0.2.png", width = 2700, height = 1750, res = 500)
nh_graph_pl + plot_layout(guides="collect")
dev.off()

png("Plots/NEW_umap_miloR_by_neighborhood_k30_prop0.2.png", width = 2900, height = 2250, res = 500)
plotNhoodGraph(seurat_milo, colour_by =  "RNA_css_snn_res.0.4", layout="UMAP_CSS")
dev.off()

saveRDS(da_results, "/home/marinani/PhD_Projects/Vascular_Organoids/Analysis/MECOM_KO/MECOM_KO_day7_day21_together/Plots/seurat_milo_da_results_k30_prop0.2.rds")