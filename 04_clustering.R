###############################################################################
# 04_clustering.R
#
#   1. Major clusters (res 0.47) and marker genes
#   2. Major cell-type annotation (cluster labels + CITE-seq ADT -> PBMC_identity)
#   3. Subclusters (res 0.8) -> published subcluster names (CD4-1 ... pDC)
#   4. Survival metadata
#   5. Figures and DE:
#        Fig 2A   UMAP of pre-ACT PBMC cell types, and coloured by OS
#        Fig 3A   UMAP of pre-ACT subclusters
#        Fig 3C   Heatmap of DEGs in subclusters dominating M2/M6 expression
#        Fig 4C   Hallmark pathway enrichment of CD8-2 DEGs
#        Fig 7B   Heatmap of DEGs across CD8 subclusters
#        Supp 3A  Heatmap of DEGs across all subclusters
#
# Input : results/rds/03_integrated.rds, data/metadata/patient_survival.csv
# Output: results/rds/04_clustered.rds (+ tables and figures)
###############################################################################

source("R/helpers.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(enrichR)
})
set.seed(42)

obj <- readRDS(rds_path("03_integrated.rds"))
DefaultAssay(obj) <- "SCT"

# =============================================================================
# 1. Major clusters
# =============================================================================
obj <- FindClusters(obj, resolution = 0.47, graph.name = "harmony_snn",
                    cluster.name = "harmony_clusters")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

markers <- FindAllMarkers(PrepSCTFindMarkers(obj, assay = "SCT"), only.pos = TRUE,
                          assay = "SCT", slot = "data") %>%
  filter(avg_log2FC > 0.5, !grepl("^RPS|^RPL|^MT-", gene))
write.csv(markers, tab_path("markers_major_clusters.csv"), row.names = FALSE)

# =============================================================================
# 2. Major cell-type annotation
# -----------------------------------------------------------------------------
# Detailed labels of the res-0.47 clusters (from marker genes + CITE-seq ADT)
# =============================================================================
cluster_map <- c(
  "0" = "Naive T",        "1" = "CD4+ TCM",        "2" = "Cytotoxic CD8+",
  "3" = "GzmK+ CD8+ TEM", "4" = "Quiescent CD4+ TCM", "5" = "CD14+ Monocyte",
  "6" = "NK",             "7" = "Activated NK/CD8+", "8" = "B",
  "9" = "CD16+ Monocyte", "10" = "Plasmablast"
)
obj$celltype <- unname(cluster_map[as.character(obj$harmony_clusters)])

# -----------------------------------------------------------------------------
# Broad cell types used throughout the paper.
# B cells, monocytes and NK cells are taken from the cluster labels above;
# the remaining (T-cell cluster) cells are split into CD4 / CD8 using the
# CITE-seq ADT signal (CLR-normalised):
#   T cell : ADT CD3 >= 0.25
#   CD8    : T cell and ADT CD8 > 1 and ADT CD4 < 1
#   CD4    : T cell and ADT CD4 > 1 and ADT CD8 < 1
# Cells that are neither (e.g. CD3-low or CD4/CD8 double-positive/negative
# cells in T-cell clusters) were labelled NK, as in the analysis for the paper.
#
#   PBMC_identity   : CD4, CD8, NK, Monocyte, B cell   (pDC added in section 3)
#   PBMC_identity_2 : same, but CD4 + CD8 merged into "CD3" (used in Fig 6C)
# -----------------------------------------------------------------------------
adt <- GetAssayData(obj, assay = "ADT", layer = "data")
obj$ADT_CD3 <- as.numeric(adt["CD3", ])
obj$ADT_CD4 <- as.numeric(adt["CD4.1", ])
obj$ADT_CD8 <- as.numeric(adt["CD8", ])

is_T <- obj$ADT_CD3 >= 0.25
obj$T_call <- "non-T"
obj$T_call[is_T & obj$ADT_CD8 > 1 & obj$ADT_CD4 < 1] <- "CD8"
obj$T_call[is_T & obj$ADT_CD8 < 1 & obj$ADT_CD4 > 1] <- "CD4"
obj$T_call <- factor(obj$T_call, levels = c("CD8", "CD4", "non-T"))
print(table(obj$T_call, useNA = "ifany"))

obj$PBMC_identity <- case_when(
  obj$celltype %in% c("B", "Plasmablast")                 ~ "B cell",
  obj$celltype %in% c("CD14+ Monocyte", "CD16+ Monocyte") ~ "Monocyte",
  obj$celltype %in% "NK"                                  ~ "NK",
  obj$T_call == "CD8"                                     ~ "CD8",
  obj$T_call == "CD4"                                     ~ "CD4",
  TRUE                                                    ~ "NK"
)
obj$PBMC_identity_2 <- ifelse(obj$PBMC_identity %in% c("CD4", "CD8"), "CD3", obj$PBMC_identity)
print(table(obj$celltype, obj$PBMC_identity, useNA = "ifany"))

p <- FeaturePlot(obj, features = c("CD19.1", "CD11c", "CD14.1", "CD16", "CD3", "CD4.1", "CD56", "CD8"),
                 reduction = "umap.harmony", cols = c("lightgrey", "#2c7fb8"), ncol = 4, order = TRUE)
save_plot(p, "ADT_markers_UMAP.png", width = 16, height = 8, dpi = 200)

# =============================================================================
# 3. Subclusters (res 0.8) and published names
# =============================================================================
obj <- FindClusters(obj, resolution = 0.8, graph.name = "harmony_snn",
                    cluster.name = "harmony_clusters2")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

new_labels <- setNames(subcluster_map$new, subcluster_map$old)
stopifnot(length(setdiff(unique(as.character(obj$harmony_clusters2)), names(new_labels))) == 0)
obj$subcluster_updated <- factor(unname(new_labels[as.character(obj$harmony_clusters2)]),
                                 levels = subcluster_order)
print(table(obj$harmony_clusters2, obj$subcluster_updated))

# Subcluster 17 expresses pDC markers (PLD4, SERPINF1, ...): relabel from B cell
obj$PBMC_identity <- as.character(obj$PBMC_identity)
obj$PBMC_identity[as.character(obj$harmony_clusters2) == "17"] <- "pDC"

# =============================================================================
# 4. Survival metadata
# =============================================================================
obj <- add_survival_metadata(obj)
stopifnot(sum(is.na(obj$OS_from_diagnosis_months)) == 0)

pre <- subset(obj, subset = timepoint == "Pre")

# =============================================================================
# 5a. Fig 2A - pre-ACT PBMC UMAPs
# =============================================================================
p <- DimPlot(pre, reduction = "umap.harmony", group.by = "PBMC_identity",
             cols = pbmc_cols, pt.size = 0.1, alpha = 0.4, raster = FALSE) +
  ggtitle(NULL)
save_plot(p, "Fig2A_UMAP_celltypes_pre.png", width = 10, height = 8)

p <- DimPlot(pre, reduction = "umap.harmony", group.by = "patient_OS_label",
             cols = patient_os_cols(pre$patient_OS_label), pt.size = 0.1, alpha = 0.5,
             shuffle = TRUE, raster = FALSE) +
  ggtitle("Patient (ordered by OS from diagnosis)")
save_plot(p, "Fig2A_UMAP_patient_OS_pre.png", width = 8, height = 7)

# =============================================================================
# 5b. Fig 3A - pre-ACT subcluster UMAP
# =============================================================================
p <- DimPlot(pre, reduction = "umap.harmony", group.by = "subcluster_updated",
             cols = subcluster_cols, label = TRUE, label.box = TRUE, repel = TRUE,
             pt.size = 0.1, alpha = 0.3, raster = FALSE) + ggtitle(NULL)
save_plot(p, "Fig3A_UMAP_subclusters_pre.png", width = 8, height = 7)

# =============================================================================
# 5c. Subcluster markers + Supp 3A heatmap
# =============================================================================
Idents(obj) <- "subcluster_updated"
obj_sct <- PrepSCTFindMarkers(obj, assay = "SCT")

sub_markers <- FindAllMarkers(obj_sct, only.pos = TRUE, assay = "SCT", slot = "data") %>%
  filter(avg_log2FC > 0.5, !grepl("^RPS|^RPL|^MT-", gene))
write.csv(sub_markers, tab_path("markers_subclusters.csv"), row.names = FALSE)

supp3a_genes <- c(
  "CCR7","LEF1","TCF7","ITGA6","BCL2","TXNIP","CD27","IL6ST","AQP3","FOS",
  "CCR4","CCR6","CCR10","RGCC","PABPC1","FTH1","DUSP4","FOXP3","HIVEP2","FOXP1",
  "KLF12","CCL5","GZMK","CXCR3","GNLY","GZMH","CCL4","IFNG","TNF","TYROBP",
  "NKG7","PTGDS","FCGR3A","S100A8","CST3","LYZ","HLA-DMA","IFI30","SPI1","SERPINF1",
  "PLD4","CD79A","HLA-DRA","ITM2C","IGHA1"
)
avg <- avg_expression_table(obj, supp3a_genes, "subcluster_updated", assay = "RNA")
write.csv(as.data.frame(avg) %>% rownames_to_column("Gene"),
          tab_path("SuppFig3A_avg_expression_all_subclusters.csv"), row.names = FALSE)
plot_avg_heatmap(avg, "SuppFig3A_heatmap_subclusters.png", width = 8, height = 10)

# =============================================================================
# 5d. Fig 3C - DEGs of the subclusters dominating M2/M6 expression
# (subclusters that together account for >=50% of module expression; see
#  05_hdWGCNA.R for the module scores)
# =============================================================================
fig3c_clusters <- c("CD4-1", "CD4-2", "CD4-3", "CD8-2", "NK-1", "NK-2")
fig3c_genes <- c(
  "CD55","CCR7","LTB","IL7R","LEF1","ITGA6","IL6ST","CD27","PLCL1","GIMAP7",
  "LMNA","AQP3","CD8B","CCL5","GZMH","GNLY","CTSW","S1PR5","NKG7",
  "TYROBP","GZMB","CD7","PTGDS","CCL4"
)
avg <- avg_expression_table(obj, fig3c_genes, "subcluster_updated", assay = "RNA",
                            groups = fig3c_clusters)
write.csv(as.data.frame(avg) %>% rownames_to_column("Gene"),
          tab_path("Fig3C_avg_expression_dominant_subclusters.csv"), row.names = FALSE)
plot_avg_heatmap(avg, "Fig3C_heatmap_dominant_subclusters.png", width = 5, height = 8)

# =============================================================================
# 5e. Fig 7B - CD8 subclusters
# =============================================================================
cd8_pairs <- list(c("CD8-2", "CD8-3"), c("CD8-2", "CD8-1"), c("CD8-3", "CD8-1"))
for (pr in cd8_pairs) {
  m <- FindMarkers(obj, ident.1 = pr[1], ident.2 = pr[2], assay = "RNA", slot = "data",
                   logfc.threshold = 0.25, min.pct = 0.25) %>%
    rownames_to_column("gene") %>%
    arrange(p_val_adj, desc(avg_log2FC))
  write.csv(m, tab_path(paste0("Fig7B_DE_", pr[1], "_vs_", pr[2], ".csv")), row.names = FALSE)
}

fig7b_genes <- c("FGFBP2", "GNLY", "LGALS1", "GZMH", "GZMK", "CXCR3",
                 "LTB", "IL7R", "CCL4L2", "CCL4", "IFNG", "TNF")
avg <- avg_expression_table(obj, fig7b_genes, "subcluster_updated", assay = "RNA",
                            groups = c("CD8-1", "CD8-2", "CD8-3"))
write.csv(as.data.frame(avg) %>% rownames_to_column("Gene"),
          tab_path("Fig7B_avg_expression_CD8_subclusters.csv"), row.names = FALSE)
plot_avg_heatmap(avg, "Fig7B_heatmap_CD8_subclusters.png", width = 4, height = 6)

# =============================================================================
# 5f. Fig 4C - CD8-2 DEGs and Hallmark enrichment
# -----------------------------------------------------------------------------
# NOTE: the code for Fig 4B/4C was not in the original scripts; this block
# follows the figure legend (Padj < 0.05, |log2FC| > 0.5, MSigDB Hallmark,
# FDR < 0.05). Check the output against the published panel.
# =============================================================================
cd8_2_de <- FindMarkers(obj_sct, ident.1 = "CD8-2", assay = "SCT", slot = "data",
                        logfc.threshold = 0.5) %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5, !grepl("^RPS|^RPL|^MT-", gene)) %>%
  arrange(p_val_adj)
write.csv(cd8_2_de, tab_path("Fig4B_CD8-2_DEGs.csv"), row.names = FALSE)

hallmark <- enrichr(cd8_2_de$gene, "MSigDB_Hallmark_2020")[["MSigDB_Hallmark_2020"]] %>%
  filter(Adjusted.P.value < 0.05) %>%
  mutate(gene_count = as.numeric(sub("/.*", "", Overlap)),
         neg_log10_fdr = -log10(Adjusted.P.value),
         Term = reorder(Term, neg_log10_fdr))
write.csv(hallmark, tab_path("Fig4C_CD8-2_Hallmark_enrichment.csv"), row.names = FALSE)

p <- ggplot(hallmark, aes(x = neg_log10_fdr, y = Term)) +
  geom_segment(aes(x = 0, xend = neg_log10_fdr, yend = Term), colour = "grey60") +
  geom_point(aes(size = gene_count, colour = neg_log10_fdr)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "red") +
  scale_colour_gradient(low = "#6CC3F4", high = "#B2182B") +
  labs(x = expression(-log[10](P[FDR])), y = NULL,
       size = "CD8-2 DEGs", colour = expression(-log[10](P[FDR]))) +
  theme_classic(base_size = 13)
save_plot(p, "Fig4C_CD8-2_Hallmark_lollipop.png", width = 7, height = 5)

saveRDS(obj, rds_path("04_clustered.rds"))
