###############################################################################
# 05_hdWGCNA.R
#
# Weighted gene co-expression network analysis (hdWGCNA) of pre-ACT PBMCs.
#   1. Metacells (per cell type x patient), soft power, network, modules
#   2. Module eigengenes (harmonised by patient), kME, module names PBMC-M1..M8
#   3. Differential module eigengenes, LTS vs STS (labels in Supp Fig 1C)
#   4. Module scores for every cell (pre + post) from genes with kME >= 0.5
#   5. Figures:
#        Fig 2D   UMAP of M2 and M6 module eigengenes (pre-ACT)
#        Fig 2E   Module gene UMAP of M2 and M6 hub genes
#        Fig 3B   Dot plot of M2/M6 expression across subclusters
#        Fig 4B   Venn diagram: CD8-2 DEGs vs M2 hub genes
#        Supp 1C  Module dendrogram
#        Supp 2   GO enrichment of M2 and M6 (EnrichR)
#
# Input : results/rds/04_clustered.rds, results/tables/Fig4B_CD8-2_DEGs.csv
# Output: results/rds/05_hdWGCNA_pre.rds, results/rds/05_scored.rds
###############################################################################

source("R/helpers.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(WGCNA)
  library(enrichR)
  library(patchwork)
  library(ggrepel)
})
future::plan("sequential")
set.seed(42)

# If EnrichR needs a proxy, configure it through environment variables
# (http_proxy / https_proxy) - never write credentials into scripts.

obj <- readRDS(rds_path("04_clustered.rds"))
celltypes_net <- c("NK", "Monocyte", "CD4", "CD8", "B cell")

# =============================================================================
# 1-2. Network construction on pre-ACT PBMCs
# =============================================================================
cl1 <- subset(obj, subset = timepoint == "Pre")
DefaultAssay(cl1) <- "RNA"

cl1 <- SetupForWGCNA(cl1, gene_select = "fraction", fraction = 0.05, wgcna_name = "PBMC_LTS")

cl1 <- MetacellsByGroups(
  seurat_obj  = cl1,
  group.by    = c("PBMC_identity", "patient code"),
  reduction   = "harmony",
  k           = 25,
  max_shared  = 10,
  ident.group = "PBMC_identity"
)
cl1 <- NormalizeMetacells(cl1)

cl1 <- SetDatExpr(cl1, group_name = celltypes_net, group.by = "PBMC_identity",
                  assay = "RNA", layer = "data")

cl1 <- TestSoftPowers(cl1, networkType = "signed")
write.csv(GetPowerTable(cl1), tab_path("hdWGCNA_soft_power_table.csv"), row.names = FALSE)
ggsave(fig_path("hdWGCNA_soft_powers.png"), wrap_plots(PlotSoftPowers(cl1), ncol = 2),
       width = 10, height = 8, dpi = 200)

cl1 <- ConstructNetwork(cl1, tom_name = "PBMC", tom_outdir = file.path(dirs$rds, "TOM"))

cl1 <- ModuleEigengenes(cl1, group.by.vars = "patient code")
cl1 <- ModuleConnectivity(cl1, group.by = "PBMC_identity", group_name = celltypes_net)
cl1 <- ResetModuleNames(cl1, new_name = "PBMC-M")
cl1 <- ResetModuleColors(cl1, new_colors = module_cols[paste0("PBMC-M", 1:8)],
                         wgcna_name = "PBMC_LTS")

modules <- GetModules(cl1)
mods <- setdiff(levels(modules$module), "grey")
write.csv(modules, tab_path("hdWGCNA_module_assignments_kME.csv"), row.names = FALSE)

# Harmonised module eigengenes -> metadata
MEs <- GetMEs(cl1, harmonized = TRUE)
cl1@meta.data <- cl1@meta.data[, !colnames(cl1@meta.data) %in% colnames(MEs)]
cl1 <- AddMetaData(cl1, MEs)

# =============================================================================
# 3. Differential module eigengenes, LTS vs STS (Supp Fig 1C labels)
# =============================================================================
DMEs <- FindDMEs(
  cl1,
  barcodes1  = colnames(cl1)[cl1$survival == "LTS"],
  barcodes2  = colnames(cl1)[cl1$survival == "STS"],
  test.use   = "wilcox",
  wgcna_name = "PBMC_LTS"
)
write.csv(DMEs, tab_path("hdWGCNA_DMEs_LTS_vs_STS.csv"), row.names = FALSE)

saveRDS(cl1, rds_path("05_hdWGCNA_pre.rds"))

# =============================================================================
# 4. Module scores for all cells (pre + post)
# -----------------------------------------------------------------------------
# Module scores (Fig 2C / Supp 1D correlations, and the "M2 signature" used in
# Fig 7C) use the genes of each module with kME >= min_kME.
# Set min_kME <- -Inf to score with all genes assigned to each module instead.
# =============================================================================
min_kME <- 0.5

module_genes <- lapply(mods, function(m) {
  g <- modules$gene_name[modules$module == m & modules[[paste0("kME_", m)]] >= min_kME]
  intersect(g, rownames(obj[["RNA"]]))
})
names(module_genes) <- mods
print(sapply(module_genes, length))
write.csv(stack(module_genes) %>% rename(gene = values, module = ind),
          tab_path("hdWGCNA_module_genes_kME_filtered.csv"), row.names = FALSE)

DefaultAssay(obj) <- "RNA"
set.seed(1)
obj <- AddModuleScore(obj, features = module_genes, name = "modscore_", assay = "RNA")
score_cols <- paste0("modscore_", gsub("-", "_", mods))          # e.g. modscore_PBMC_M2
colnames(obj@meta.data)[match(paste0("modscore_", seq_along(mods)), colnames(obj@meta.data))] <- score_cols

# Hub-gene signatures (top genes by kME) used for the subcluster dot plot (Fig 3B)
top_hubs <- function(m, n) {
  modules %>% filter(module == m) %>% arrange(desc(.data[[paste0("kME_", m)]])) %>%
    slice_head(n = n) %>% pull(gene_name)
}
hub_sets <- list(PBMC_M2 = top_hubs("PBMC-M2", 20), PBMC_M6 = top_hubs("PBMC-M6", 5))
# The published hub-gene sets were:
#   M2: EFHD2 GNLY CD8A ZEB2 FGFBP2 GZMA HLA-B AOAH PRF1 CTSW KLRD1 GZMH GZMB CCL5
#       METRNL GZMM TGFB1 NKG7 CST7 APOBEC3G
#   M6: LEF1 IL7R LTB CCR7 TSHZ2
print(hub_sets)
message("M2 hub genes identical to kME>=", min_kME, " set: ",
        setequal(hub_sets$PBMC_M2, module_genes[["PBMC-M2"]]))
message("M6 hub genes identical to kME>=", min_kME, " set: ",
        setequal(hub_sets$PBMC_M6, module_genes[["PBMC-M6"]]))

for (nm in names(hub_sets)) {
  obj <- AddModuleScore(obj, features = list(intersect(hub_sets[[nm]], rownames(obj))),
                        name = paste0("hub_", nm, "_score"))
  obj@meta.data[[paste0("hub_", nm)]] <- obj@meta.data[[paste0("hub_", nm, "_score1")]]
  obj@meta.data[[paste0("hub_", nm, "_score1")]] <- NULL
}

saveRDS(obj, rds_path("05_scored.rds"))

# =============================================================================
# 5a. Supp Fig 1C - module dendrogram
# =============================================================================
pdf(fig_path("SuppFig1C_module_dendrogram.pdf"), width = 10, height = 7)
PlotDendrogram(cl1, wgcna_name = "PBMC_LTS", main = "PBMC co-expression module dendrogram")
dev.off()

# =============================================================================
# 5b. Fig 2D - M2 and M6 eigengenes on the pre-ACT UMAP (shared colour scale)
# =============================================================================
scale_mods <- c("PBMC-M2", "PBMC-M5", "PBMC-M6")
all_vals <- unlist(cl1@meta.data[, scale_mods])
for (m in c("PBMC-M2", "PBMC-M6")) {
  p <- FeaturePlot(cl1, features = m, reduction = "umap.harmony", order = TRUE,
                   min.cutoff = min(all_vals, na.rm = TRUE),
                   max.cutoff = max(all_vals, na.rm = TRUE),
                   cols = c("grey95", module_cols[m]), raster = FALSE) +
    ggtitle(m) +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  save_plot(p, paste0("Fig2D_UMAP_", m, "_hME.png"), width = 5, height = 4, dpi = 300)
}

# =============================================================================
# 5c. Fig 2E - module gene UMAP (M2 and M6)
# =============================================================================
set.seed(123)
cl1 <- RunModuleUMAP(cl1, n_hubs = 10, n_neighbors = 15, min_dist = 0.1)
umap_df <- GetModuleUMAP(cl1)

genes_to_label <- c(hub_sets$PBMC_M2, hub_sets$PBMC_M6)
plot_df  <- umap_df %>% filter(module %in% c("PBMC-M2", "PBMC-M6"))
label_df <- plot_df %>% filter(gene %in% genes_to_label)

p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(colour = module, size = kME), alpha = 0.7) +
  geom_label_repel(data = label_df, aes(label = gene), colour = "black", fill = "white",
                   label.size = 0.5, size = 5, max.overlaps = Inf, box.padding = 1.2,
                   point.padding = 0.5, force = 5, force_pull = 0.3, max.iter = 50000,
                   max.time = 5, min.segment.length = 0, segment.colour = "black",
                   segment.alpha = 0.6, label.padding = unit(0.25, "lines"), seed = 42) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  scale_colour_manual(values = module_cols[c("PBMC-M2", "PBMC-M6")]) +
  scale_size(range = c(0.4, 3)) +
  coord_equal() +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
save_plot(p, "Fig2E_module_gene_UMAP_M2_M6.png", width = 20, height = 10, dpi = 200)

saveRDS(cl1, rds_path("05_hdWGCNA_pre.rds"))   # now includes the module UMAP

# =============================================================================
# 5d. Supp Fig 2 - GO enrichment (EnrichR) of M2 and M6
# =============================================================================
dbs <- c("GO_Biological_Process_2023", "GO_Biological_Process_2025", "KEGG_2026",
         "KEGG_2019_Human", "KEGG_2021_Human", "MSigDB_Hallmark_2020", "NCI_Nature_2016",
         "WikiPathways_2021_Human", "WikiPathways_2019_Human", "WikiPathways_2024_Human")
cl1 <- RunEnrichr(cl1, dbs = dbs, max_genes = 100)
enrich_df <- GetEnrichrTable(cl1)
write.csv(enrich_df, tab_path("SuppFig2_EnrichR_all_modules.csv"), row.names = FALSE)

# Top five significant GO BP terms per module, selected from the table above
go_selected <- bind_rows(
  data.frame(
    module = "M2",
    Term = c("Positive Regulation Of CD8-positive, Alpha-Beta T Cell Activation (GO:2001187)",
             "Antigen Processing And Presentation Of Peptide Antigen Via MHC Class Ib (GO:0002428)",
             "T-helper 1 Cell Differentiation (GO:0045063)",
             "Natural Killer Cell Mediated Immunity (GO:0002228)",
             "Regulation Of T Cell Mediated Cytotoxicity (GO:0001914)"),
    Adjusted.P.value = c(1.51E-05, 4.53E-05, 0.007191374, 1.94E-05, 1.29E-05)),
  data.frame(
    module = "M6",
    Term = c("Regulation Of Myeloid Leukocyte Differentiation (GO:0002761)",
             "Response To Interleukin-4 (GO:0070670)",
             "T-helper Cell Differentiation (GO:0042093)",
             "Regulation Of T Cell Differentiation In Thymus (GO:0033081)",
             "Regulation Of T Cell Proliferation (GO:0042129)"),
    Adjusted.P.value = c(0.011277417, 0.014581867, 0.002654448, 0.021284467, 0.004374614))
) %>%
  mutate(log_adjP = -log10(Adjusted.P.value),
         Term_wrapped = str_wrap(Term, width = 45),
         module = factor(module, levels = c("M2", "M6")))
go_selected$Term_wrapped <- factor(go_selected$Term_wrapped,
                                   levels = rev(unique(go_selected$Term_wrapped)))

p <- ggplot(go_selected, aes(x = module, y = Term_wrapped)) +
  geom_point(aes(size = log_adjP, fill = module), colour = "black", shape = 21,
             stroke = 1, alpha = 0.9) +
  scale_fill_manual(values = c(M2 = unname(module_cols["PBMC-M2"]),
                               M6 = unname(module_cols["PBMC-M6"]))) +
  scale_size(range = c(4, 14), name = expression(-log[10](P[FDR]))) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold"))
save_plot(p, "SuppFig2_GO_dotplot_M2_M6.png", width = 9, height = 9, dpi = 200)

# =============================================================================
# 5e. Fig 3B - M2 / M6 hub signatures across pre-ACT subclusters
# =============================================================================
dot_df <- obj@meta.data %>%
  filter(timepoint == "Pre", !is.na(subcluster_updated)) %>%
  select(subcluster_updated, hub_PBMC_M2, hub_PBMC_M6) %>%
  pivot_longer(-subcluster_updated, names_to = "module", values_to = "score") %>%
  mutate(module = recode(module, hub_PBMC_M2 = "PBMC-M2", hub_PBMC_M6 = "PBMC-M6")) %>%
  group_by(module, subcluster_updated) %>%
  summarise(mean_score = mean(score, na.rm = TRUE),
            pct_high   = mean(score > 0, na.rm = TRUE), .groups = "drop") %>%
  group_by(module) %>%
  mutate(mean_score_scaled = as.numeric(scale(mean_score))) %>%
  ungroup() %>%
  mutate(subcluster_updated = factor(subcluster_updated, levels = subcluster_order),
         module = factor(module, levels = c("PBMC-M2", "PBMC-M6")))
write.csv(dot_df, tab_path("Fig3B_module_scores_by_subcluster.csv"), row.names = FALSE)

p <- ggplot(dot_df, aes(x = subcluster_updated, y = module)) +
  geom_point(aes(size = pct_high, colour = mean_score_scaled), alpha = 0.9) +
  scale_colour_gradientn(colours = c("#011f4b", "#005b96", "#6497b1", "#b3cde0",
                                     "#ffefea", "#ffb09c", "#fe5757", "#900000"),
                         limits = c(-2, 2), oob = scales::squish) +
  scale_size(range = c(2, 8), labels = scales::percent_format()) +
  labs(x = "Subcluster", y = NULL, colour = "Scaled mean\nmodule score", size = "% cells > 0") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "Fig3B_module_dotplot_subclusters.png", width = 8, height = 4, dpi = 300)

# =============================================================================
# 5f. Fig 4B - Venn diagram: CD8-2 DEGs vs M2 hub genes (kME > 0.5)
# (reconstructed from the figure legend; the original code was not available)
# =============================================================================
cd8_2_degs <- read.csv(tab_path("Fig4B_CD8-2_DEGs.csv"))$gene
m2_hubs <- modules$gene_name[modules$module == "PBMC-M2" & modules$`kME_PBMC-M2` > 0.5]
overlap <- intersect(cd8_2_degs, m2_hubs)
write.csv(data.frame(gene = overlap), tab_path("Fig4B_overlap_CD8-2_DEGs_M2_hubs.csv"),
          row.names = FALSE)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
VennDiagram::venn.diagram(
  x = list(`CD8-2 DEGs` = cd8_2_degs, `PBMC-M2 hub genes` = m2_hubs),
  filename = fig_path("Fig4B_venn_CD8-2_M2.png"), imagetype = "png",
  fill = c("#5c59ff", module_cols[["PBMC-M2"]]), alpha = 0.5,
  cex = 2, cat.cex = 1.4, height = 2000, width = 2000, resolution = 300
)
