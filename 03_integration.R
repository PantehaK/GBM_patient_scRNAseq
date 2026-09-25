###############################################################################
# 03_integration.R
#
# SCTransform (per sequencing run) + Harmony integration on sequencing run.
# Integrating on sequencing run rather than GEM batch avoided over-correction.
# Contaminating / damaged clusters are removed and the data re-integrated:
#   Round 1 (res 0.4): remove clusters 8, 9, 12 (contaminants; 12 = heat-shocked)
#   Round 2 (res 0.5): remove cluster 12
#   Round 3          : final integration used for all downstream analyses
#
# Cluster numbers refer to the run used for the paper; if you re-run from
# scratch, check the marker tables written at each round before removing.
#
# Input : results/rds/02_merged_qc_tcr.rds
# Output: results/rds/03_integrated.rds
#         results/tables/markers_round{1,2}.csv
###############################################################################

source("R/helpers.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})
set.seed(42)

obj <- readRDS(rds_path("02_merged_qc_tcr.rds"))

obj$sequencing_run <- factor(case_when(
  obj$batch %in% c("GBM3", "GBM6")                 ~ "run_3",
  obj$batch %in% c("GBM2", "GBM5", "GBM9")         ~ "run_2",
  obj$batch %in% c("GBM1", "GBM4", "GBM7", "GBM8") ~ "run_1",
  TRUE ~ NA_character_
))
stopifnot(!any(is.na(obj$sequencing_run)))

integrate <- function(obj, return_only_var_genes = TRUE) {
  DefaultAssay(obj) <- "RNA"
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sequencing_run)
  obj <- SCTransform(obj, assay = "RNA", new.assay.name = "SCT",
                     vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
                     return.only.var.genes = return_only_var_genes, verbose = FALSE)
  obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- RunHarmony(obj, group.by.vars = "sequencing_run", assay.use = "SCT", verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, graph.name = "harmony_snn")
  obj
}

cluster_and_markers <- function(obj, resolution, out_csv) {
  obj <- FindClusters(obj, resolution = resolution, graph.name = "harmony_snn",
                      cluster.name = "harmony_clusters")
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
  markers <- FindAllMarkers(PrepSCTFindMarkers(obj, assay = "SCT"),
                            only.pos = TRUE, assay = "SCT", slot = "data") %>%
    filter(avg_log2FC > 0.5, !grepl("^RPS|^RPL|^MT-", gene))
  write.csv(markers, tab_path(out_csv), row.names = FALSE)
  obj
}

# Round 1
obj <- integrate(obj, return_only_var_genes = TRUE)
obj <- cluster_and_markers(obj, 0.4, "markers_round1.csv")
obj <- subset(obj, subset = !(harmony_clusters %in% c(8, 9, 12)))

# Round 2
obj <- integrate(obj, return_only_var_genes = FALSE)
obj <- cluster_and_markers(obj, 0.5, "markers_round2.csv")
obj <- subset(obj, subset = !(harmony_clusters %in% 12))

# Round 3 (final)
obj <- integrate(obj, return_only_var_genes = TRUE)

DefaultAssay(obj) <- "ADT"
obj[["ADT"]] <- JoinLayers(obj[["ADT"]])
obj <- NormalizeData(obj, assay = "ADT", normalization.method = "CLR")
DefaultAssay(obj) <- "SCT"

saveRDS(obj, rds_path("03_integrated.rds"))
