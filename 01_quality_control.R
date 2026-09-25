###############################################################################
# 01_quality_control.R
#
# Per-sample preprocessing of Cell Ranger multi outputs (GEX + ADT):
#   1. Read each sample and compute QC metrics
#   2. Filter low-quality cells, remove TCR/BCR genes, SCTransform
#   3. Remove doublets (DoubletFinder)
#   4. Cell-cycle scoring
#   5. Merge all samples and add sample-level metadata
#
# Input : data/raw/GBM*/**/sample_filtered_feature_bc_matrix
#         data/metadata/Sample_information.xlsx
# Output: results/rds/01_merged_qc.rds
#         results/tables/qc_summary_{raw,post_filter,post_doublet}.csv
###############################################################################

source("R/helpers.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(readxl)
})
set.seed(42)

# QC thresholds
qc_thresholds <- list(
  percent.mt       = 15,
  percent.ribo     = 60,
  percent.hb       = 0.1,
  nFeature_RNA_min = 200,
  nCount_RNA_min   = 200
)
doublet_rate <- 0.08   # expected doublet rate used by DoubletFinder

qc_summary <- function(objs) {
  bind_rows(lapply(names(objs), function(nm) {
    o <- objs[[nm]]
    tibble(sample = nm, nCells = ncol(o),
           median_nCount_RNA = median(o$nCount_RNA),
           median_nFeature_RNA = median(o$nFeature_RNA),
           median_log10_UMI = median(log10(o$nCount_RNA + 1)))
  }))
}

# =============================================================================
# 1. Read data
# =============================================================================
batch_folders <- list.dirs(dirs$raw, recursive = FALSE, full.names = TRUE) %>%
  grep("^.*/GBM", ., value = TRUE)
sample_paths <- unlist(lapply(batch_folders, function(b) {
  list.dirs(b, recursive = TRUE, full.names = TRUE) %>%
    grep("sample_filtered_feature_bc_matrix$", ., value = TRUE)
}))

read_sample <- function(data_dir, sample_name) {
  data <- Read10X(data.dir = data_dir)
  obj <- CreateSeuratObject(counts = data$`Gene Expression`, project = sample_name)
  obj[["RNA"]] <- CreateAssay5Object(counts = data$`Gene Expression`, min.cells = 3)
  if ("Antibody Capture" %in% names(data)) {
    obj[["ADT"]] <- CreateAssay5Object(counts = data$`Antibody Capture`)
  }
  obj$sample <- sample_name
  obj$sample_id <- sample_name
  obj
}

seurat_objects <- list()
for (path in sample_paths) {
  sample_name <- basename(dirname(dirname(path)))
  message("Reading: ", sample_name)
  seurat_objects[[sample_name]] <- read_sample(path, sample_name)
}
write.csv(qc_summary(seurat_objects), tab_path("qc_summary_raw.csv"), row.names = FALSE)

# =============================================================================
# 2. QC filtering, removal of immune receptor genes, SCTransform
# =============================================================================
qc_filter <- function(obj) {
  obj <- PercentageFeatureSet(obj, pattern = "^MT-",        col.name = "percent.mt")
  obj <- PercentageFeatureSet(obj, pattern = "^RPS|^RPL",   col.name = "percent.ribo")
  obj <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]",   col.name = "percent.hb")
  subset(obj, subset = percent.ribo < qc_thresholds$percent.ribo &
                       percent.hb   < qc_thresholds$percent.hb &
                       percent.mt   < qc_thresholds$percent.mt &
                       nFeature_RNA > qc_thresholds$nFeature_RNA_min &
                       nCount_RNA   > qc_thresholds$nCount_RNA_min)
}

# TCR/BCR variable genes are removed so that clonotype identity does not drive
# clustering. Constant IGH genes (IGHM/G/D/E) are removed as well.
remove_receptor_genes <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  receptor_genes <- c(grep("^TR[AB]|^IGL|^IGK|^IGHV", rownames(obj), value = TRUE),
                      grep("^IGH[MGDHE]$", rownames(obj), value = TRUE))
  keep <- setdiff(rownames(obj[["RNA"]]), receptor_genes)
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[keep, ]
  obj[["RNA"]] <- CreateAssay5Object(counts = counts)
  obj
}

for (nm in names(seurat_objects)) {
  message("QC: ", nm)
  obj <- seurat_objects[[nm]]
  obj <- qc_filter(obj)
  obj <- remove_receptor_genes(obj)
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", variable.features.n = 2000,
                     ncells = 3000, return.only.var.genes = FALSE, verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)

  p <- VlnPlot(obj, features = c("percent.mt", "percent.ribo", "percent.hb",
                                 "nCount_RNA", "nFeature_RNA"), group.by = "sample", ncol = 3)
  ggsave(file.path(dirs$qc, paste0(nm, "_qc_violin.png")), p, width = 12, height = 8, dpi = 150)

  seurat_objects[[nm]] <- obj
}
write.csv(qc_summary(seurat_objects), tab_path("qc_summary_post_filter.csv"), row.names = FALSE)

# =============================================================================
# 3. Doublet removal (DoubletFinder, homotypic-adjusted)
# =============================================================================
for (nm in names(seurat_objects)) {
  obj <- seurat_objects[[nm]]
  if (ncol(obj) < 100) { message("Skipping DoubletFinder for ", nm, " (<100 cells)"); next }
  message("DoubletFinder: ", nm)

  DefaultAssay(obj) <- "SCT"
  obj <- FindNeighbors(obj, dims = 1:15, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)

  sweep <- paramSweep(obj, PCs = 1:10, sct = TRUE)
  bcmvn <- find.pK(summarizeSweep(sweep, GT = FALSE))
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  homotypic_prop <- modelHomotypic(obj$seurat_clusters)
  nExp <- round(round(doublet_rate * ncol(obj)) * (1 - homotypic_prop))

  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp, sct = TRUE)
  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[df_col]] == "Singlet"])
  message("  retained ", ncol(obj), " singlets (pK = ", pK, ")")
  seurat_objects[[nm]] <- obj
}
write.csv(qc_summary(seurat_objects), tab_path("qc_summary_post_doublet.csv"), row.names = FALSE)

# =============================================================================
# 4. Cell-cycle scoring
# =============================================================================
for (nm in names(seurat_objects)) {
  obj <- seurat_objects[[nm]]
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(obj, s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  seurat_objects[[nm]] <- obj
}

# =============================================================================
# 5. Merge samples and add sample-level metadata
# =============================================================================
for (nm in names(seurat_objects)) {
  obj <- RenameCells(seurat_objects[[nm]], add.cell.id = nm)   # barcodes become <sample>_<barcode>
  obj$id <- nm
  seurat_objects[[nm]] <- obj
}
merged <- merge(x = seurat_objects[[1]], y = seurat_objects[-1])

# Sample sheet: one row per sample (column `orig.ident` = sample name).
# Only the columns needed downstream are added here.
sample_metadata <- as.data.frame(read_excel(file.path(dirs$meta, "Sample_information.xlsx")))
rownames(sample_metadata) <- sample_metadata$orig.ident
for (col in c("patient code", "timepoint", "survival")) {
  merged[[col]] <- sample_metadata[merged$id, col]
}

merged@meta.data <- merged@meta.data %>%
  select(-starts_with("pANN_"), -starts_with("DF.classifications_"))

saveRDS(merged, rds_path("01_merged_qc.rds"))
message("Done: ", ncol(merged), " cells from ", length(unique(merged$id)), " samples")
