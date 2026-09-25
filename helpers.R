###############################################################################
# helpers.R
# Shared paths, colour palettes and helper functions used by every script.
# Source this from the repository root:  source("R/helpers.R")
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
# Set GBM_PROJECT_DIR to the folder that contains data/ and results/,
# or run the scripts from the repository root.
PROJECT_DIR <- Sys.getenv("GBM_PROJECT_DIR", unset = normalizePath("."))

dirs <- list(
  raw     = file.path(PROJECT_DIR, "data", "raw"),        # Cell Ranger multi outputs (GBM*/...)
  meta    = file.path(PROJECT_DIR, "data", "metadata"),   # sample sheet, survival table, product TCR keys
  rds     = file.path(PROJECT_DIR, "results", "rds"),
  tables  = file.path(PROJECT_DIR, "results", "tables"),
  figures = file.path(PROJECT_DIR, "results", "figures"),
  qc      = file.path(PROJECT_DIR, "results", "figures", "qc")
)
invisible(lapply(dirs[-(1:2)], dir.create, recursive = TRUE, showWarnings = FALSE))

rds_path <- function(x) file.path(dirs$rds, x)
tab_path <- function(x) file.path(dirs$tables, x)
fig_path <- function(x) file.path(dirs$figures, x)

options(future.globals.maxSize = 8 * 1024^3)

# -----------------------------------------------------------------------------
# 2. Colour palettes
# -----------------------------------------------------------------------------
pbmc_cols <- c(
  "CD8"      = "#3531ff",
  "CD4"      = "#ff3434",
  "NK"       = "#904ad2",
  "Monocyte" = "#F5A56B",
  "B cell"   = "#31c7ff",
  "pDC"      = "#00C2A8"
)

module_cols <- c(
  "PBMC-M1" = "#6CC3F4", "PBMC-M2" = "#B560DD", "PBMC-M3" = "#f467ba",
  "PBMC-M4" = "#9BE599", "PBMC-M5" = "#DD9560", "PBMC-M6" = "#E65757",
  "PBMC-M7" = "#72BEB7", "PBMC-M8" = "#49d10b", "grey" = "grey80"
)

# Subcluster map: harmony_clusters2 (resolution 0.8) -> published subcluster
# name and colour. Numbering restarts within each lineage; old cluster 17 = pDC.
subcluster_map <- data.frame(
  old    = c("0","4","3","5","12","11","7",   # CD4
             "9","1","2",                     # CD8
             "6","14",                        # NK
             "8","15","13",                   # Monocytes
             "10","16",                       # B cells
             "17"),                           # pDC
  subset = c(rep("CD4", 7), rep("CD8", 3), rep("NK", 2), rep("Mono", 3), rep("B", 2), "pDC"),
  colour = c("#ff3434","#e31a1c","#fb6a4a","#fc9272","#fcbba1","#f768a1","#ff9ec4",
             "#3531ff","#5c59ff","#8a87ff",
             "#904ad2","#b07be0",
             "#F28E2B","#E07B1F","#C76A14",
             "#6fd8ff","#66b3cf",
             "#00C2A8"),
  stringsAsFactors = FALSE
)
subcluster_map$idx <- ave(seq_along(subcluster_map$subset), subcluster_map$subset, FUN = seq_along)
subcluster_map$new <- ifelse(subcluster_map$subset == "pDC", "pDC",
                             paste0(subcluster_map$subset, "-", subcluster_map$idx))
subcluster_cols  <- setNames(subcluster_map$colour, subcluster_map$new)
subcluster_order <- subcluster_map$new

# -----------------------------------------------------------------------------
# 3. Survival metadata
# -----------------------------------------------------------------------------
surv_cols <- c("OS_from_diagnosis_months", "OS_from_consent_months",
               "PFS_from_diagnosis_months", "PFS_from_consent_months")

# data/metadata/patient_survival.csv must contain `patient_num` (the digits after
# the dash in the `patient code` column, e.g. "2105-01" -> "01"), surv_cols and
# `survival_group` (LTS / STS; only used for the hdWGCNA DME test).
read_survival <- function(file = file.path(dirs$meta, "patient_survival.csv")) {
  df <- read.csv(file, colClasses = c(patient_num = "character"), stringsAsFactors = FALSE)
  stopifnot(all(c("patient_num", "survival_group", surv_cols) %in% colnames(df)))
  df
}

# Adds OS/PFS columns, `patient_code` and an OS-ordered factor `patient_OS_label`
# ("<patient code> (<OS> mo)") to a Seurat object.
add_survival_metadata <- function(obj, survival_df = read_survival()) {
  md <- obj@meta.data
  md$patient_code <- trimws(as.character(md$`patient code`))
  idx <- match(sub("^.*-", "", md$patient_code), survival_df$patient_num)
  if (any(is.na(idx))) stop("Some patients were not found in patient_survival.csv")
  for (col in surv_cols) md[[col]] <- survival_df[[col]][idx]
  md$survival <- survival_df$survival_group[idx]   # final LTS/STS grouping

  pt_tbl <- md %>%
    distinct(patient_code, OS_from_diagnosis_months) %>%
    arrange(OS_from_diagnosis_months) %>%
    mutate(label = paste0(patient_code, " (", OS_from_diagnosis_months, " mo)"))
  md$patient_OS_label <- factor(pt_tbl$label[match(md$patient_code, pt_tbl$patient_code)],
                                levels = pt_tbl$label)
  obj@meta.data <- md
  obj
}

# Rainbow (Spectral) palette used for "coloured by OS" panels.
patient_os_cols <- function(labels) {
  labels <- levels(droplevels(factor(labels)))
  setNames(hcl.colors(length(labels), "Spectral", rev = TRUE), labels)
}

# One row per patient with survival columns.
patient_info_table <- function(meta) {
  meta %>%
    distinct(patient_OS_label, patient_code, across(all_of(surv_cols))) %>%
    arrange(OS_from_diagnosis_months)
}

# -----------------------------------------------------------------------------
# 4. Composition tables (proportion of `group_col` per patient)
# -----------------------------------------------------------------------------
make_prop_table <- function(meta, group_col, levels_keep = NULL) {
  df <- meta %>% filter(!is.na(.data[[group_col]]))
  if (!is.null(levels_keep)) {
    missing_lv <- setdiff(levels_keep, unique(as.character(df[[group_col]])))
    if (length(missing_lv) > 0) stop("Not found in ", group_col, ": ", paste(missing_lv, collapse = ", "))
    df <- df %>% filter(.data[[group_col]] %in% levels_keep)
    df[[group_col]] <- factor(df[[group_col]], levels = levels_keep)
  }
  out <- df %>%
    count(patient_OS_label, group = .data[[group_col]], name = "n_cells") %>%
    complete(patient_OS_label, group, fill = list(n_cells = 0)) %>%
    group_by(patient_OS_label) %>%
    mutate(total_cells = sum(n_cells), proportion = n_cells / total_cells) %>%
    ungroup() %>%
    as.data.frame()
  names(out)[names(out) == "group"] <- group_col
  out %>%
    left_join(patient_info_table(meta), by = "patient_OS_label") %>%
    select(patient_OS_label, patient_code, all_of(group_col),
           n_cells, total_cells, proportion, all_of(surv_cols)) %>%
    arrange(OS_from_diagnosis_months, .data[[group_col]])
}

# -----------------------------------------------------------------------------
# 5. TCR helpers
# -----------------------------------------------------------------------------
# Clonotype used for diversity / expansion / sharing: TRB V + J + CDR3 (aa).
make_clone_id <- function(v, j, cdr3) {
  ok <- !is.na(v) & v != "" & !is.na(j) & j != "" & !is.na(cdr3) & cdr3 != ""
  ifelse(ok, paste(v, j, cdr3, sep = "_"), NA_character_)
}

# Normalised TRB key used to match single-cell TCRs to the ACT product TCRs
# (format "CDR3|TRBVx-y|TRBJx-y", allele suffixes and leading zeros removed).
norm_trbv <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[x == ""] <- NA_character_
  x <- x %>% str_replace("^TCRBV", "TRBV") %>% str_remove("\\*\\d+$")
  m <- str_match(x, "TRBV(\\d+)-(\\d+)")
  ifelse(is.na(m[, 1]), x, paste0("TRBV", as.integer(m[, 2]), "-", as.integer(m[, 3])))
}
norm_trbj <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[x == ""] <- NA_character_
  x <- x %>% str_replace("^TCRBJ", "TRBJ") %>% str_remove("\\*\\d+$")
  m <- str_match(x, "TRBJ(\\d+)-(\\d+)")
  ifelse(is.na(m[, 1]), x, paste0("TRBJ", as.integer(m[, 2]), "-", as.integer(m[, 3])))
}
norm_cdr3 <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x[x == ""] <- NA_character_
  x
}
make_trb_key <- function(cdr3, v_gene, j_gene) {
  str_c(norm_cdr3(cdr3), norm_trbv(v_gene), norm_trbj(j_gene), sep = "|")
}

# D50: percentage of unique clonotypes that together account for 50% of cells.
calc_d50 <- function(clones) {
  clones <- clones[!is.na(clones) & clones != ""]
  if (length(clones) == 0) {
    return(tibble(n_tcr_cells = 0L, n_clonotypes = 0L,
                  d50_clone_count = NA_integer_, d50_percent = NA_real_))
  }
  clone_counts <- sort(table(clones), decreasing = TRUE)
  cumulative_fraction <- cumsum(clone_counts) / sum(clone_counts)
  d50_clone_count <- which(cumulative_fraction >= 0.5)[1]
  tibble(
    n_tcr_cells     = sum(clone_counts),
    n_clonotypes    = length(clone_counts),
    d50_clone_count = d50_clone_count,
    d50_percent     = 100 * d50_clone_count / length(clone_counts)
  )
}

# -----------------------------------------------------------------------------
# 6. Misc
# -----------------------------------------------------------------------------
save_plot <- function(p, file, width, height, dpi = 400, ...) {
  ggsave(filename = fig_path(file), plot = p, width = width, height = height,
         dpi = dpi, bg = "white", ...)
}

# Average expression table (genes x groups) in a fixed gene order.
avg_expression_table <- function(obj, genes, group_by, assay = "RNA", groups = NULL) {
  if (!is.null(groups)) {
    obj <- obj[, obj@meta.data[[group_by]] %in% groups]
  }
  genes_present <- intersect(genes, rownames(obj[[assay]]))
  missing <- setdiff(genes, genes_present)
  if (length(missing) > 0) message("Genes not found: ", paste(missing, collapse = ", "))

  grp <- droplevels(factor(obj@meta.data[[group_by]]))
  obj@meta.data[[group_by]] <- grp
  avg <- Seurat::AverageExpression(obj, features = genes_present, group.by = group_by,
                                   assays = assay, layer = "data")[[1]]
  avg <- as.matrix(avg)
  # Seurat v5 replaces "_" with "-" and prefixes numeric names with "g"
  seurat_names <- gsub("_", "-", levels(grp))
  seurat_names <- ifelse(grepl("^[0-9]", seurat_names), paste0("g", seurat_names), seurat_names)
  avg <- avg[genes_present, seurat_names, drop = FALSE]
  colnames(avg) <- levels(grp)
  if (!is.null(groups)) avg <- avg[, intersect(groups, colnames(avg)), drop = FALSE]
  avg
}

# Row z-scored heatmap of an average-expression matrix.
plot_avg_heatmap <- function(avg, file, width = 6, height = 8) {
  z <- t(scale(t(avg)))
  z[is.na(z)] <- 0
  cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  pheatmap::pheatmap(z, color = cols, breaks = seq(-2, 2, length.out = 102),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     border_color = NA, filename = fig_path(file),
                     width = width, height = height)
}
