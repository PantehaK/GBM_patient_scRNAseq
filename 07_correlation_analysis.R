###############################################################################
# 07_correlation_analysis.R
#
# Patient-level summaries (pre- and post-ACT) and their association with
# overall survival (OS from diagnosis). Pearson correlation; for the scatter
# plots the simple linear-regression slope test gives the same P value.
#
#   Fig 2B, Supp 1A/1B   major cell-type proportions vs OS
#   Fig 2C, Supp 1D      hdWGCNA module scores vs OS
#   Fig 4A, Supp 3B/4    subcluster proportions vs OS
#   Fig 5B               D50 and clonal expansion (all T cells) vs OS
#   Fig 6B               proportion of ACT-derived clones vs OS
#   Fig 7C               M2 ratio (non-ACT / ACT-derived clones) vs OS
#   Supp 5A/5B           subcluster D50 vs OS
#
# NOTE: the plotting code in this script was written when the repository was
# assembled (the original scripts exported the tables only). Check the
# correlation coefficients / P values against the published panels.
#
# Input : results/rds/06_final.rds, tables written by 06_TCR_repertoire.R
# Output: results/tables/correlations_*.csv, results/figures/Fig*/Supp*
###############################################################################

source("R/helpers.R")
suppressPackageStartupMessages(library(Seurat))

obj  <- readRDS(rds_path("06_final.rds"))
meta <- obj@meta.data %>% mutate(timepoint = factor(timepoint, levels = c("Pre", "Post")))
os_col <- "OS_from_diagnosis_months"
tp_cols <- c(Pre = "#2986cc", Post = "#F28E2B")

# =============================================================================
# Helpers
# =============================================================================
# per_patient: one row per patient x timepoint x feature with columns
#   patient_code, timepoint, feature, value, OS_from_diagnosis_months
cor_with_os <- function(per_patient) {
  per_patient %>%
    filter(!is.na(value)) %>%
    group_by(feature, timepoint) %>%
    filter(n() >= 3) %>%
    summarise(
      n = n(),
      r = suppressWarnings(cor.test(value, .data[[os_col]], method = "pearson")$estimate),
      p = suppressWarnings(cor.test(value, .data[[os_col]], method = "pearson")$p.value),
      .groups = "drop"
    )
}

# Butterfly plot: bar length = -log10(P), left = pre-ACT, right = post-ACT,
# fill = Pearson r (red negative, blue positive), dashed lines at P = 0.05.
butterfly_plot <- function(cor_df, feature_order = NULL, title = NULL) {
  if (is.null(feature_order)) feature_order <- unique(cor_df$feature)
  df <- cor_df %>%
    mutate(feature = factor(feature, levels = rev(feature_order)),
           x = ifelse(timepoint == "Pre", -1, 1) * -log10(p),
           star = ifelse(p < 0.05, "*", ""))
  lim <- max(abs(df$x), -log10(0.05), na.rm = TRUE) * 1.15
  ggplot(df, aes(x = x, y = feature, fill = r)) +
    geom_col(colour = "black", linewidth = 0.3, width = 0.75) +
    geom_text(aes(label = star, hjust = ifelse(x < 0, 1.3, -0.3)), size = 6, vjust = 0.75) +
    geom_vline(xintercept = 0, colour = "black") +
    geom_vline(xintercept = c(-1, 1) * -log10(0.05), linetype = "dashed", colour = "red") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                         limits = c(-1, 1), name = "Pearson r") +
    scale_x_continuous(limits = c(-lim, lim), labels = function(b) abs(b)) +
    annotate("text", x = -lim * 0.6, y = length(feature_order) + 0.8, label = "Pre-ACT", fontface = "bold") +
    annotate("text", x =  lim * 0.6, y = length(feature_order) + 0.8, label = "Post-ACT", fontface = "bold") +
    coord_cartesian(clip = "off") +
    labs(x = expression(-log[10](P)), y = NULL, title = title) +
    theme_classic(base_size = 13)
}

# Scatter vs OS: pre = blue circles + dashed fit, post = orange triangles + solid
# fit, grey lines join the two samples of each patient.
prepost_scatter <- function(per_patient, ylab = "Value", ncol = 4) {
  ggplot(per_patient, aes(x = .data[[os_col]], y = value)) +
    geom_line(aes(group = patient_code), colour = "grey75") +
    geom_smooth(aes(colour = timepoint, linetype = timepoint), method = "lm",
                se = FALSE, formula = y ~ x) +
    geom_point(aes(colour = timepoint, shape = timepoint), size = 2.5) +
    scale_colour_manual(values = tp_cols) +
    scale_shape_manual(values = c(Pre = 16, Post = 17)) +
    scale_linetype_manual(values = c(Pre = "dashed", Post = "solid")) +
    facet_wrap(~feature, scales = "free_y", ncol = ncol) +
    labs(x = "OS from diagnosis (months)", y = ylab, colour = NULL, shape = NULL, linetype = NULL) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
}

# Single-feature regression, pre and post side by side (Fig 5B, 6B, 7C)
regression_plot <- function(per_patient, ylab) {
  stats <- cor_with_os(per_patient) %>%
    mutate(label = sprintf("r = %.2f, P = %.3g", r, p))
  ggplot(per_patient, aes(x = .data[[os_col]], y = value)) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey85") +
    geom_point(aes(colour = timepoint), size = 3) +
    geom_text(data = stats, aes(label = label), x = -Inf, y = Inf,
              hjust = -0.1, vjust = 1.5, inherit.aes = FALSE) +
    scale_colour_manual(values = tp_cols, guide = "none") +
    facet_grid(feature ~ timepoint, scales = "free_y") +
    labs(x = "OS from diagnosis (months)", y = ylab) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
}

run_block <- function(per_patient, name, ylab, feature_order = NULL, ncol = 4,
                      butterfly_file = NULL, scatter_file = NULL, w = 12, h = 8) {
  write.csv(per_patient, tab_path(paste0("patient_level_", name, ".csv")), row.names = FALSE)
  cors <- cor_with_os(per_patient)
  write.csv(cors, tab_path(paste0("correlations_", name, ".csv")), row.names = FALSE)
  if (!is.null(butterfly_file)) {
    save_plot(butterfly_plot(cors, feature_order), butterfly_file,
              width = 7, height = 1 + 0.35 * length(unique(cors$feature)))
  }
  if (!is.null(scatter_file)) {
    save_plot(prepost_scatter(per_patient, ylab, ncol), scatter_file, width = w, height = h)
  }
  invisible(cors)
}

stacked_patient_bar <- function(long, group_col, cols, file, w = 10, h = 5) {
  p <- ggplot(long, aes(x = patient_OS_label, y = 100 * proportion, fill = .data[[group_col]])) +
    geom_col(width = 0.85, colour = "black", linewidth = 0.2) +
    facet_wrap(~timepoint, ncol = 1) +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = "% of PBMCs", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, file, width = w, height = h)
}

# =============================================================================
# 1. Major cell types (Fig 2B, Supp 1A, Supp 1B)
# =============================================================================
keep_types <- c("NK", "CD4", "CD8", "Monocyte", "B cell", "pDC")
ct_long <- bind_rows(lapply(c("Pre", "Post"), function(tp) {
  make_prop_table(meta %>% filter(timepoint == tp), "PBMC_identity", keep_types) %>%
    mutate(timepoint = tp)
})) %>% mutate(timepoint = factor(timepoint, levels = c("Pre", "Post")))

stacked_patient_bar(ct_long, "PBMC_identity", pbmc_cols, "SuppFig1A_celltype_composition.png")
run_block(ct_long %>% transmute(patient_code, timepoint, feature = PBMC_identity,
                                value = proportion, OS_from_diagnosis_months),
          "celltype_proportion", "Proportion of PBMCs", feature_order = keep_types, ncol = 3,
          butterfly_file = "Fig2B_butterfly_celltypes.png",
          scatter_file = "SuppFig1B_scatter_celltypes.png", w = 10, h = 6)

# =============================================================================
# 2. hdWGCNA module scores (Fig 2C, Supp 1D)
# =============================================================================
score_cols <- grep("^modscore_PBMC_M[0-9]+$", colnames(meta), value = TRUE)
mod_pp <- meta %>%
  select(patient_code, timepoint, all_of(os_col), all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "feature", values_to = "score") %>%
  group_by(patient_code, timepoint, feature, OS_from_diagnosis_months) %>%
  summarise(value = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(feature = sub("^modscore_PBMC_", "PBMC-", feature))
run_block(mod_pp, "module_scores", "Mean module score",
          feature_order = paste0("PBMC-M", 1:8), ncol = 4,
          butterfly_file = "Fig2C_butterfly_modules.png",
          scatter_file = "SuppFig1D_scatter_modules.png", w = 12, h = 6)

# =============================================================================
# 3. Subclusters (Fig 4A, Supp 3B, Supp 4)
# =============================================================================
sub_long <- bind_rows(lapply(c("Pre", "Post"), function(tp) {
  make_prop_table(meta %>% filter(timepoint == tp), "subcluster_updated") %>%
    mutate(timepoint = tp)
})) %>% mutate(timepoint = factor(timepoint, levels = c("Pre", "Post")))

stacked_patient_bar(sub_long, "subcluster_updated", subcluster_cols,
                    "SuppFig3B_subcluster_composition.png", h = 7)
run_block(sub_long %>% transmute(patient_code, timepoint, feature = subcluster_updated,
                                 value = proportion, OS_from_diagnosis_months),
          "subcluster_proportion", "Proportion of PBMCs", feature_order = subcluster_order,
          ncol = 6, butterfly_file = "Fig4A_butterfly_subclusters.png",
          scatter_file = "SuppFig4_scatter_subclusters.png", w = 16, h = 10)

# =============================================================================
# 4. TCR diversity and expansion, all T cells (Fig 5B)
# =============================================================================
d50 <- read.csv(tab_path("Fig5B_D50_by_patient_timepoint.csv")) %>%
  filter(T_cell_group == "All_T_cells") %>%
  transmute(patient_code, timepoint = factor(timepoint, levels = c("Pre", "Post")),
            feature = "D50 (%)", value = d50_percent, OS_from_diagnosis_months)
expn <- read.csv(tab_path("Fig5B_expansion_by_patient_timepoint.csv")) %>%
  transmute(patient_code, timepoint = factor(timepoint, levels = c("Pre", "Post")),
            feature = "Expanded clonotypes (%)", value = percent_unique_clones_expanded,
            OS_from_diagnosis_months)
fig5b <- bind_rows(d50, expn)
write.csv(cor_with_os(fig5b), tab_path("correlations_Fig5B_D50_expansion.csv"), row.names = FALSE)
save_plot(regression_plot(fig5b, NULL), "Fig5B_regression_D50_expansion.png", width = 8, height = 7)

# =============================================================================
# 5. ACT-derived clones (Fig 6B) - CD8 version = table used for the paper
# =============================================================================
fig6b <- read.csv(tab_path("Fig6B_ACT_clone_proportion.csv")) %>%
  transmute(patient_code, timepoint = factor(timepoint, levels = c("Pre", "Post")),
            feature = paste0("ACT-derived clones (% of ", population, ")"),
            value = percent_product, OS_from_diagnosis_months)
write.csv(cor_with_os(fig6b), tab_path("correlations_Fig6B_ACT_clones.csv"), row.names = FALSE)
save_plot(regression_plot(fig6b, "% of cells"), "Fig6B_regression_ACT_clones.png", width = 8, height = 7)

# =============================================================================
# 6. M2 signature ratio, non-ACT / ACT-derived T-cell clones (Fig 7C)
# =============================================================================
fig7c <- read.csv(tab_path("Fig7C_M2_ACT_vs_nonACT_clones.csv")) %>%
  select(patient_code, timepoint, product_clone_label, mean_M2, all_of(os_col)) %>%
  pivot_wider(names_from = product_clone_label, values_from = mean_M2) %>%
  filter(!is.na(`CMV-product`), !is.na(`Non-CMV`)) %>%
  transmute(patient_code, timepoint = factor(timepoint, levels = c("Pre", "Post")),
            feature = "M2 ratio (non-ACT / ACT-derived)",
            value = `Non-CMV` / `CMV-product`, OS_from_diagnosis_months)
write.csv(fig7c, tab_path("patient_level_Fig7C_M2_ratio.csv"), row.names = FALSE)
write.csv(cor_with_os(fig7c), tab_path("correlations_Fig7C_M2_ratio.csv"), row.names = FALSE)
save_plot(regression_plot(fig7c, "M2 signature ratio"), "Fig7C_regression_M2_ratio.png",
          width = 8, height = 4)

# =============================================================================
# 7. Subcluster D50 (Supp 5A, 5B)
# =============================================================================
d50_sub <- read.csv(tab_path("SuppFig5_D50_by_subcluster.csv")) %>%
  filter(!is.na(subcluster_updated)) %>%
  transmute(patient_code, timepoint = factor(timepoint, levels = c("Pre", "Post")),
            feature = subcluster_updated, value = d50_percent, OS_from_diagnosis_months)
run_block(d50_sub, "subcluster_D50", "D50 (%)",
          feature_order = intersect(subcluster_order, unique(d50_sub$feature)), ncol = 6,
          butterfly_file = "SuppFig5A_butterfly_subcluster_D50.png",
          scatter_file = "SuppFig5B_scatter_subcluster_D50.png", w = 16, h = 10)

sessionInfo()
