###############################################################################
# 06_TCR_repertoire.R
#
# Clonotype = TRB V gene + J gene + CDR3 (aa).
#   1. Clonotype IDs and ACT-product (CMV-specific) clone flags
#   2. Diversity (D50) per patient and timepoint          -> Fig 5B, Supp 5
#   3. Clonal expansion per patient and timepoint          -> Fig 5B
#   4. Clonal expansion UMAP                               -> Fig 5A
#   5. Pre/post clonotype sharing (circos)                 -> Fig 5C
#   6. ACT-derived clones: UMAP, abundance, persistence    -> Fig 6A-C
#   7. ACT-derived clones across subclusters, M2 in ACT vs
#      non-ACT clones                                      -> Fig 7A, 7C, Supp 6
#
# Input : results/rds/05_scored.rds
#         data/metadata/product_clones_pre.csv, product_clones_post.csv
#           (TRB keys "CDR3|TRBVx-y|TRBJx-y" of single-cell clonotypes that
#            matched the TCR repertoire of the patient's ACT product; see README)
# Output: results/rds/06_final.rds, tables/figures listed below
###############################################################################

source("R/helpers.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(circlize)
})
set.seed(1)

obj <- readRDS(rds_path("05_scored.rds"))

# =============================================================================
# 1. Clonotypes and ACT-product clones
# =============================================================================
obj$clone_id <- make_clone_id(obj$TRB_v_gene, obj$TRB_j_gene, obj$TRB_cdr3)
obj$TRB_key  <- make_trb_key(obj$TRB_cdr3, obj$TRB_v_gene, obj$TRB_j_gene)

# Export of single-cell TRB clonotypes, used as input for matching against the
# ACT product TCR repertoire (matching was done outside this repository).
obj@meta.data %>%
  filter(!is.na(TRB_key)) %>%
  distinct(patient_code, timepoint, TRB_key, TRB_cdr3, TRB_v_gene, TRB_j_gene) %>%
  write.csv(tab_path("single_cell_TRB_clonotypes.csv"), row.names = FALSE)

read_keys <- function(f) {
  k <- unique(read.csv(file.path(dirs$meta, f))$key)
  k[!is.na(k) & k != ""]
}
pre_product_keys  <- read_keys("product_clones_pre.csv")
post_product_keys <- read_keys("product_clones_post.csv")

obj$pre_ACT  <- obj$TRB_key %in% pre_product_keys
obj$post_ACT <- obj$TRB_key %in% post_product_keys
# timepoint-matched flag: pre cells vs pre list, post cells vs post list
obj$product_clone <- case_when(obj$timepoint == "Pre"  ~ obj$pre_ACT,
                               obj$timepoint == "Post" ~ obj$post_ACT,
                               TRUE ~ FALSE)
obj$product_clone_label <- ifelse(obj$product_clone, "CMV-product", "Non-CMV")
print(table(obj$timepoint, obj$product_clone_label))

meta <- obj@meta.data %>% rownames_to_column("cell")
tcells <- meta %>% filter(PBMC_identity %in% c("CD4", "CD8"))
pt_info <- patient_info_table(meta)

# =============================================================================
# 2. D50 diversity
# =============================================================================
# Fig 5B (left): per patient and timepoint, all T cells (and CD4 / CD8 separately)
d50_input <- tcells %>% filter(!is.na(clone_id))
d50_T <- bind_rows(d50_input %>% mutate(T_cell_group = PBMC_identity),
                   d50_input %>% mutate(T_cell_group = "All_T_cells")) %>%
  group_by(patient_OS_label, patient_code, timepoint, T_cell_group) %>%
  summarise(calc_d50(clone_id), .groups = "drop") %>%
  left_join(pt_info, by = c("patient_OS_label", "patient_code")) %>%
  arrange(timepoint, T_cell_group, OS_from_diagnosis_months)
write.csv(d50_T, tab_path("Fig5B_D50_by_patient_timepoint.csv"), row.names = FALSE)

# Supp Fig 5: per subcluster (all cells of the subcluster with a TRB)
d50_sub <- meta %>%
  filter(!is.na(clone_id)) %>%
  group_by(patient_OS_label, patient_code, timepoint, subcluster_updated) %>%
  summarise(calc_d50(clone_id), .groups = "drop") %>%
  left_join(pt_info, by = c("patient_OS_label", "patient_code"))
write.csv(d50_sub, tab_path("SuppFig5_D50_by_subcluster.csv"), row.names = FALSE)

# =============================================================================
# 3. Clonal expansion (Fig 5B right)
# -----------------------------------------------------------------------------
# NOTE: expanded = clone detected in MORE than 3 cells (>= 4), as in the code
# used for the paper. The figure legend says ">= 3 cells" - make these agree.
# =============================================================================
expansion_min_cells <- 4

expansion <- tcells %>%
  filter(!is.na(TRB_key)) %>%
  count(patient_OS_label, patient_code, timepoint, TRB_key, name = "clone_cell_count") %>%
  mutate(expanded_clone = clone_cell_count >= expansion_min_cells) %>%
  group_by(patient_OS_label, patient_code, timepoint) %>%
  summarise(
    n_unique_TRB_clones = n(),
    n_expanded_TRB_clones = sum(expanded_clone),
    percent_unique_clones_expanded = 100 * n_expanded_TRB_clones / n_unique_TRB_clones,
    n_T_cells_with_TRB = sum(clone_cell_count),
    n_T_cells_in_expanded_clones = sum(clone_cell_count[expanded_clone]),
    percent_T_cells_in_expanded_clones = 100 * n_T_cells_in_expanded_clones / n_T_cells_with_TRB,
    .groups = "drop"
  ) %>%
  left_join(pt_info, by = c("patient_OS_label", "patient_code"))
write.csv(expansion, tab_path("Fig5B_expansion_by_patient_timepoint.csv"), row.names = FALSE)

# =============================================================================
# 4. Fig 5A - clonal expansion UMAP (pre-ACT T cells)
# -----------------------------------------------------------------------------
# Clone size = % of the patient's TRB+ CD4 (or CD8) T cells belonging to the clone
# =============================================================================
umap <- as.data.frame(Embeddings(obj, "umap.harmony"))
colnames(umap) <- c("UMAP_1", "UMAP_2")
umap$cell <- rownames(umap)

pre_T <- tcells %>% filter(timepoint == "Pre") %>% left_join(umap, by = "cell")
clone_size <- pre_T %>%
  filter(!is.na(clone_id)) %>%
  count(patient_code, PBMC_identity, clone_id, name = "clone_size") %>%
  group_by(patient_code, PBMC_identity) %>%
  mutate(clone_percent = 100 * clone_size / sum(clone_size)) %>%
  ungroup()

fg <- pre_T %>%
  inner_join(clone_size, by = c("patient_code", "PBMC_identity", "clone_id")) %>%
  mutate(UMAP_1 = UMAP_1 + rnorm(n(), 0, 0.02),   # tiny jitter so stacked clones are visible
         UMAP_2 = UMAP_2 + rnorm(n(), 0, 0.02)) %>%
  arrange(clone_percent)                          # largest clones drawn on top
os_cols <- patient_os_cols(meta$patient_OS_label)

p <- ggplot() +
  geom_point(data = pre_T, aes(UMAP_1, UMAP_2), colour = "grey90", size = 0.15) +
  geom_point(data = fg, aes(UMAP_1, UMAP_2, fill = patient_OS_label, size = clone_percent),
             shape = 21, colour = "black", stroke = 0.1, alpha = 0.7) +
  scale_fill_manual(values = os_cols, name = "Patient\n(OS from diagnosis)", drop = FALSE) +
  scale_size_area(name = "Clone size\n(% of patient's\nCD4/CD8 cells)", max_size = 1.8,
                  breaks = c(1, 5, 10, 20)) +
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1)),
         size = guide_legend(override.aes = list(fill = "grey40"))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 10) +
  theme(legend.key.size = unit(0.4, "cm"))
save_plot(p, "Fig5A_UMAP_clonal_expansion_pre.png", width = 7, height = 6)

# =============================================================================
# 5. Fig 5C - pre/post clonotype overlap (circos, one per patient)
# =============================================================================
circos_dir <- fig_path("Fig5C_circos")
dir.create(circos_dir, showWarnings = FALSE)

clones_by_tp <- meta %>%
  filter(!is.na(clone_id)) %>%
  distinct(patient_code, OS_from_diagnosis_months, timepoint, clone_id)

paired <- clones_by_tp %>%
  group_by(patient_code, OS_from_diagnosis_months) %>%
  filter(all(c("Pre", "Post") %in% timepoint)) %>%
  distinct(patient_code, OS_from_diagnosis_months) %>%
  arrange(OS_from_diagnosis_months)

overlap_tables <- list()
for (i in seq_len(nrow(paired))) {
  pid <- paired$patient_code[i]
  pre_ids  <- clones_by_tp %>% filter(patient_code == pid, timepoint == "Pre")  %>% pull(clone_id)
  post_ids <- clones_by_tp %>% filter(patient_code == pid, timepoint == "Post") %>% pull(clone_id)
  n_shared <- length(intersect(pre_ids, post_ids))

  tbl <- tibble(
    patient_code = pid, OS_from_diagnosis_months = paired$OS_from_diagnosis_months[i],
    total_pre = length(pre_ids), total_post = length(post_ids), n_shared = n_shared,
    pre_unique_prop  = 1 - n_shared / length(pre_ids),
    post_unique_prop = 1 - n_shared / length(post_ids),
    post_shared_prop = n_shared / length(post_ids),
    post_shared_percent = 100 * n_shared / length(post_ids)
  )
  overlap_tables[[pid]] <- tbl

  # shared ribbon = proportion of POST clonotypes also found pre
  edges <- data.frame(
    from  = c("Pre ACT", "Pre ACT", "Post ACT"),
    to    = c("Post ACT", "Pre ACT", "Post ACT"),
    value = c(tbl$post_shared_prop, tbl$pre_unique_prop, tbl$post_unique_prop),
    col   = c("#0b5394", "#d9d9d9", "#f4a261"),
    z     = c(10, 1, 2)
  ) %>% filter(value > 0) %>% arrange(z)

  png(file.path(circos_dir, sprintf("%02d_%s.png", i, gsub("[^A-Za-z0-9_-]", "_", pid))),
      width = 2400, height = 2400, res = 300)
  circos.clear()
  circos.par(gap.after = c(2, 10), track.margin = c(0.001, 0.001), cell.padding = c(0, 0, 0, 0))
  chordDiagram(edges[, c("from", "to", "value")], order = c("Pre ACT", "Post ACT"),
               grid.col = c("Pre ACT" = "#bcbcbc", "Post ACT" = "#ce7e00"),
               col = edges$col, link.zindex = edges$z, transparency = 0.15,
               self.link = 2, reduce = -1, annotationTrack = "grid",
               annotationTrackHeight = mm_h(c(12, 4)))
  title(main = paste0(pid, " (OS ", tbl$OS_from_diagnosis_months, " mo) | post clonotypes found pre: ",
                      round(tbl$post_shared_percent, 1), "%"), cex.main = 1.2)
  dev.off()
  circos.clear()
}
write.csv(bind_rows(overlap_tables), tab_path("Fig5C_pre_post_clonotype_overlap.csv"), row.names = FALSE)

# =============================================================================
# 6a. Fig 6A - ACT-derived clones on the pre-ACT T-cell UMAP, coloured by OS
# =============================================================================
p <- ggplot(pre_T, aes(UMAP_1, UMAP_2)) +
  geom_point(data = pre_T %>% filter(!product_clone), colour = "grey88", size = 0.2) +
  geom_point(data = pre_T %>% filter(product_clone), aes(fill = patient_OS_label),
             shape = 21, colour = "black", stroke = 0.1, size = 1.8, alpha = 0.85) +
  scale_fill_manual(values = os_cols, name = "Patient\n(OS from diagnosis)", drop = FALSE) +
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 10) +
  theme(legend.key.size = unit(0.4, "cm"))
save_plot(p, "Fig6A_UMAP_ACT_clones_pre.png", width = 7, height = 5.5)

# =============================================================================
# 6b. Fig 6B - proportion of ACT-derived clones among T cells
# -----------------------------------------------------------------------------
# The table used for the paper was calculated on CD8 T cells; the legend says
# "across all T cells". Both are exported - use the one that matches the figure.
# =============================================================================
act_prop <- bind_rows(
  tcells %>% filter(PBMC_identity == "CD8") %>% mutate(population = "CD8"),
  tcells %>% mutate(population = "All_T_cells")
) %>%
  group_by(patient_OS_label, patient_code, timepoint, population) %>%
  summarise(n_cells = n(), n_product = sum(product_clone),
            percent_product = 100 * n_product / n_cells, .groups = "drop") %>%
  left_join(pt_info, by = c("patient_OS_label", "patient_code"))
write.csv(act_prop, tab_path("Fig6B_ACT_clone_proportion.csv"), row.names = FALSE)

# =============================================================================
# 6c. Fig 6C - persistence of ACT-derived clonotypes pre vs post
# =============================================================================
cmv_clones <- meta %>%
  filter(PBMC_identity_2 %in% c("NK", "CD3"), !is.na(clone_id), product_clone) %>%
  distinct(patient_code, OS_from_diagnosis_months, timepoint, clone_id)

cmv_overlap <- cmv_clones %>%
  group_by(patient_code, OS_from_diagnosis_months) %>%
  filter(all(c("Pre", "Post") %in% timepoint)) %>%
  summarise(
    n_pre_cmv  = n_distinct(clone_id[timepoint == "Pre"]),
    n_post_cmv = n_distinct(clone_id[timepoint == "Post"]),
    n_shared   = length(intersect(clone_id[timepoint == "Pre"], clone_id[timepoint == "Post"])),
    .groups = "drop"
  ) %>%
  mutate(pct_pre_cmv_shared_with_post = 100 * n_shared / n_pre_cmv,
         pct_post_cmv_shared_with_pre = 100 * n_shared / n_post_cmv) %>%
  arrange(OS_from_diagnosis_months)
write.csv(cmv_overlap, tab_path("Fig6C_ACT_clone_persistence.csv"), row.names = FALSE)

plot_6c <- cmv_overlap %>%
  transmute(patient_code, OS_from_diagnosis_months,
            Pre_Shared = pct_pre_cmv_shared_with_post, Post_Shared = pct_post_cmv_shared_with_pre) %>%
  mutate(Pre_Unique = 100 - Pre_Shared, Post_Unique = 100 - Post_Shared) %>%
  pivot_longer(-c(patient_code, OS_from_diagnosis_months),
               names_to = c("timepoint", "status"), names_sep = "_", values_to = "percent") %>%
  mutate(timepoint = factor(timepoint, levels = c("Pre", "Post")),
         status = factor(status, levels = c("Unique", "Shared")),
         patient = factor(paste0(patient_code, "\n(", OS_from_diagnosis_months, " mo)"),
                          levels = unique(paste0(cmv_overlap$patient_code, "\n(",
                                                 cmv_overlap$OS_from_diagnosis_months, " mo)"))))
p <- ggplot(plot_6c, aes(x = timepoint, y = percent, fill = status)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.8) +
  geom_text(aes(label = ifelse(percent > 0, round(percent), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  facet_wrap(~patient, nrow = 1) +
  scale_fill_manual(values = c(Shared = "#8e7cc3", Unique = "white")) +
  labs(x = NULL, y = "ACT-derived clonotypes (%)", fill = NULL) +
  theme_classic(base_size = 11)
save_plot(p, "Fig6C_ACT_clone_persistence.png", width = 11, height = 4)

# =============================================================================
# 7a. Fig 7A / Supp 6 - subcluster distribution of ACT-derived clones
# =============================================================================
comp <- meta %>%
  filter(timepoint %in% c("Pre", "Post"), !is.na(TRB_key), !is.na(subcluster_updated)) %>%
  mutate(timepoint = factor(timepoint, levels = c("Pre", "Post"))) %>%
  count(patient_OS_label, patient_code, timepoint, product_clone_label, subcluster_updated,
        name = "n_cells") %>%
  complete(nesting(patient_OS_label, patient_code, timepoint, product_clone_label),
           subcluster_updated, fill = list(n_cells = 0)) %>%
  group_by(patient_OS_label, timepoint, product_clone_label) %>%
  mutate(total_cells = sum(n_cells), proportion = n_cells / total_cells) %>%
  ungroup() %>%
  left_join(pt_info, by = c("patient_OS_label", "patient_code"))
write.csv(comp, tab_path("Fig7A_ACT_clone_subcluster_composition.csv"), row.names = FALSE)

act_comp <- comp %>% filter(product_clone_label == "CMV-product", total_cells > 0)

p <- act_comp %>%
  group_by(timepoint, subcluster_updated) %>%
  summarise(mean_pct = 100 * mean(proportion), .groups = "drop") %>%
  ggplot(aes(x = timepoint, y = mean_pct, fill = subcluster_updated)) +
  geom_col(colour = "white", linewidth = 0.2, width = 0.7) +
  scale_fill_manual(values = subcluster_cols) +
  labs(x = NULL, y = "ACT-derived clones (%, mean of patients)", fill = "Subcluster") +
  theme_classic(base_size = 14)
save_plot(p, "Fig7A_ACT_clone_subclusters_mean.png", width = 5, height = 6)

p <- ggplot(act_comp, aes(x = timepoint, y = 100 * proportion, fill = subcluster_updated)) +
  geom_col(colour = "white", linewidth = 0.2, width = 0.8) +
  facet_wrap(~patient_OS_label, nrow = 1) +
  scale_fill_manual(values = subcluster_cols) +
  labs(x = NULL, y = "ACT-derived clones (%)", fill = "Subcluster") +
  theme_classic(base_size = 11)
save_plot(p, "SuppFig6_ACT_clone_subclusters_per_patient.png", width = 14, height = 5)

# =============================================================================
# 7b. Fig 7C input - M2 score in ACT-derived vs other T-cell clones
# =============================================================================
m2_by_clone <- tcells %>%
  filter(timepoint %in% c("Pre", "Post"), !is.na(TRB_key)) %>%
  group_by(patient_OS_label, patient_code, timepoint, product_clone_label) %>%
  summarise(mean_M2 = mean(modscore_PBMC_M2, na.rm = TRUE),
            n_cells = n(), .groups = "drop") %>%
  left_join(pt_info, by = c("patient_OS_label", "patient_code"))
write.csv(m2_by_clone, tab_path("Fig7C_M2_ACT_vs_nonACT_clones.csv"), row.names = FALSE)

saveRDS(obj, rds_path("06_final.rds"))
