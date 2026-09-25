###############################################################################
# 02_TCR_preprocessing.R
#
#   1. Collect filtered_contig_annotations.csv from every sample (vdj_t)
#   2. Collapse contigs to one row per cell; resolve cells with two TRB or two
#      TRA chains by keeping the chain with the most UMIs (ties are removed);
#      drop cells without a TRB chain (TRB is used for clonotype calling)
#   3. Add the cleaned TCR table and sequencing batch to the Seurat object
#
# Input : data/raw/GBM*/*/outs/per_sample_outs/<sample>/vdj_t/filtered_contig_annotations.csv
#         results/rds/01_merged_qc.rds
# Output: results/tables/GBM_TCR_contigs.csv, results/tables/GBM_TCR_clean.csv
#         results/rds/02_merged_qc_tcr.rds
###############################################################################

source("R/helpers.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
})

merged <- readRDS(rds_path("01_merged_qc.rds"))
sample_names <- unique(merged$orig.ident)

# =============================================================================
# 1. Collect VDJ contigs
# =============================================================================
gem_dirs <- list.dirs(dirs$raw, recursive = FALSE, full.names = TRUE)
gem_dirs <- gem_dirs[grepl("^GBM", basename(gem_dirs))]

vdj_entries <- list()
for (gem_dir in gem_dirs) {
  vdj_t_dirs <- Sys.glob(file.path(gem_dir, "*", "outs", "per_sample_outs", "*", "vdj_t"))
  for (vdj_dir in vdj_t_dirs) {
    sample_name <- basename(dirname(vdj_dir))
    if (!sample_name %in% sample_names) next
    vdj_file <- file.path(vdj_dir, "filtered_contig_annotations.csv")
    if (!file.exists(vdj_file)) next
    dat <- read_csv(vdj_file, show_col_types = FALSE)
    if (nrow(dat) == 0) next
    vdj_entries[[paste(sample_name, vdj_dir)]] <- dat %>%
      mutate(barcode = paste0(sample_name, "_", barcode),   # matches Seurat cell names
             id      = sample_name,
             batch   = basename(gem_dir))                   # GEM run / sequencing batch
  }
}
df <- bind_rows(vdj_entries)
write.csv(df, tab_path("GBM_TCR_contigs.csv"), row.names = FALSE)

# =============================================================================
# 2. One row per cell, resolve double chains
# =============================================================================
tcr_columns  <- c("cdr3", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt")
meta_columns <- setdiff(colnames(df), c(tcr_columns, "chain", "barcode", "contig_id",
                                        "umis", "reads", "length", "productive"))

df_wide <- df %>%
  filter(chain %in% c("TRA", "TRB")) %>%
  mutate(across(all_of(tcr_columns), as.character)) %>%
  group_by(barcode, chain) %>%
  mutate(chain_label = paste0(chain, row_number())) %>%   # TRA1, TRA2, TRB1, TRB2
  ungroup() %>%
  select(barcode, chain_label, all_of(tcr_columns), umis, reads) %>%
  pivot_wider(names_from = chain_label,
              values_from = c(all_of(tcr_columns), umis, reads),
              names_glue = "{chain_label}_{.value}")

df_collapsed <- df_wide %>%
  left_join(df %>% select(barcode, all_of(meta_columns)) %>% distinct(barcode, .keep_all = TRUE),
            by = "barcode")

# For cells with two chains of the same type keep the chain with more UMIs;
# drop the cell if UMIs are tied (true chain cannot be identified).
resolve_double_chain <- function(d, chain) {
  c1 <- grep(paste0("^", chain, "1_"), colnames(d), value = TRUE)
  c2 <- grep(paste0("^", chain, "2_"), colnames(d), value = TRUE)
  if (length(c2) == 0) return(d)
  cdr3_1 <- paste0(chain, "1_cdr3"); cdr3_2 <- paste0(chain, "2_cdr3")
  umi_1  <- paste0(chain, "1_umis"); umi_2  <- paste0(chain, "2_umis")

  double <- !is.na(d[[cdr3_1]]) & !is.na(d[[cdr3_2]])
  tie    <- double & (is.na(d[[umi_1]]) | is.na(d[[umi_2]]) | d[[umi_1]] == d[[umi_2]])
  use_2  <- double & !tie & d[[umi_2]] > d[[umi_1]]

  d[use_2, c1] <- d[use_2, c2]     # move chain 2 into slot 1 (whole block, genes stay aligned)
  d[double & !tie, c2] <- NA       # blank slot 2
  message(chain, ": ", sum(double), " cells with 2 chains, ", sum(tie), " removed (tied UMIs)")
  d[!tie, ]
}

df_clean <- df_collapsed %>%
  resolve_double_chain("TRB") %>%
  filter(!is.na(TRB1_cdr3)) %>%                 # TRB required for clonotype calling
  resolve_double_chain("TRA")

# Collapse to single TRA_* / TRB_* columns
for (chain in c("TRA", "TRB")) {
  for (f in c(tcr_columns, "umis", "reads")) {
    c1 <- paste0(chain, "1_", f); c2 <- paste0(chain, "2_", f)
    v1 <- if (c1 %in% colnames(df_clean)) df_clean[[c1]] else NA
    v2 <- if (c2 %in% colnames(df_clean)) df_clean[[c2]] else NA
    df_clean[[paste0(chain, "_", f)]] <- coalesce(v1, v2)
  }
}
df_clean <- df_clean %>%
  select(-matches("^TR[AB][12]_"), -any_of("sample")) %>%
  filter(!is.na(TRB_cdr3))

write.csv(df_clean, tab_path("GBM_TCR_clean.csv"), row.names = FALSE)

# =============================================================================
# 3. Add TCR + batch to the Seurat object
# -----------------------------------------------------------------------------
# NOTE: the original script for this step (which produced
# 8_VDJ_merged_metadata_all_samples.rds) was not among the files used to build
# this repository. This is a reconstruction: a left join on the cell barcode,
# and the GEM batch assigned per sample so that every cell (with or without a
# TCR) has a batch for the Harmony integration in 03_integration.R.
# =============================================================================
tcr_meta <- df_clean %>%
  select(barcode, starts_with("TRA_"), starts_with("TRB_")) %>%
  distinct(barcode, .keep_all = TRUE) %>%
  as.data.frame()
rownames(tcr_meta) <- tcr_meta$barcode
merged <- AddMetaData(merged, tcr_meta[intersect(colnames(merged), rownames(tcr_meta)), -1])

sample_batch <- df %>% distinct(id, batch)
stopifnot(!any(duplicated(sample_batch$id)))
merged$batch <- sample_batch$batch[match(merged$id, sample_batch$id)]
merged$barcode <- colnames(merged)

message("Cells with TRB: ", sum(!is.na(merged$TRB_cdr3)), " / ", ncol(merged))
table(merged$batch, useNA = "ifany")

saveRDS(merged, rds_path("02_merged_qc_tcr.rds"))
