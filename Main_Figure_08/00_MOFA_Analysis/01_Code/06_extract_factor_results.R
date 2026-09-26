.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("MOFA2", "dplyr", "tidyr"))
suppressPackageStartupMessages({library(MOFA2); library(dplyr); library(tidyr)})
root <- normalizePath(file.path(.scf_start_dir, "../../.."), winslash = "/")
input_dir <- file.path(root, "00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA")
source_dir <- file.path(root, "Main_Figure_08/00_MOFA_Analysis/02_Data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
original_model_file <- file.path(input_dir, "Fetal_Immune_Atlas_MOFA.hdf5")
original_enrichment_all <- file.path(.scf_shared, "Main_Figure_08/00_MOFA/MOFA_GO_BP_All.csv")
original_enrichment_top <- file.path(input_dir, "MOFA_GO_BP_Selected.csv")
model <- load_model(original_model_file, verbose = FALSE)
data_list <- readRDS(file.path(input_dir, "MOFA_View_Inputs.rds"))
factor_scores <- get_factors(model, factors = "all", as.data.frame = TRUE)
weights <- get_weights(model, factors = "all", views = "all", as.data.frame = TRUE)
factor_scores$sample <- as.character(factor_scores$sample)
factor_scores$factor <- as.character(factor_scores$factor)
factor_scores$pcw <- as.numeric(substring(factor_scores$sample, 2, 5))
factor_scores$factor_number <- as.integer(sub("Factor", "", factor_scores$factor))
factor_scores <- factor_scores[order(factor_scores$factor_number, factor_scores$pcw, factor_scores$sample), ]
weights$view <- as.character(weights$view)
weights$factor <- as.character(weights$factor)

sample_order <- factor_scores %>%
  distinct(sample, pcw) %>%
  arrange(pcw, sample) %>%
  pull(sample)

view_info <- data.frame(
  view = c(
    "HSC_MPP", "MEMP", "MEP", "Naive_CD4_T", "Treg", "Naive_CD8_T",
    "Gamma_Delta_V2_T", "Th17like_INNATE_T", "NK_T", "CX3CR1x_NK", "CXCR6x_NK",
    "CXCR5low_Naive_B", "CXCR5high_Naive_B", "MyeloidxCD177",
    "CD14xPPBPx_Monocytes", "Classical_Monocytes", "DC2", "pDC",
    "Cell_stimulate_BulkRNA", "Cell_stimulate_Olink", "Plama_Olink",
    "Sflow_Freq_inParent", "CD4TCR", "CD8TCR", "BCR"
  ),
  label = c(
    "HSC/MPP", "MEMP", "MEP", "Naïve CD4+ T", "Treg", "Naïve CD8+ T",
    "Vδ2", "Th17-like Innate T", "NKT", "CX3CR1+ NK", "CXCR6+ NK",
    "CXCR5− Naïve B", "CXCR5+ Naïve B", "Myeloid−CD177",
    "CD14+PPBP+ Monocytes", "Classical Monocytes", "DC2", "pDC",
    "Stimulation bulk RNA-seq", "Stimulation proteomics", "Plasma proteomics",
    "Spectral flow cytometry", "CD4 TCR", "CD8 TCR", "BCR"
  ),
  modality = c(
    rep("Single-cell RNA-seq", 18), "Stimulation bulk RNA-seq",
    "Olink", "Olink", "Spectral flow cytometry", "TCR", "TCR", "BCR"
  )
)
if (!setequal(view_info$view, unique(weights$view)) || !setequal(view_info$view, names(data_list))) {
  stop("View mapping does not match the MOFA model and input data list")
}

view_order <- view_info$view
label_order <- view_info$label
view_to_label <- setNames(view_info$label, view_info$view)

clean_feature_name <- function(feature, view) {
  output <- sub(".*__", "", feature)
  if (view %in% c("Cell_stimulate_Olink", "Plama_Olink")) {
    output <- sub(paste0("_", view, ".*$"), "", output)
  }
  if (view == "Cell_stimulate_BulkRNA") {
    output <- sub("_Cell_stimulate_BulkRNA$", "", output)
  }
  output <- sub("_CD4TCR$", "", output)
  output <- sub("_CD8TCR$", "", output)
  output <- gsub("CXCR5.B.1", "CXCR5+ B", output, fixed = TRUE)
  output <- gsub("CXCR5.B", "CXCR5− B", output, fixed = TRUE)
  output <- gsub("Inhibitory.CD8.T", "Inhibitory CD8+T", output, fixed = TRUE)
  if (view %in% c("CD4TCR", "CD8TCR") && grepl("^CD[48]_CDR3b_length_", output)) {
    output <- gsub("_", "\u2212", output, fixed = TRUE)
  }
  output
}

variance <- plot_variance_explained(model, x = "view", y = "factor")$data %>%
  dplyr::select(-any_of("group")) %>%
  mutate(
    view = as.character(view), factor = as.character(factor),
    factor_number = as.integer(sub("Factor", "", factor))
  )
factor_levels <- paste0("Factor", sort(unique(variance$factor_number)))

variance_source <- variance %>%
  left_join(view_info, by = "view") %>%
  arrange(match(view, view_order), factor_number)
write.csv(variance_source, file.path(root, "Main_Figure_08/b_M08b/03_PlotData/M08b_variance_explained_plot_data.csv"), row.names = FALSE)

feature_counts <- data.frame(
  view = names(data_list),
  feature_count = vapply(data_list, nrow, integer(1))
)
view_summary <- variance %>%
  group_by(view) %>%
  summarise(total_variance = min(sum(value, na.rm = TRUE), 100), .groups = "drop") %>%
  left_join(feature_counts, by = "view") %>%
  left_join(view_info, by = "view") %>%
  arrange(match(view, view_order))
write.csv(view_summary, file.path(root, "Main_Figure_08/b_M08b/03_PlotData/M08b_view_summary.csv"), row.names = FALSE)

age_stats <- factor_scores %>%
  group_by(factor, factor_number) %>%
  summarise(
    rho = cor(value, pcw, method = "spearman"),
    p = cor.test(value, pcw, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  arrange(rho)
write.csv(age_stats, file.path(root, "Main_Figure_08/c_M08c/03_PlotData/M08c_factor_age_spearman_statistics.csv"), row.names = FALSE)

factor3_data <- factor_scores %>% filter(factor == "Factor3") %>% arrange(pcw, sample)
factor3_stat <- age_stats %>% filter(factor == "Factor3")
if (!nrow(factor3_data) || nrow(factor3_stat) != 1L) stop("Factor3 source data are incomplete")
write.csv(factor3_data, file.path(root, "Main_Figure_08/d_M08d/03_PlotData/M08d_Factor3_scores.csv"), row.names = FALSE)

factor3_weights <- weights %>% filter(factor == "Factor3")
top2_weights <- bind_rows(
  factor3_weights %>% filter(value > 0) %>% group_by(view) %>% slice_max(value, n = 2, with_ties = FALSE) %>% mutate(sign = "positive"),
  factor3_weights %>% filter(value < 0) %>% group_by(view) %>% slice_min(value, n = 2, with_ties = FALSE) %>% mutate(sign = "negative")
) %>%
  ungroup() %>%
  mutate(
    view_label = view_to_label[view],
    feature_label = mapply(clean_feature_name, feature, view, USE.NAMES = FALSE)
  ) %>%
  arrange(match(view, view_order), sign, desc(value))
write.csv(top2_weights, file.path(source_dir, "Figure8e_factor3_top2_features.csv"), row.names = FALSE)

get_view_matrix <- function(view_name) {
  view_data <- model@data[[view_name]]
  if (is.list(view_data)) view_data <- view_data[[1]]
  as.matrix(view_data)
}

heat_rows <- list()
for (i in seq_len(nrow(top2_weights))) {
  row <- top2_weights[i, ]
  view_name <- as.character(row$view[[1]])
  feature_name <- as.character(row$feature[[1]])
  view_label <- as.character(row$view_label[[1]])
  sign_name <- as.character(row$sign[[1]])
  feature_label <- as.character(row$feature_label[[1]])
  matrix_data <- get_view_matrix(view_name)
  feature_index <- match(feature_name, rownames(matrix_data))
  sample_indices <- match(sample_order, colnames(matrix_data))
  if (is.na(feature_index)) {
    stop("Selected feature is absent from model data: ", view_name, " / ", feature_name)
  }
  if (any(is.na(sample_indices))) {
    stop("Sample order is incomplete in model data for view ", view_name)
  }
  heat_rows[[i]] <- data.frame(
    view = view_name,
    view_label = view_label,
    sign = sign_name,
    feature = feature_name,
    feature_label = feature_label,
    feature_key = paste(view_name, sign_name, feature_name, sep = "___"),
    sample = sample_order,
    pcw = as.numeric(substring(sample_order, 2, 5)),
    expression_value = as.numeric(matrix_data[feature_index, sample_indices, drop = TRUE]),
    factor3_weight = as.numeric(row$value[[1]])
  )
}
heat_data <- bind_rows(heat_rows) %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    view_label = factor(view_label, levels = label_order),
    sign = factor(sign, levels = c("positive", "negative"))
  )
feature_key_levels <- heat_data %>%
  distinct(view, sign, feature_key) %>%
  arrange(match(view, view_order), sign) %>%
  pull(feature_key)
heat_data$feature_key <- factor(heat_data$feature_key, levels = rev(unique(feature_key_levels)))
write.csv(heat_data, file.path(root, "Main_Figure_08/e_M08e/03_PlotData/M08e_Factor3_top2_heatmap_values.csv"), row.names = FALSE)

all_enrichment <- read.csv(original_enrichment_all, check.names = FALSE)
top_enrichment <- read.csv(original_enrichment_top, check.names = FALSE)
top_enrichment$sign <- sub("_.*$", "", top_enrichment$Group)
top_enrichment$view <- sub("^(positive|negative)_", "", top_enrichment$Group)
top_enrichment <- top_enrichment %>%
  filter(view %in% view_order, sign %in% c("positive", "negative")) %>%
  mutate(view_label = view_to_label[view]) %>%
  arrange(sign, match(view, view_order), p.adjust)
write.csv(top_enrichment, file.path(root, "Main_Figure_08/f_M08f/03_PlotData/M08f_Factor3_top1_GO_BP_by_view.csv"), row.names = FALSE)

s18_rows <- list()
for (view_name in view_order) {
  view_weights <- factor3_weights %>% filter(view == view_name)
  denominator <- max(abs(view_weights$value), na.rm = TRUE)
  if (!is.finite(denominator) || denominator <= 0) stop("Invalid Factor3 weight scale for ", view_name)
  selected <- bind_rows(
    view_weights %>% filter(value < 0) %>% slice_min(value, n = 6, with_ties = FALSE) %>% mutate(sign = "negative"),
    view_weights %>% filter(value > 0) %>% slice_max(value, n = 6, with_ties = FALSE) %>% mutate(sign = "positive")
  ) %>%
    mutate(
      normalized_weight = value / denominator,
      feature_label = mapply(clean_feature_name, feature, view, USE.NAMES = FALSE),
      view_label = view_to_label[view]
    )
  s18_rows[[view_name]] <- selected
}
s18_data <- bind_rows(s18_rows) %>% arrange(match(view, view_order), normalized_weight)
write.csv(s18_data, file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_18/a_S18a/03_PlotData/S18a_Factor3_Top6_Weights_source.csv"), row.names = FALSE)

write.csv(factor_scores, file.path(source_dir, "MOFA_factor_scores.csv"), row.names = FALSE)
write.csv(factor3_weights, file.path(source_dir, "MOFA_factor3_weights.csv"), row.names = FALSE)
write.csv(factor3_stat, file.path(root, "Main_Figure_08/d_M08d/03_PlotData/M08d_Factor3_statistics.csv"), row.names = FALSE)
