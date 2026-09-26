.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "patchwork"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

panel_dir <- normalizePath(file.path(.scf_start_dir, ".."))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0L && length(args) != 2L) stop("Usage: Rscript 01_plot_mofa_factor_feature_weights.R [source_csv output_directory]")
source_csv <- if (length(args)) args[[1]] else file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_18/a_S18a/03_PlotData/S18a_Factor3_Top6_Weights_source.csv")
output_dir <- if (length(args)) args[[2]] else file.path(panel_dir, "02_Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

view_order <- c(
  "HSC_MPP", "MEMP", "MEP", "Naive_CD4_T", "Treg", "Naive_CD8_T",
  "Gamma_Delta_V2_T", "Th17like_INNATE_T", "NK_T", "CX3CR1x_NK", "CXCR6x_NK",
  "CXCR5low_Naive_B", "CXCR5high_Naive_B", "MyeloidxCD177",
  "CD14xPPBPx_Monocytes", "Classical_Monocytes", "DC2", "pDC",
  "Cell_stimulate_BulkRNA", "Cell_stimulate_Olink", "Plama_Olink",
  "Sflow_Freq_inParent", "CD4TCR", "CD8TCR", "BCR"
)
view_labels <- c(
  "HSC/MPP", "MEMP", "MEP", "Naïve CD4+ T", "Treg", "Naïve CD8+ T",
  "Vδ2", "Th17-like Innate T", "NKT", "CX3CR1+ NK", "CXCR6+ NK",
  "CXCR5− Naïve B", "CXCR5+ Naïve B", "Myeloid−CD177",
  "CD14+PPBP+ Monocytes", "Classical Monocytes", "DC2", "pDC",
  "Stimulation bulk RNA-seq", "Stimulation proteomics", "Plasma proteomics",
  "Spectral flow cytometry", "CD4 TCR", "CD8 TCR", "BCR"
)
view_to_label <- setNames(view_labels, view_order)
mop_cells_palette <- c(
  "#E41A1C", "#FF7F00", "#FFD92F", "#E6AB02", "#FFFF99", "#BC80BD",
  "#1B9E77", "#66C2A5", "#8DA0CB", "#A6CEE3", "#B3E2CD", "#FFFFB3",
  "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
  "#BC80BD", "#CCEBC5", "#377EB8", "#11A579", "#F2B701", "#66C5CC", "#80BA5A"
)
view_colors <- setNames(mop_cells_palette, view_order)

s18_data <- read.csv(source_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("view", "feature_label", "normalized_weight")
if (!all(required_columns %in% names(s18_data))) {
  stop("S18 source table lacks required columns")
}
if (any(!is.finite(s18_data$normalized_weight)) || max(abs(s18_data$normalized_weight)) > 1 + 1e-12) {
  stop("S18 normalized weights are invalid")
}

s18_plots <- lapply(view_order, function(view_name) {
  plot_data <- s18_data %>% filter(view == view_name) %>% arrange(normalized_weight)
  plot_data$feature_label <- factor(plot_data$feature_label, levels = plot_data$feature_label)
  ggplot(plot_data, aes(feature_label, normalized_weight)) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_col(fill = view_colors[[view_name]], width = 0.92) +
    geom_text(
      aes(y = ifelse(normalized_weight > 0, -0.02, 0.02), label = feature_label),
      angle = 90, hjust = ifelse(plot_data$normalized_weight > 0, 1, 0), vjust = 0.5,
      size = 3.8, family = "Arial"
    ) +
    scale_y_continuous(limits = c(-1.08, 1.08), breaks = c(-1, -0.5, 0, 0.5, 1)) +
    labs(x = NULL, y = NULL, title = view_to_label[[view_name]]) +
    theme_classic(base_size = 10.5, base_family = "Arial") +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.margin = margin(3, 3, 3, 3)
    )
})
p_s18 <- wrap_plots(s18_plots, ncol = 5)

pdf_path <- file.path(output_dir, "S18a_Factor3_Top6_Weights.pdf")
png_path <- file.path(output_dir, "S18a_Factor3_Top6_Weights.png")
ggsave(pdf_path, p_s18, width = 13.333, height = 12, units = "in", device = cairo_pdf, limitsize = FALSE)
ggsave(png_path, p_s18, width = 13.333, height = 12, units = "in", dpi = 300, bg = "white", limitsize = FALSE)
