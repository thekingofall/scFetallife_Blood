.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "magrittr", "nichenetr", "tibble", "tidyr"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(magrittr)
  library(nichenetr)
})

panel_dir <- normalizePath(file.path(.scf_start_dir, ".."))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0L && length(args) != 2L) stop("Usage: Rscript 01_plot_nk_ligand_receptor_heatmap.R [matrix_csv output_prefix]")
matrix_csv <- if (length(args)) args[[1]] else file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_08/03_PlotData/S08a_CX3CR1_NK_Interaction_Matrix.csv")
output_prefix <- if (length(args)) args[[2]] else file.path(panel_dir, "02_Figures/S08a_CX3CR1_ligand_receptor_heatmap")

matrix_df <- read.csv(
  matrix_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)
if (!"receptor" %in% colnames(matrix_df)) {
  stop("matrix CSV lacks required receptor column: ", matrix_csv)
}
if (ncol(matrix_df) < 2L || nrow(matrix_df) < 1L) {
  stop("matrix CSV is empty or has no ligand columns: ", matrix_csv)
}

excluded_ligands <- c("HLA-A", "HLA-G")
missing_excluded_ligands <- setdiff(excluded_ligands, colnames(matrix_df))
if (length(missing_excluded_ligands) > 0L) {
  stop("requested ligand columns are absent: ", paste(missing_excluded_ligands, collapse = ", "))
}
matrix_df <- matrix_df[, !colnames(matrix_df) %in% excluded_ligands, drop = FALSE]

receptors <- matrix_df$receptor
lr_matrix <- as.matrix(matrix_df[, setdiff(colnames(matrix_df), "receptor"), drop = FALSE])
storage.mode(lr_matrix) <- "double"
rownames(lr_matrix) <- receptors

heatmap_plot <- make_heatmap_ggplot(
  t(lr_matrix),
  y_name = "Ligands",
  x_name = "Receptors",
  color = "mediumvioletred",
  legend_title = "Prior interaction potential",
  legend_position = "right",
  x_axis_position = "top"
) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x.top = element_text(angle = 90, hjust = 0, size = 7),
    axis.text.y = element_text(size = 7)
  )

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
pdf_path <- paste0(output_prefix, ".pdf")
png_path <- paste0(output_prefix, ".png")

ggsave(
  pdf_path,
  heatmap_plot,
  width = 8,
  height = 6,
  units = "in",
  device = cairo_pdf
)
ggsave(
  png_path,
  heatmap_plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

