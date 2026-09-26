.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggforce", "ggplot2", "ggrastr", "ggrepel", "grid", "jsonlite", "ragg"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(
  Sys.getenv("SCF_PROJECT_ROOT", unset = file.path(.scf_code_dir, "../../..")),
  winslash = "/",
  mustWork = FALSE
)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(jsonlite)
  library(grid)
})

repo <- .scf_project_root
runner_dir <- Sys.getenv(
  "SCF_RASTER_RUNNER_DIR",
  unset = .scf_code_dir
)
output_root <- Sys.getenv(
  "SCF_RASTER_OUTPUT_ROOT",
  unset = .scf_project_root
)
original_umap <- file.path(.scf_shared, "Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv")
original_colors <- file.path(.scf_shared, "Main_Figure_01/01_Cell_Annotation/M01_Cell_Type_Colors.json")
source(file.path(.scf_code_dir, "01_define_umap_plot_functions.R"))
source(file.path(.scf_code_dir, "02_define_umap_arrow_layer.R"))

umap_data <- read.csv(original_umap, row.names = 1)

colors_dict <- unlist(fromJSON(original_colors))

celltype_dir <- file.path(output_root, "Main_Figure_01/c_M01c/02_Figures")
organ_dir <- file.path(output_root, "Main_Figure_01/d_M01d/02_Figures")
dir.create(celltype_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(organ_dir, recursive = TRUE, showWarnings = FALSE)

umap_plot <- create_umap_plot(umap_data = umap_data, colors_dict = colors_dict) +
  geom_markArrow(
    data = umap_data, rel.pos = 0.06, rel.len = 0.1,
    label.shift = c(0.025, -0.025)
  )
ggsave(
  file.path(celltype_dir, "M01c_celltype_UMAP.pdf"),
  plot = umap_plot, width = 22, height = 12,
  device = cairo_pdf, limitsize = FALSE
)
ggsave(
  file.path(celltype_dir, "M01c_celltype_UMAP.png"),
  plot = umap_plot, width = 22, height = 12,
  device = ragg::agg_png, dpi = 180, bg = "white", limitsize = FALSE
)

umap_data$organ <- factor(umap_data$organ, levels = c("PBMC", "Liver", "Thymus", "Spleen"))
umap_organ <- ggplot() +
  ggrastr::geom_point_rast(
    data = umap_data,
    aes(x = UMAP1, y = UMAP2, color = organ),
    size = 0.5, alpha = 0.8, key_glyph = draw_number_circle,
    raster.dpi = 300
  ) +
  scale_color_manual(values = c("#C71000FF", "#f49128", "#023f75", "#5A9599FF")) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 10), ncol = 1)) +
  geom_markArrow(data = umap_data, rel.pos = 0.06, rel.len = 0.1) +
  theme_minimal() +
  theme(
    panel.spacing.y = unit(0, "mm"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) +
  labs(color = "")
ggsave(
  file.path(organ_dir, "M01d_organ_UMAP.pdf"),
  plot = umap_organ, width = 9, height = 9,
  device = cairo_pdf, limitsize = FALSE
)
ggsave(
  file.path(organ_dir, "M01d_organ_UMAP.png"),
  plot = umap_organ, width = 9, height = 9,
  device = ragg::agg_png, dpi = 320, bg = "white", limitsize = FALSE
)

cat("M01c-d raster-point PDFs complete; retained cells=", nrow(umap_data), "\n", sep = "")
