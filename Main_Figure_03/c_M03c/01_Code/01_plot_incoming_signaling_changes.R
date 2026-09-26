.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("CellChat", "circlize", "ComplexHeatmap", "dplyr", "ggplot2", "grid", "jsonlite", "tidyr"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(ggplot2); library(grid)})
root <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
figures <- file.path(root, "02_Figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
selected_signaling <- read.csv(file.path(.scf_shared, "Main_Figure_03/c_M03c/03_PlotData/M03c_signaling_changes_data.csv"), check.names = FALSE)
selected_signaling$Group <- factor(selected_signaling$Group, levels = c("CX3CR1+ NK", "CXCR6+ NK", "CD56highCD16low NK", "Naïve CD8 T", "NK T"))
pathway_levels <- unique(as.character(selected_signaling$labels))
palette_table <- read.delim(file.path(.scf_project_root, "00_Settings/palette_registry.tsv"))
pathway_table <- palette_table[palette_table$scope == "M03c", ]
pathway_colors <- setNames(pathway_table$hex, pathway_table$element)
stopifnot(all(pathway_levels %in% names(pathway_colors)))

p_c <- ggplot(selected_signaling, aes(x = reorder(labels, -incoming), y = incoming, fill = labels)) +
  geom_col(width = 0.78, color = "black", linewidth = 0.35) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pathway_colors, guide = "none") +
  labs(x = NULL, y = "Differential incoming interaction strength\n(PBMC Early vs. PBMC Late)") +
  theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.75),
    axis.ticks = element_line(color = "black", linewidth = 0.65),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11, color = "black"),
    axis.text.y = element_text(size = 11.5, color = "black"),
    axis.title.y = element_text(size = 13, color = "black", margin = margin(r = 8)),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.75),
    strip.text = element_text(size = 12.5, color = "black", margin = margin(4, 5, 4, 5)),
    panel.spacing.x = unit(4, "mm"),
    plot.margin = margin(8, 10, 8, 8)
  )

save_gg_pair <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(figures, paste0(stem, ".pdf")), plot = plot, device = cairo_pdf, family = "Arial", width = width, height = height, units = "in", limitsize = FALSE)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(file.path(figures, paste0(stem, ".png")), plot = plot, device = ragg::agg_png, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE, bg = "white")
  } else {
    ggsave(file.path(figures, paste0(stem, ".png")), plot = plot, device = "png", type = "cairo", width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE, bg = "white")
  }
}
save_gg_pair(p_c, "M03c_signaling_changes", 12, 5.6, 300)

