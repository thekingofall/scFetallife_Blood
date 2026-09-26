.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "grid", "patchwork", "ragg", "readr", "tidyr"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(readr)
  library(tidyr)
})

repo <- .scf_project_root
input_dir <- file.path(.scf_start_dir, "../03_PlotData")

save_pair <- function(plot, directory, stem, width, height, dpi = 320) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(directory, paste0(stem, ".png")), plot, width = width, height = height, device = ragg::agg_png, dpi = dpi, bg = "white", limitsize = FALSE)
}


variance <- read_csv(file.path(input_dir, "M08b_variance_explained_plot_data.csv"), show_col_types = FALSE)
view_summary <- read_csv(file.path(input_dir, "M08b_view_summary.csv"), show_col_types = FALSE)
view_info <- variance %>% distinct(view, label, modality)
view_order <- view_info$view
label_order <- view_info$label
factor_order <- paste0("Factor", sort(unique(variance$factor_number)))

palette <- read.delim(file.path(.scf_project_root, "00_Settings/palette_registry.tsv"))
palette <- subset(palette, scope == "M08_views")
view_colors <- setNames(palette$hex, palette$element)
modality_palette <- read.delim(file.path(.scf_project_root, "00_Settings/palette_registry.tsv"))
modality_palette <- subset(modality_palette, scope == "M08_modality")
modality_colors <- setNames(modality_palette$hex, modality_palette$element)

variance <- variance %>%
  mutate(
    label = factor(label, levels = rev(label_order)),
    factor = factor(factor, levels = factor_order)
  )
view_info <- view_info %>% mutate(label = factor(label, levels = rev(label_order)))
view_summary <- view_summary %>% mutate(label = factor(label, levels = rev(label_order)))

p8b_type <- ggplot(view_info, aes(1, label, fill = modality)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = modality_colors, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "View", y = NULL) +
  theme_minimal(base_size = 11, base_family = "Arial") +
  theme(
    panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_text(color = "black", size = 11),
    panel.border = element_blank(),
    plot.margin = margin(4, 0, 24, 4)
  )

p8b_heat <- ggplot(variance, aes(factor, label, fill = value)) +
  geom_tile(color = "black", linewidth = 0.20) +
  scale_fill_gradientn(
    colours = c("#FAFCFB", "#CBE7E2", "#99D7CA", "#54C3C4", "#14AFBE", "#0092AE", "#00789F", "#183D59", "#334555"),
    name = "Variance (%)",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
  ) +
  scale_x_discrete(labels = function(x) sub("Factor", "", x)) +
  labs(
    x = "Factor", y = NULL,
    title = "Percentage variance in each\ndata type (view) explained by each factor"
  ) +
  theme_minimal(base_size = 11, base_family = "Arial") +
  theme(
    panel.grid = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
    panel.border = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "plain"),
    legend.text = element_text(size = 11),
    plot.margin = margin(4, 2, 4, 0)
  )

p8b_total <- ggplot(view_summary, aes(total_variance, label, fill = view)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = view_colors, guide = "none") +
  labs(
    x = "Variance (%)", y = NULL,
    title = "Percentage of\nexplained variance per view"
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  )

p8b_count <- ggplot(view_summary, aes(feature_count, label, fill = view)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = view_colors, guide = "none") +
  labs(
    x = "Number of features", y = NULL,
    title = "Feature count\nin each view"
  ) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  )

p8b <- wrap_plots(
  p8b_type, p8b_heat, p8b_total, p8b_count, guide_area(),
  widths = c(0.11, 2.4, 1.216, 1.216, 0.55),
  guides = "collect", nrow = 1
)
out8b <- file.path(.scf_project_root, "Main_Figure_08", "b_M08b", "02_Figures")
save_pair(p8b, out8b, "M08b_variance_explained", 10, 5, dpi = 300)

