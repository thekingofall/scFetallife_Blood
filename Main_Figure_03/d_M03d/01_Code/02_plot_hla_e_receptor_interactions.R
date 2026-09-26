.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "grid"))

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

root <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
plot_ready <- file.path(root, "03_PlotData")
figures <- file.path(root, "02_Figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

signals <- c(
  "HLA-E_CD94:NKG2E",
  "HLA-E_CD94:NKG2C",
  "HLA-E_CD94:NKG2A",
  "HLA-E_KLRK1",
  "HLA-E_KLRC2",
  "HLA-E_KLRC1",
  "CLEC1B_KLRB1",
  "CLEC2C_KLRB1",
  "CLEC2B_KLRB1",
  "CLEC2D_KLRB1",
  "CD99_CD99"
)
targets <- c("CX3CR1+ NK", "CXCR6+ NK", "CD56highCD16low NK")
data <- read.csv(file.path(.scf_shared, "Main_Figure_03/d_M03d/03_PlotData/M03d_diffcell_pathwayNK3lateEarly_data.csv"), check.names = FALSE, fileEncoding = "UTF-8")
source_order <- unique(data$source[order(data$source_order)])
data$source <- factor(data$source, levels = source_order)
data$target <- factor(data$target, levels = targets)
data$dataset <- factor(data$dataset, levels = c("Early", "Late"))
data$interaction_name <- factor(data$interaction_name, levels = rev(signals))
data$p_value <- factor(data$p_value, levels = c("p < 0.01", "0.01 <= p < 0.05"))
prob_limits <- range(data$prob, na.rm = TRUE)
palette_table <- read.delim(file.path(.scf_project_root, "00_Settings/palette_registry.tsv"))
probability_levels <- paste("Communication probability", c("very low", "low", "midpoint", "high", "very high"))
probability_table <- palette_table[palette_table$scope == "M03d", ]
probability_colors <- probability_table$hex[match(probability_levels, probability_table$element)]
stopifnot(!anyNA(probability_colors))


plot <- ggplot(data, aes(x = source, y = interaction_name, color = prob, size = p_value)) +
  geom_point(alpha = 1) +
  facet_grid(dataset ~ target, scales = "fixed", switch = "y") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_color_gradientn(
    colors = probability_colors,
    limits = prob_limits,
    name = "Communication\nprobability"
  ) +
  scale_size_manual(
    values = c("p < 0.01" = 4.8, "0.01 <= p < 0.05" = 3.2),
    drop = FALSE,
    name = "p-value"
  ) +
  guides(
    color = guide_colorbar(order = 1, title.position = "top", barheight = unit(30, "mm"), barwidth = unit(5, "mm")),
    size = guide_legend(order = 2, title.position = "top", override.aes = list(color = "#202020", alpha = 1))
  ) +
  labs(x = NULL, y = "Interaction\n(ligand-receptor)") +
  theme_bw(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    panel.grid.major = element_line(color = "#DEDEDE", linewidth = 0.38),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.9),
    axis.line = element_line(color = "black", linewidth = 0.75),
    axis.ticks = element_line(color = "black", linewidth = 0.7),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 10.5, color = "black"),
    axis.text.y = element_text(size = 11.5, color = "black"),
    axis.title.y = element_text(size = 13, face = "bold", color = "black", margin = margin(r = 8)),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.9),
    strip.text.x = element_text(size = 13, face = "bold", color = "black", margin = margin(4, 5, 4, 5)),
    strip.text.y.left = element_text(size = 12.5, face = "bold", color = "black", angle = 0),
    strip.placement = "outside",
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    legend.box.spacing = unit(5, "mm"),
    panel.spacing = unit(2.5, "mm"),
    plot.margin = margin(8, 10, 8, 8)
  )

stem <- "M03d_diffcell_pathwayNK3lateEarly"
ggsave(file.path(figures, paste0(stem, ".pdf")), plot = plot, device = cairo_pdf, family = "Arial", width = 16.5, height = 8.2, units = "in", limitsize = FALSE)
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(file.path(figures, paste0(stem, ".png")), plot = plot, device = ragg::agg_png, width = 16.5, height = 8.2, units = "in", dpi = 240, limitsize = FALSE, bg = "white")
} else {
  ggsave(file.path(figures, paste0(stem, ".png")), plot = plot, device = "png", type = "cairo", width = 16.5, height = 8.2, units = "in", dpi = 240, limitsize = FALSE, bg = "white")
}
