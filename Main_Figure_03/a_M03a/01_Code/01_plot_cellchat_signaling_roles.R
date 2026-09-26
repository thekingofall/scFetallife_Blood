.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "ggrepel", "grid", "jsonlite", "magick", "patchwork", "ragg"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")
.scf_panel_root <- normalizePath(file.path(.scf_code_dir, ".."), winslash = "/", mustWork = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(patchwork)
  library(ragg)
  library(magick)
})

table_dir <- file.path(.scf_panel_root, "03_Source_Data")
figure_dir <- file.path(.scf_panel_root, "02_Figures")
for (path in c(table_dir, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

signaling_dir <- file.path(.scf_shared, "Main_Figure_03/a_M03a")
signaling_roles_path <- file.path(.scf_shared, "Main_Figure_03/a_M03a/M03a_Signaling_Roles.csv")
color_path <- file.path(signaling_dir, "M03_Cell_Type_Colors.json")

required_inputs <- c(signaling_roles_path, color_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "))

group_order <- c("PBMC_Early", "PBMC_Late", "Liver_Early", "Thymus_Early", "Spleen_Early")
display_titles <- c(
  PBMC_Early = "PBMC\n(early)",
  PBMC_Late = "PBMC\n(late)",
  Liver_Early = "Liver",
  Thymus_Early = "Thymus",
  Spleen_Early = "Spleen"
)
focus_labels <- c(
  "DPQ T", "DPP T", "CD56highCD16low NK", "CX3CR1+ NK",
  "CXCR6+ NK", "Naïve CD8 T", "NK T"
)
size_breaks <- c(0, 250, 500, 750, 1250)
size_labels <- c("250", "500", "750", "1000")
size_values <- c("250" = 1, "500" = 4, "750" = 7, "1000" = 10)

roles <- read.csv(signaling_roles_path, stringsAsFactors = FALSE, check.names = FALSE)
required_role_columns <- c("x", "y", "labels", "Count", "group")
if (!all(required_role_columns %in% names(roles))) stop("role source-data schema mismatch")
roles <- roles[, required_role_columns]
roles$x <- as.numeric(roles$x)
roles$y <- as.numeric(roles$y)
roles$Count <- as.numeric(roles$Count)
if (any(!is.finite(as.matrix(roles[, c("x", "y", "Count")]))) || any(roles$Count < 0)) {
  stop("role data contain invalid numeric values")
}
if (anyDuplicated(roles[, c("group", "labels")])) stop("Duplicate group/cell-type role row")
if (!setequal(unique(roles$group), group_order)) stop("Unexpected five-group set")
if (any(roles$x < 0 | roles$x > 20 | roles$y < 0 | roles$y > 60)) {
  stop("Signaling-role coordinates exceed the plot axis limits")
}
if (any(roles$Count > max(size_breaks))) stop("Cell count exceeds the point-size bins")

roles$group <- factor(roles$group, levels = group_order)
roles$display_title <- unname(display_titles[as.character(roles$group)])
roles$label_selected <- roles$labels %in% focus_labels
roles$size_class <- cut(
  roles$Count,
  breaks = size_breaks,
  labels = size_labels,
  include.lowest = TRUE
)
roles$point_size_mm <- unname(size_values[as.character(roles$size_class)])
if (anyNA(roles$size_class)) stop("Unassigned point-size class")


colors_dict <- unlist(fromJSON(color_path))
missing_colors <- setdiff(unique(roles$labels), names(colors_dict))
if (length(missing_colors)) stop("Missing colors for: ", paste(missing_colors, collapse = ", "))

plot_source <- roles[, c(
  "group", "display_title", "labels", "x", "y", "Count",
  "size_class", "point_size_mm", "label_selected"
)]
plot_source$group <- as.character(plot_source$group)
write.csv(
  plot_source,
  file.path(table_dir, "Figure3a_signaling_roles_plot_source.csv"),
  row.names = FALSE
)

set.seed(1)
role_plots <- lapply(group_order, function(group) {
  data <- roles[roles$group == group, , drop = FALSE]
  label_data <- data[data$label_selected, , drop = FALSE]
  ggplot(data, aes(x = x, y = y)) +
    geom_point(aes(size = size_class, colour = labels), alpha = 1) +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = labels, colour = labels),
      size = 5,
      seed = 1,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.2,
      segment.alpha = 0.5,
      box.padding = 0.25,
      point.padding = 0.15,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = colors_dict, guide = "none", drop = FALSE) +
    scale_size_manual(
      values = size_values,
      breaks = size_labels,
      labels = size_labels,
      name = "Count",
      drop = FALSE,
      guide = guide_legend(override.aes = list(colour = "black", alpha = 1))
    ) +
    scale_x_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60), expand = c(0, 0)) +
    labs(
      title = unname(display_titles[[group]]),
      x = "Outgoing interaction strength",
      y = "Incoming interaction strength"
    ) +
    theme_classic(base_size = 10, base_family = "Arial") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 30, face = "plain", hjust = 0.5, lineheight = 0.95),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10, colour = "black"),
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.height = grid::unit(0.18, "in"),
      plot.margin = margin(8, 18, 8, 8)
    )
})
names(role_plots) <- group_order
combined <- wrap_plots(role_plots, ncol = 5, guides = "collect")

pdf_path <- file.path(figure_dir, "Figure3a_CellChat_signaling_roles.pdf")
png_path <- file.path(figure_dir, "Figure3a_CellChat_signaling_roles.png")
ggsave(
  pdf_path, plot = combined, width = 40, height = 10, units = "in",
  device = cairo_pdf, bg = "white", limitsize = FALSE
)
ggsave(
  png_path, plot = combined, width = 40, height = 10, units = "in", dpi = 300,
  device = ragg::agg_png, bg = "white", limitsize = FALSE
)


