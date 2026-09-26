.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "grDevices", "grid"))

panel_root <- normalizePath(file.path(.scf_start_dir,".."),mustWork=TRUE)
input_path <- file.path(
  panel_root,
  "03_PlotData",
  "S11d_TCR_top15_plot_source_data.csv"
)
output_dir <- file.path(panel_root, "02_Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("missing R package: ggplot2")
}
suppressPackageStartupMessages(library(ggplot2))

plot_data <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("clone_id", "MainID", "Last_cell_type", "cells")
if (!all(required_columns %in% names(plot_data))) {
  stop("input is missing required columns")
}

plot_data$clone_id <- as.character(plot_data$clone_id)
plot_data$cells <- as.integer(plot_data$cells)

sample_levels <- c(
  "B11.6_P24", "B22.4_P9", "B37.9_P20", "T10.0_P1",
  "T10.1_P25", "T18.6_P5", "T24.6_P11"
)
cell_levels <- c(
  "DP(P) T", "DN(Q) T", "Treg", "Naïve CD4 T",
  "Th17like_INNATE_T", "Cycling Treg"
)
clone_levels <- unique(plot_data$clone_id)


plot_data$MainID <- factor(plot_data$MainID, levels = sample_levels)
plot_data$clone_id <- factor(plot_data$clone_id, levels = clone_levels)
plot_data$Last_cell_type <- factor(plot_data$Last_cell_type, levels = cell_levels)

display_data <- plot_data
display_data$plot_x <- as.character(display_data$clone_id)
spacer_data <- do.call(rbind, lapply(sample_levels, function(sample_id) {
  sample_clone_count <- length(unique(plot_data$clone_id[plot_data$MainID == sample_id]))
  left_count <- if (sample_clone_count == 1L) 2L else 1L
  right_count <- if (sample_clone_count == 1L || (sample_clone_count == 2L && nchar(sample_id) >= 9L)) 2L else 1L
  data.frame(
    clone_id = NA_character_,
    MainID = sample_id,
    Last_cell_type = "DP(P) T",
    cells = 0L,
    plot_x = c(
      paste0(sample_id, "__spacer_left_", seq_len(left_count)),
      paste0(sample_id, "__spacer_right_", seq_len(right_count))
    ),
    stringsAsFactors = FALSE
  )
}))
display_data$clone_id <- as.character(display_data$clone_id)
display_data$MainID <- as.character(display_data$MainID)
display_data$Last_cell_type <- as.character(display_data$Last_cell_type)
display_data <- rbind(display_data, spacer_data)
display_data$MainID <- factor(display_data$MainID, levels = sample_levels)
display_data$Last_cell_type <- factor(display_data$Last_cell_type, levels = cell_levels)
plot_x_levels <- unlist(lapply(sample_levels, function(sample_id) {
  sample_clones <- clone_levels[
    clone_levels %in% as.character(plot_data$clone_id[plot_data$MainID == sample_id])
  ]
  left_count <- if (length(sample_clones) == 1L) 2L else 1L
  right_count <- if (length(sample_clones) == 1L || (length(sample_clones) == 2L && nchar(sample_id) >= 9L)) 2L else 1L
  c(
    paste0(sample_id, "__spacer_left_", seq_len(left_count)),
    sample_clones,
    paste0(sample_id, "__spacer_right_", seq_len(right_count))
  )
}), use.names = FALSE)
display_data$plot_x <- factor(display_data$plot_x, levels = plot_x_levels)

cell_palette <- c(
  "DP(P) T" = "#E06E3A",
  "DN(Q) T" = "#C73934",
  "Treg" = "#008EA0",
  "Naïve CD4 T" = "#5A9599",
  "Th17like_INNATE_T" = "#DE624E",
  "Cycling Treg" = "#88C9D6"
)
cell_labels <- c(
  "DP(P) T" = "DP(P) T",
  "DN(Q) T" = "DN(Q) T",
  "Treg" = "Treg",
  "Naïve CD4 T" = "Naïve CD4 T",
  "Th17like_INNATE_T" = "Th17-like innate T",
  "Cycling Treg" = "Cycling Treg"
)

figure <- ggplot(
  display_data,
  aes(x = plot_x, y = cells, fill = Last_cell_type)
) +
  geom_col(width = 0.82, linewidth = 0) +
  facet_grid(
    cols = vars(MainID),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = cell_palette,
    breaks = cell_levels,
    labels = unname(cell_labels[cell_levels]),
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = function(x) ifelse(grepl("__spacer_", x), "", x),
    expand = expansion(add = c(0.15, 0.15))
  ) +
  scale_y_continuous(
    breaks = seq(0, 18, 3),
    limits = c(0, 18.8),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Clone ID",
    y = "Clone number of cells\nwith the same clonotype",
    fill = "Cell type"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0,
      byrow = TRUE,
      override.aes = list(alpha = 1)
    )
  ) +
  theme_bw(base_size = 18, base_family = "") +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    panel.border = element_rect(
      fill = NA,
      colour = "#181818",
      linewidth = 0.9
    ),
    panel.grid.major = element_line(colour = "#EBEBEB", linewidth = 0.55),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(0.08, "cm"),
    strip.background = element_rect(
      fill = "#DAD9D9",
      colour = "#181818",
      linewidth = 0.9
    ),
    strip.text.x = element_text(
      colour = "#181818",
      size = 19.5,
      face = "bold",
      margin = margin(t = 5, r = 4, b = 5, l = 4)
    ),
    axis.text.x = element_text(
      colour = "#181818",
      size = 18,
      angle = 52,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(colour = "#181818", size = 18),
    axis.title.x = element_text(colour = "#181818", size = 23, face = "bold", margin = margin(t = 9)),
    axis.title.y = element_text(colour = "#181818", size = 23, face = "bold", lineheight = 0.95, margin = margin(r = 10)),
    axis.ticks = element_line(colour = "#181818", linewidth = 0.65),
    axis.ticks.length = grid::unit(2.5, "pt"),
    legend.position = "right",
    legend.justification = "center",
    legend.title = element_text(size = 20, face = "bold", colour = "#181818"),
    legend.text = element_text(size = 18.5, colour = "#181818"),
    legend.key.size = grid::unit(22, "pt"),
    legend.spacing.y = grid::unit(3, "pt"),
    legend.background = element_blank(),
    plot.margin = margin(t = 7, r = 8, b = 5, l = 5)
  )

pdf_path <- file.path(output_dir, "S11d_TCR_top15_clones_by_7_samples.pdf")
png_path <- file.path(output_dir, "S11d_TCR_top15_clones_by_7_samples.png")

grDevices::cairo_pdf(pdf_path, width = 15.6, height = 4.6, family = "Arial")
print(figure)
grDevices::dev.off()

grDevices::png(
  png_path,
  width = 15.6,
  height = 4.6,
  units = "in",
  res = 250,
  type = "cairo",
  family = "Arial",
  bg = "white"
)
print(figure)
grDevices::dev.off()

