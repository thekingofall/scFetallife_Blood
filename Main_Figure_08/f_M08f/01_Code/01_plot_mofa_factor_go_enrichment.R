.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "patchwork", "ragg", "readr"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

code_dir <- .scf_start_dir
panel_dir <- normalizePath(file.path(code_dir, ".."), winslash = "/")
input_file <- file.path(panel_dir, "03_PlotData", "M08f_Factor3_top1_GO_BP_by_view.csv")
output_dir <- file.path(panel_dir, "02_Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_pair <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = width, height = height, device = ragg::agg_png, dpi = dpi, bg = "white", limitsize = FALSE)
}

view_order <- read_csv(file.path(panel_dir, "../b_M08b/03_PlotData/M08b_view_summary.csv"), show_col_types = FALSE)$label
top1_by_padjust <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(G1 = sign, View = factor(view_label, levels = view_order))

positive_data <- subset(top1_by_padjust, G1 == "positive") %>% arrange(View)
negative_data <- subset(top1_by_padjust, G1 == "negative")

positive_data$Description <- factor(
  positive_data$Description,
  levels = rev(unique(positive_data$Description))
)
negative_data$Description <- factor(
  negative_data$Description,
  levels = rev(unique(negative_data$Description[order(negative_data$View)]))
)

color_midpoint <- median(negative_data$p.adjust, na.rm = TRUE)

make_panel <- function(data, title_text) {
  ggplot(data, aes(x = Description, y = View, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.7) +
    coord_flip() +
    labs(title = title_text, x = NULL, y = NULL) +
    facet_grid(. ~ G1, scales = "free", space = "free", switch = "y",
      labeller = labeller(G1 = c(positive = "Increase with pcw", negative = "Decrease with pcw"))) +
    scale_color_gradient2(
      low = "#058786",
      mid = "#7EB5B4",
      high = "#224767",
      midpoint = color_midpoint
    ) +
    scale_size(
      breaks = function(x) pretty(x, n = 4),
      labels = function(x) as.integer(x), range = c(1, 4)
    ) +
    theme_linedraw(base_size = 12, base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 9),
      axis.text.y = element_text(color = "black", size = 9),
      plot.title = element_text(size = 12, hjust = 0.5, margin = margin(b = 6)),
      strip.text = element_text(size = 12, color = "white"),
      strip.background = element_rect(fill = "black", color = "black", linewidth = 0.8),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.35),
      panel.grid.minor = element_line(color = "#ECECEC", linewidth = 0.25),
      axis.ticks = element_line(color = "black", linewidth = 0.7),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      plot.margin = margin(8, 10, 8, 8)
    )
}

positive_plot <- make_panel(
  positive_data,
  "Top1 Pathways in Top100 Positive Features\nin Each View"
)
negative_plot <- make_panel(
  negative_data,
  "Top1 Pathways in Top100 Negative Features\nin Each View"
)
combined_plot <- positive_plot / negative_plot + plot_layout(ncol = 1)

save_pair(positive_plot, "M08f_positive_GO_BP", 9, 5)
save_pair(negative_plot, "M08f_negative_GO_BP", 9, 5)
save_pair(combined_plot, "M08f_combined_GO_BP", 9, 10)

