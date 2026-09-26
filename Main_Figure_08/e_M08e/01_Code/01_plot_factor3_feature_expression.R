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


Weight_long <- read_csv(file.path(input_dir, "M08e_Factor3_top2_heatmap_values.csv"), show_col_types = FALSE)

Weight_long <- Weight_long %>%
  mutate(
    variable = factor(sample, levels = unique(sample[order(pcw)])),
    features2 = feature_label,
    View = factor(view_label, levels = unique(view_label)),
    np = factor(sign, levels = c("positive", "negative")),
    value = expression_value
  )

sample_ages <- unique(Weight_long[, c("sample", "pcw")])
age_boundary <- sum(sample_ages$pcw < 32) + 0.5
heatmap_original <- function(data, show_legend = TRUE) {
  ggplot(data, aes(x = variable, y = features2, fill = value)) +
    scale_fill_gradient2(low = "#334555", mid = "white", high = "#8E328A", midpoint = 0) +
    geom_tile() +
    geom_vline(xintercept = age_boundary, color = "#D93B35", linetype = "dashed", linewidth = 0.3) +
    facet_grid(View ~ np, scales = "free", space = "free", switch = "y",
      labeller = labeller(np = c(positive = "Increase with pcw", negative = "Decrease with pcw"))) +
    theme_minimal(base_size = 10, base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(hjust = 1, vjust = 1),
      strip.background = element_rect(color = "black", fill = "white", linetype = "solid"),
      strip.placement = "outside", strip.switch.pad.grid = unit(0.4, "cm"),
      strip.text.y.left = element_text(angle = 0),
      legend.position = if (show_legend) "right" else "none"
    ) +
    labs(x = "Sample ID", y = "Gene") +
    ylab("")
}

all_heatmap2up3 <- heatmap_original(subset(Weight_long, np == "positive"), FALSE) +
  theme(strip.text.y.left = element_blank(), strip.background.y = element_blank())
all_heatmap2down3 <- heatmap_original(subset(Weight_long, np == "negative"), TRUE)
combined_plot3 <- all_heatmap2up3 + all_heatmap2down3 +
  plot_layout(ncol = 2)
out8e <- file.path(.scf_project_root, "Main_Figure_08", "e_M08e", "02_Figures")
save_pair(combined_plot3, out8e, "M08e_Factor3_top2_heatmap", 10, 10)


