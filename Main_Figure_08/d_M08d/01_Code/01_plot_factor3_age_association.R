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


factor3_data <- read_csv(file.path(input_dir, "M08d_Factor3_scores.csv"), show_col_types = FALSE)

factor3_result <- read_csv(file.path(input_dir, "M08d_Factor3_statistics.csv"), show_col_types = FALSE)
age_boundary <- 32
factor_boundary <- predict(lm(value ~ pcw, data = factor3_data), data.frame(pcw = age_boundary))
F3scatter_plot <- factor3_data %>%
  ggplot(aes(x = pcw, y = value)) +
  geom_vline(xintercept = age_boundary, color = "#32348B", linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = factor_boundary, color = "#32348B", linetype = "dashed", linewidth = 0.35) +
  geom_point(color = "#272E6A", shape = 21, size = 2, stroke = 1.2) +
  geom_smooth(method = "lm", se = TRUE, color = "#272E6A", fill = "lightblue", alpha = 0.3) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf(
      "Spearman R = %.2f\np = %.5f\nq = %.5f",
      factor3_result$rho, factor3_result$p, factor3_result$q
    ),
    hjust = -0.1, vjust = 1, size = 3, color = "black", lineheight = 1.2, family = "Arial"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10), legend.position = "none", aspect.ratio = 1
  ) +
  ylab("Factor 3 value")
out8d <- file.path(.scf_project_root, "Main_Figure_08", "d_M08d", "02_Figures")
save_pair(F3scatter_plot, out8d, "M08d_Factor3_vs_pcw", 3, 3)


