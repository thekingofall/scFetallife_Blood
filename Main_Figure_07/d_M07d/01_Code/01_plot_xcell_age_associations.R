.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "grid", "ragg", "readr", "scales"))

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
  library(readr)
  library(scales)
})

repo <- .scf_project_root
input_dir <- file.path(.scf_shared, "Main_Figure_07")
panel_root <- normalizePath(file.path(.scf_start_dir, ".."), winslash="/")
plot_data <- file.path(panel_root,"03_PlotData")
dir.create(plot_data,recursive=TRUE,showWarnings=FALSE)

save_plot_pair <- function(plot, directory, stem, width, height, dpi = 320) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".png")), plot, width = width, height = height, device = ragg::agg_png, dpi = dpi, bg = "white")
}


colorname2 <- c(
  "#46A040", "#00AF99", "#FFC179", "#98D9E9", "#F6313E", "#FFA300", "#333366", "#FF5A00", "#663366", "#FF6666",
  "#8F1336", "#0081C9", "#001588", "#CC0033", "#CC9966", "#CC0033", "#999933", "#009966", "#CCCC33", "#333399", "#993333",
  "#490C65", "#BA7FD0", "#A6CEE3", "#1F78B4", "#DE77AE", "#006D2C", "#868686", "#B5AD64", "#9DA8E2", "#91C392", "#FF9900", "#339966"
)

f7d <- read_csv(file.path(input_dir, "M07d_xCell_Scores.csv"), show_col_types = FALSE)

display_order <- c("Naïve CD4+ T cell", "Th2", "CD8+ T cell", "Central memory CD8+ T cell", "ImmuneScore")
f7d$display_name <- factor(f7d$display_name, levels = display_order)
f7d_stats <- f7d %>%
  distinct(display_name, notebook_r, notebook_p, notebook_q_BH_39)
score_colors <- setNames(colorname2[c(18, 30, 31, 12, 7)], display_order)

p7d <- ggplot(f7d, aes(x = actual_pcw, y = fraction, color = display_name)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  theme_classic() +
  xlab("pcw") +
  facet_wrap(~ display_name, scales = "free", labeller = label_value, nrow = 1) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8),
    aspect.ratio = 0.9
  ) +
  ylab("Proportion") +
  geom_text(
    data = f7d_stats,
    aes(
      x = -Inf, y = Inf,
      label = paste0(
        "r = ", round(notebook_r, 2),
        "\np ", ifelse(notebook_p < 0.001, "< 0.001", paste0("= ", round(notebook_p, 3))),
        "\nq ", ifelse(notebook_q_BH_39 < 0.001, "< 0.001", paste0("= ", round(notebook_q_BH_39, 3)))
      )
    ),
    hjust = 0, vjust = 1, color = "black", inherit.aes = FALSE
  ) +
  scale_color_manual(values = score_colors)

out7d <- file.path(.scf_project_root, "Main_Figure_07", "d_M07d", "02_Figures")
save_plot_pair(p7d, out7d, "M07d_Stimulated_PBMC_xCell5", 13, 4.0)
write_csv(f7d, file.path(.scf_shared, "Main_Figure_07/d_M07d/03_PlotData/M07d_Stimulated_PBMC_xCell5_plot_data.csv"))
write_csv(f7d_stats, file.path(.scf_shared, "Main_Figure_07/d_M07d/03_PlotData/M07d_Stimulated_PBMC_xCell5_statistics.csv"))

