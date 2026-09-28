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


f7b <- read_csv(file.path(input_dir, "M07b_Selected_Plasma_Protein_NPX.csv"), show_col_types = FALSE)
plasma_stats <- read_tsv(file.path(input_dir, "M07a_Plasma_Protein_Age_Statistics.tsv"), show_col_types = FALSE)
stat_index <- match(f7b$protein, plasma_stats$feature_display)
stopifnot(!anyDuplicated(plasma_stats$feature_display), !anyNA(stat_index))
f7b$rho <- plasma_stats$rho[stat_index]
f7b$p_value <- plasma_stats$p_value[stat_index]
f7b$q_value_BH_76 <- plasma_stats$q_value_BH_family[stat_index]

f7b_stats <- f7b %>%
  distinct(protein, rho, p_value, q_value_BH_76) %>%
  arrange(desc(rho))
f7b$protein <- factor(f7b$protein, levels = f7b_stats$protein)
f7b_stats$protein <- factor(f7b_stats$protein, levels = f7b_stats$protein)
palette_table <- read.delim(file.path(.scf_project_root, "00_Settings/palette_registry.tsv"))
protein_palette <- palette_table[palette_table$scope == "M07b", ]
protein_colors <- setNames(protein_palette$hex, protein_palette$element)
stopifnot(all(as.character(f7b_stats$protein) %in% names(protein_colors)))

p7b <- ggplot(f7b, aes(x = pcw, y = NPX, color = protein)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  theme_classic(base_family = "Arial") +
  xlab("pcw") +
  facet_wrap(~ protein, scales = "free", labeller = label_value, ncol = 7) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    aspect.ratio = 0.8
  ) +
  ylab("NPX") +
  geom_text(
    data = f7b_stats,
    aes(
      x = Inf, y = Inf,
      label = paste0(
        "r = ", round(rho, 2),
        "\np ", ifelse(p_value < 0.001, "< 0.001", paste0("= ", round(p_value, 3))),
        "\nq ", ifelse(q_value_BH_76 < 0.001, "< 0.001", paste0("= ", round(q_value_BH_76, 3)))
      )
    ),
    hjust = 1, vjust = 1, color = "black", inherit.aes = FALSE
  ) +
  scale_color_manual(values = protein_colors)

out7b <- file.path(.scf_project_root, "Main_Figure_07", "b_M07b", "02_Figures")
save_plot_pair(p7b, out7b, "M07b_Plasma_Olink_selected17", 15, 6.2)
write_csv(f7b, file.path(.scf_shared, "Main_Figure_07/b_M07b/03_PlotData/M07b_Plasma_Olink_selected17_plot_data.csv"))
write_csv(f7b_stats, file.path(.scf_shared, "Main_Figure_07/b_M07b/03_PlotData/M07b_Plasma_Olink_selected17_statistics.csv"))


