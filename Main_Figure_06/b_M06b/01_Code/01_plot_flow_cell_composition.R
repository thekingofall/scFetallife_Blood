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
input_dir <- file.path(.scf_shared, "Main_Figure_06/bc_M06bc")
panel_root <- normalizePath(file.path(.scf_start_dir, ".."), winslash="/")
plot_data <- file.path(panel_root,"03_PlotData")
dir.create(plot_data,recursive=TRUE,showWarnings=FALSE)

save_plot_pair <- function(plot, directory, stem, width, height, dpi = 320) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".png")), plot, width = width, height = height, device = ragg::agg_png, dpi = dpi, bg = "white")
}


f6b <- read_tsv(
  file.path(input_dir, "M06b_Flow_Cell_Proportions_by_Sample.tsv"),
  show_col_types = FALSE
)


cell_levels <- c(
  "CXCR5- B", "CXCR5+ B", "Immature NK", "Mature NK",
  "Naïve CD8+ T", "Central memory-like CD8+ T",
  "Effective memory-like CD8+ T", "Terminal differentiated CD8+ T",
  "Naïve CD4+ T", "Central memory-like CD4+ T",
  "Effective memory-like CD4+ T", "Terminal differentiated CD4+ T",
  "DNT", "LDP", "EDP", "CD8+ ISP", "CD4+ ISP", "Ungate"
)
f6b$Cell_type <- factor(f6b$Cell_type, levels = cell_levels)
f6b$Tissue <- factor(f6b$Tissue, levels = c("PBMC", "Liver", "Thymus", "Spleen"))
sample_levels <- f6b %>%
  distinct(Tissue, MainID) %>%
  arrange(Tissue, readr::parse_number(MainID)) %>%
  pull(MainID)
f6b$MainID <- factor(f6b$MainID, levels = unique(sample_levels))

palette18 <- c(
  "CXCR5- B" = "#DC3B35", "CXCR5+ B" = "#E78F96",
  "Immature NK" = "#43B3CB", "Mature NK" = "#12859A",
  "Naïve CD8+ T" = "#C66B39", "Central memory-like CD8+ T" = "#DF5B39",
  "Effective memory-like CD8+ T" = "#EB893E", "Terminal differentiated CD8+ T" = "#F1B574",
  "Naïve CD4+ T" = "#87499A", "Central memory-like CD4+ T" = "#96549B",
  "Effective memory-like CD4+ T" = "#B57DB5", "Terminal differentiated CD4+ T" = "#D8C3D8",
  "DNT" = "#EAE638", "LDP" = "#A2A044", "EDP" = "#DD88B1",
  "CD8+ ISP" = "#D32F74", "CD4+ ISP" = "#5EBD54", "Ungate" = "#8888C5"
)

f6b_mean <- f6b %>%
  group_by(Tissue, Cell_type) %>%
  summarise(Percentage = mean(Percentage), .groups = "drop") %>%
  mutate(MainID = factor(as.character(Tissue), levels = levels(f6b$Tissue)))

p6b_samples <- ggplot(f6b, aes(x = MainID, y = Percentage, fill = Cell_type)) +
  geom_col(width = 0.92) +
  facet_grid(. ~ Tissue, space = "free", scales = "free", switch = "y") +
  scale_fill_manual(values = palette18, drop = FALSE) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.35, "cm")
  ) +
  ylab("Composition (%)") +
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = c(0, 0)) +
  xlab("Sample ID") +
  labs(fill = NULL) +
  guides(fill = guide_legend(ncol = 5, byrow = TRUE))

p6b_mean <- ggplot(f6b_mean, aes(x = MainID, y = Percentage, fill = Cell_type)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = palette18, drop = FALSE) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  ylab("Composition (%)") +
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = c(0, 0)) +
  xlab("")

out6b <- file.path(.scf_project_root, "Main_Figure_06", "b_M06b", "02_Figures")
save_plot_pair(p6b_samples, out6b, "M06b_SFC_18type_per_sample", 14, 7.2)
save_plot_pair(p6b_mean, out6b, "M06b_SFC_18type_tissue_mean", 5.2, 6.2)
write_csv(f6b, file.path(.scf_shared, "Main_Figure_06/b_M06b/03_PlotData/M06b_SFC_18type_per_sample_plot_data.csv"))
write_csv(f6b_mean, file.path(.scf_shared, "Main_Figure_06/b_M06b/03_PlotData/M06b_SFC_18type_tissue_mean_plot_data.csv"))


