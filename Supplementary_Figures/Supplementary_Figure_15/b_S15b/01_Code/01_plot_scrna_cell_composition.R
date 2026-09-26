.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "readr"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

repo <- .scf_project_root
input_path <- file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_15/b_S15b/S15b_scRNA_Cell_Composition.csv")
output_dir <- file.path(.scf_project_root, "Supplementary_Figures/Supplementary_Figure_15", "b_S15b", "02_Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

S1 <- read_csv(input_path, show_col_types = FALSE)


stopifnot(all(abs(S1 %>% group_by(MainID) %>% summarise(total = sum(percentage), .groups = "drop") %>% pull(total) - 100) < 1e-8))

S1$new_Body <- S1$Main_Organ
S1$Last_cell_type <- S1$display_cell_type
color_mapping <- c(
  "CXCR5-B" = "#FF0000",
  "CXCR5+B" = "#FF9B9B",
  "Immature NK" = "#13C0DF",
  "Mature NK" = "#077E97",
  "Innate T" = "#00C000",
  "Gamma Delta T" = "#0000FF",
  "Naïve CD8 T" = "#FF6000",
  "CD8 T" = "#FF6000",
  "Naïve CD4 T" = "#B856D7",
  "CD4 T" = "#B856D7",
  "abT(entry)" = "#6a73cf",
  "CD3+double positive T" = "#A2A200",
  "CD3+double negative T" = "#E8E800",
  "Others" = "#8C8ED2"
)
S1$new_Body <- factor(S1$new_Body, levels = c("PBMC", "Liver", "Thymus", "Spleen"))
S1$Last_cell_type <- factor(
  S1$Last_cell_type,
  levels = c(
    "CXCR5-B", "CXCR5+B", "Immature NK", "Mature NK", "CD8 T", "CD4 T",
    "Innate T", "CD3+double positive T", "CD3+double negative T", "abT(entry)", "Others"
  )
)

P6 <- ggplot(S1, aes(x = MainID, y = count)) +
  geom_bar(aes(fill = Last_cell_type), stat = "identity", position = "fill") +
  facet_grid(. ~ new_Body, space = "free", scales = "free", switch = "y") +
  scale_fill_manual(values = color_mapping) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  ) +
  ylab("Composition (%)") +
  scale_y_continuous(labels = seq(0, 100, by = 25)) +
  xlab("") +
  labs(fill = "Cell_type")

ggsave(file.path(output_dir, "S15b_scRNAseq_Cell_Composition.pdf"), P6, width = 12, height = 6, useDingbats = FALSE)
ggsave(file.path(output_dir, "S15b_scRNAseq_Cell_Composition.png"), P6, width = 12, height = 6, dpi = 300)
write_csv(S1, file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_15/b_S15b/03_PlotData/S15b_scRNAseq_Cell_Composition_source.csv"))

cat("Saved S15b cell composition.\n")
