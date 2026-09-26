.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "ggsci"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
})

clean_root <- .scf_project_root
input_path <- file.path(
  .scf_shared,
  "Main_Figure_04_05/01_Chain_Pairing_Clonotypes/M04i_TCR_clone_size_by_sample.csv"
)
output_dir <- file.path(.scf_code_dir, "../02_Figures")
plot_data_dir <- file.path(.scf_code_dir, "../03_PlotData")
dir.create(plot_data_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(input_path, check.names = FALSE)

stopifnot(all(c("MainID", "Main_Organ", "clone_size", "proportion") %in% names(meta)))

organ_levels <- c("PBMC", "Liver", "Thymus", "Spleen")
clone_levels <- c("1", "2", "3", "4", "5", "6", ">=7")
clone_labels <- c("1", "2", "3", "4", "5", "6", "\u22657")

meta$MainID <- factor(meta$MainID, levels = sort(unique(meta$MainID)))
meta$Main_Organ <- factor(meta$Main_Organ, levels = organ_levels)
meta$clone_size_plot <- factor(
  meta$clone_size,
  levels = clone_levels,
  labels = clone_labels
)


P7 <- ggplot(
  meta,
  aes(x = MainID, y = proportion, fill = clone_size_plot)
) +
  geom_col(width = 0.9) +
  labs(x = "", y = "Proportion", fill = "Clonotype size") +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_grid(. ~ Main_Organ, space = "free", scales = "free", switch = "y") +
  scale_fill_lancet(drop = FALSE)

ggsave(
  plot = P7,
  filename = file.path(output_dir, "M04i_TCR_clone_size_by_sample.pdf"),
  width = 12,
  height = 4,
  device = cairo_pdf
)
ggsave(
  plot = P7,
  filename = file.path(output_dir, "M04i_TCR_clone_size_by_sample.png"),
  width = 12,
  height = 4,
  dpi = 300
)

meta$MainID <- as.character(meta$MainID)
meta$Main_Organ <- as.character(meta$Main_Organ)
meta$clone_size_plot <- NULL
write.csv(
  meta,
  file.path(.scf_shared, "Main_Figure_04/i_M04i/03_PlotData/M04i_TCR_clone_size_by_sample_data.csv"),
  row.names = FALSE
)
message("M04i complete")
