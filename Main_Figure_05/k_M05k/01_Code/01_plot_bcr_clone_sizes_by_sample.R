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
  "Main_Figure_04_05/01_Chain_Pairing_Clonotypes/M05k_BCR_clone_size_by_sample.csv"
)
output_dir <- file.path(.scf_code_dir, "../02_Figures")
plot_data_dir <- file.path(.scf_code_dir, "../03_PlotData")
dir.create(plot_data_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(input_path, check.names = FALSE)

meta$Main_Organ <- factor(meta$Main_Organ, levels = c("PBMC", "Liver", "Thymus", "Spleen"))
meta$clone_size <- factor(meta$clone_size, levels = c("1", "2", "3", "4", ">=6"))

colorname2 <- c(
  "#006D2C", "#B5AD64", "#9DA8E2", "#91C392", "#FF9900", "#46A040",
  "#FFC179", "#98D9E9", "#F6313E", "#FFA300", "#333366", "#FF5A00",
  "#663366", "#FF6666", "#8F1336", "#0081C9", "#001588", "#CC0033",
  "#CC9966", "#CC0033", "#999933", "#009966", "#CCCC33", "#CCFF99",
  "#333399", "#993333"
)


P10 <- ggplot(meta, aes(fill = clone_size, x = MainID, y = proportion)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Adjusted ID", y = "Count", fill = "clone_id_size_max_6") +
  facet_grid(. ~ Main_Organ, space = "free", scales = "free", switch = "y") +
  scale_fill_manual(values = colorname2) +
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
  labs(color = "Cell_type") +
  scale_fill_jama()

ggsave(
  plot = P10,
  filename = file.path(output_dir, "M05k_BCR_clone_size_by_sample.pdf"),
  width = 12,
  height = 5
)
ggsave(
  plot = P10,
  filename = file.path(output_dir, "M05k_BCR_clone_size_by_sample.png"),
  width = 12,
  height = 5,
  dpi = 300
)
write.csv(meta, file.path(.scf_shared, "Main_Figure_05/k_M05k/03_PlotData/M05k_BCR_clone_size_by_sample_data.csv"), row.names = FALSE)
message("M05k complete")
