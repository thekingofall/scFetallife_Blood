.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "ggsci", "tidyverse"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(tidyverse)
})

clean_root <- .scf_project_root
input_path <- file.path(.scf_shared, "Main_Figure_04_05/01_Chain_Pairing_Clonotypes/M05j_BCR_Isotypes_by_Sample.csv")
output_dir <- file.path(.scf_code_dir, "../02_Figures")
plot_data_dir <- file.path(.scf_code_dir, "../03_PlotData")
dir.create(plot_data_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(input_path, check.names = FALSE)

meta <- subset(meta, isotype_status != "Multi") %>%
  subset(isotype_status != "None") %>%
  subset(isotype_status != "IgL")
meta$Main_Organ <- factor(meta$Main_Organ, levels = c("PBMC", "Liver", "Thymus", "Spleen"))


Px1 <- ggplot(meta, aes(fill = isotype_status, x = MainID, y = proportion)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Adjusted ID", y = "Count", fill = "Isotype Status") +
  facet_grid(. ~ Main_Organ, space = "free", scales = "free", switch = "y") +
  scale_fill_nejm() +
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
  labs(color = "Cell_type")

ggsave(
  plot = Px1,
  filename = file.path(output_dir, "M05j_BCR_isotype_by_sample.pdf"),
  width = 12,
  height = 5
)
ggsave(
  plot = Px1,
  filename = file.path(output_dir, "M05j_BCR_isotype_by_sample.png"),
  width = 12,
  height = 5,
  dpi = 300
)
message("M05j complete")
