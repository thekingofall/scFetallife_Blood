.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "RColorBrewer", "tidyr"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
})

input_path <- file.path(.scf_shared, "Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv")
output_root <- file.path(.scf_project_root, "Supplementary_Figures/Supplementary_Figure_01", "b_S01b", "02_Figures")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

adata1_obs_all <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
adata1_obs_all <- adata1_obs_all[, !duplicated(names(adata1_obs_all)), drop = FALSE]
adata1_obs <- adata1_obs_all


process_grouped_data <- function(adata1_obs) {
  grouped_data <- adata1_obs %>%
    group_by(MainID, Last_cell_type_num, Cell_lineage) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    group_by(MainID) %>%
    mutate(week_total = sum(count)) %>%
    ungroup()

  grouped_data <- grouped_data %>%
    separate(MainID, into = c("Part1", "Part2"), sep = "\\.", remove = FALSE) %>%
    separate(Part1, into = c("Body", "Post_Conception_Age_Weeks"), sep = "(?<=\\D)(?=\\d)", remove = FALSE) %>%
    mutate(Post_Conception_Age_Weeks = as.numeric(sub("^[A-Za-z]+([0-9]+([.][0-9]+)?).*", "\\1", MainID)), percentage = count / week_total, week_number = row_number())

  adjustedID_levels <- unique(grouped_data$MainID)
  adjustedID_mapping <- setNames(seq_along(adjustedID_levels), adjustedID_levels)
  grouped_data <- grouped_data %>%
    mutate(MainID_numeric = unname(adjustedID_mapping[MainID])) %>%
    mutate(MainID_letter = substr(MainID, 1, 1)) %>%
    group_by(MainID_letter) %>%
    mutate(MainID_numeric = match(MainID, unique(MainID))) %>%
    ungroup()
  grouped_data
}

grouped_data_result <- process_grouped_data(adata1_obs)
prefixes <- as.numeric(gsub("_.*$", "", levels(factor(adata1_obs$Last_cell_type_num))))
sorted_indices <- order(prefixes)
sorted_cell_types <- levels(factor(adata1_obs$Last_cell_type_num))[sorted_indices]
grouped_data_result$Last_cell_type_num <- factor(grouped_data_result$Last_cell_type_num, levels = sorted_cell_types)
grouped_data_result$Cell_lineage <- factor(
  grouped_data_result$Cell_lineage,
  levels = c("PRECURSOR", "B_CELL", "T/ILC", "NK", "MYELOID", "DC", "MK/ERY", "OTHERS")
)
grouped_data_result$Body <- factor(grouped_data_result$Body, levels = c("B", "L", "T", "S"))

resultw <- grouped_data_result %>%
  group_by(Body, Last_cell_type_num, Cell_lineage) %>%
  summarise(WeightedAverage = weighted.mean(count, week_total), .groups = "drop")
resultw$Body <- factor(resultw$Body, levels = c("B", "L", "T", "S"), labels = c("PBMC", "Liver", "Thymus", "Spleen"))

wegitcountplot <- ggplot(resultw, aes(x = Last_cell_type_num, y = Body, fill = WeightedAverage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(WeightedAverage)), size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ Cell_lineage, space = "free", scales = "free", switch = "y") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  xlab("") +
  ylab("Weight Count") +
  scale_fill_gradientn(colours = rev(colorRampPalette(brewer.pal(9, "Spectral"))(100)))

source_dir <- file.path(.scf_code_dir, "../03_Source_Data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_root, "S01b_weighted_cell_abundance.pdf"), plot = wegitcountplot, width = 20, height = 4, useDingbats = FALSE)
ggsave(file.path(output_root, "S01b_weighted_cell_abundance.png"), plot = wegitcountplot, width = 20, height = 4, dpi = 300)
write.csv(resultw, file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_01/b_S01b/03_PlotData/S01b_weighted_cell_abundance_source.csv"), row.names = FALSE)
write.csv(grouped_data_result, file.path(source_dir, "S01b_sample_cell_counts.csv"), row.names = FALSE)
