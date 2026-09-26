.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "ggsci", "scales", "tidyverse"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)

clean_root <- .scf_project_root
input_dir <- file.path(.scf_shared, "Main_Figure_04_05/01_Chain_Pairing_Clonotypes")
source(file.path(.scf_code_dir, "01_define_repertoire_palettes.R"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
})

build_abundance <- function(data, group_column) {
  counts <- table(data[[group_column]], data$chain_pairing)
  as.data.frame.matrix(counts)
}

plot_tcr <- function(data) {
  m04a_dir <- file.path(.scf_project_root, "Main_Figure_04", "a_M04a", "02_Figures")
  m04b_dir <- file.path(.scf_project_root, "Main_Figure_04", "b_M04b", "02_Figures")
  dir.create(m04a_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(m04b_dir, recursive = TRUE, showWarnings = FALSE)

  abundance_MainID <- build_abundance(data, "MainID")
  data_prop <- abundance_MainID %>%
    rownames_to_column(var = "tissue") %>%
    mutate(total = rowSums(.[, -1])) %>%
    mutate(across(
      starts_with("single") | starts_with("orphan") | starts_with("extra") | starts_with("two"),
      ~ . / total,
      .names = "prop_{.col}"
    ))
  data_long <- data_prop %>%
    select(tissue, starts_with("prop")) %>%
    pivot_longer(
      cols = starts_with("prop"),
      names_to = "category",
      values_to = "proportion"
    ) %>%
    mutate(category = str_remove(category, "prop_"))
  data_long$Organ <- str_sub(data_long$tissue, 1, 1)
  data_long <- data_long %>%
    mutate(Organ = case_when(
      Organ == "B" ~ "PBMC",
      Organ == "L" ~ "Liver",
      Organ == "S" ~ "Spleen",
      Organ == "T" ~ "Thymus"
    ))
  data_long$Organ <- factor(data_long$Organ, levels = c("PBMC", "Liver", "Thymus", "Spleen"))


  Ppcen_sample <- ggplot(data_long, aes(x = tissue, y = proportion, fill = category)) +
    geom_col() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Tissue", y = "Proportion", fill = "Category") +
    theme_classic() +
    xlab("") +
    theme(
      text = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 16),
      axis.text.x = element_text(angle = 90)
    ) +
    scale_fill_d3() +
    scale_fill_manual(values = col1) +
    facet_grid(. ~ Organ, space = "free", scales = "free", switch = "y")

  ggsave(
    plot = Ppcen_sample,
    filename = file.path(m04a_dir, "M04a_TCR_pairing_by_sample.pdf"),
    width = 12,
    height = 5
  )
  ggsave(
    plot = Ppcen_sample,
    filename = file.path(m04a_dir, "M04a_TCR_pairing_by_sample.png"),
    width = 12,
    height = 5,
    dpi = 300
  )
  write.csv(data_long, file.path(.scf_shared, "Main_Figure_04/a_M04a/03_PlotData/M04a_TCR_pairing_by_sample_data.csv"), row.names = FALSE)

  testabundance <- build_abundance(data, "Main_Organ")
  data_prop <- testabundance %>%
    rownames_to_column(var = "tissue") %>%
    mutate(total = rowSums(.[, -1])) %>%
    mutate(across(
      starts_with("single") | starts_with("orphan") | starts_with("extra") | starts_with("two"),
      ~ . / total,
      .names = "prop_{.col}"
    ))
  organ_long <- data_prop %>%
    select(tissue, starts_with("prop")) %>%
    pivot_longer(
      cols = starts_with("prop"),
      names_to = "category",
      values_to = "proportion"
    ) %>%
    mutate(category = str_remove(category, "prop_"))
  organ_long$tissue <- factor(organ_long$tissue, levels = c("PBMC", "Liver", "Thymus", "Spleen"))


  Ppcen <- ggplot(organ_long, aes(x = tissue, y = proportion, fill = category)) +
    geom_col() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Tissue", y = "Proportion", fill = "Category") +
    theme_classic() +
    xlab("") +
    theme(
      text = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 16),
      axis.text.x = element_text(angle = 90),
      aspect.ratio = 2
    ) +
    scale_fill_d3() +
    scale_fill_manual(values = col1)

  ggsave(
    plot = Ppcen,
    filename = file.path(m04b_dir, "M04b_TCR_pairing_by_organ.pdf"),
    width = 6,
    height = 5
  )
  ggsave(
    plot = Ppcen,
    filename = file.path(m04b_dir, "M04b_TCR_pairing_by_organ.png"),
    width = 6,
    height = 5,
    dpi = 300
  )
  write.csv(organ_long, file.path(.scf_shared, "Main_Figure_04/b_M04b/03_PlotData/M04b_TCR_pairing_by_organ_data.csv"), row.names = FALSE)
}

plot_bcr <- function(data) {
  m05a_dir <- file.path(.scf_project_root, "Main_Figure_05", "a_M05a", "02_Figures")
  m05b_dir <- file.path(.scf_project_root, "Main_Figure_05", "b_M05b", "02_Figures")
  dir.create(m05a_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(m05b_dir, recursive = TRUE, showWarnings = FALSE)

  abundance_MainID <- build_abundance(data, "MainID")
  data_prop <- abundance_MainID %>%
    rownames_to_column(var = "tissue") %>%
    mutate(total = rowSums(.[, -1])) %>%
    mutate(across(
      starts_with("single") | starts_with("orphan") | starts_with("extra") |
        starts_with("ambiguous") | starts_with("two"),
      ~ . / total,
      .names = "prop_{.col}"
    ))
  data_long <- data_prop %>%
    select(tissue, starts_with("prop")) %>%
    pivot_longer(
      cols = starts_with("prop"),
      names_to = "category",
      values_to = "proportion"
    ) %>%
    mutate(category = str_remove(category, "prop_"))
  data_long$Organ <- str_sub(data_long$tissue, 1, 1)
  data_long <- data_long %>%
    mutate(Organ = case_when(
      Organ == "B" ~ "PBMC",
      Organ == "L" ~ "Liver",
      Organ == "S" ~ "Spleen",
      Organ == "T" ~ "Thymus"
    ))
  data_long$Organ <- factor(data_long$Organ, levels = c("PBMC", "Liver", "Thymus", "Spleen"))
  data_long$category <- factor(
    data_long$category,
    levels = c("ambiguous", "extra VJ", "extra VDJ", "orphan VDJ", "orphan VJ", "single pair", "two full chains")
  )


  Ppcen_sample <- ggplot(data_long, aes(x = tissue, y = proportion, fill = category)) +
    geom_col() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Tissue", y = "Proportion", fill = "Category") +
    theme_classic() +
    xlab("") +
    theme(
      text = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 16),
      axis.text.x = element_text(angle = 90)
    ) +
    scale_fill_d3() +
    scale_fill_manual(values = col11) +
    facet_grid(. ~ Organ, space = "free", scales = "free", switch = "y")

  ggsave(
    plot = Ppcen_sample,
    filename = file.path(m05a_dir, "M05a_BCR_pairing_by_sample.pdf"),
    width = 12,
    height = 5
  )
  ggsave(
    plot = Ppcen_sample,
    filename = file.path(m05a_dir, "M05a_BCR_pairing_by_sample.png"),
    width = 12,
    height = 5,
    dpi = 300
  )
  write.csv(data_long, file.path(.scf_shared, "Main_Figure_05/a_M05a/03_PlotData/M05a_BCR_pairing_by_sample_data.csv"), row.names = FALSE)

  testabundance <- build_abundance(data, "Main_Organ")
  data_prop <- testabundance %>%
    rownames_to_column(var = "tissue") %>%
    mutate(total = rowSums(.[, -1])) %>%
    mutate(across(
      starts_with("single") | starts_with("orphan") | starts_with("extra") |
        starts_with("ambiguous") | starts_with("two"),
      ~ . / total,
      .names = "prop_{.col}"
    ))
  organ_long <- data_prop %>%
    select(tissue, starts_with("prop")) %>%
    pivot_longer(
      cols = starts_with("prop"),
      names_to = "category",
      values_to = "proportion"
    ) %>%
    mutate(category = str_remove(category, "prop_"))
  organ_long$tissue <- factor(organ_long$tissue, levels = c("PBMC", "Liver", "Thymus", "Spleen"))
  organ_long$category <- factor(
    organ_long$category,
    levels = c("ambiguous", "extra VJ", "extra VDJ", "orphan VDJ", "orphan VJ", "single pair", "two full chains")
  )


  Ppcen <- ggplot(organ_long, aes(x = tissue, y = proportion, fill = category)) +
    geom_col() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Tissue", y = "Proportion", fill = "Category") +
    theme_classic() +
    xlab("") +
    theme(
      text = element_text(size = 16),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 16),
      axis.text.x = element_text(angle = 90),
      aspect.ratio = 2
    ) +
    scale_fill_d3() +
    scale_fill_manual(values = col11)

  ggsave(
    plot = Ppcen,
    filename = file.path(m05b_dir, "M05b_BCR_pairing_by_organ.pdf"),
    width = 6,
    height = 5
  )
  ggsave(
    plot = Ppcen,
    filename = file.path(m05b_dir, "M05b_BCR_pairing_by_organ.png"),
    width = 6,
    height = 5,
    dpi = 300
  )
  write.csv(organ_long, file.path(.scf_shared, "Main_Figure_05/b_M05b/03_PlotData/M05b_BCR_pairing_by_organ_data.csv"), row.names = FALSE)
}

tcr <- read.csv(file.path(.scf_shared, "Main_Figure_04_05/01_Chain_Pairing_Clonotypes/M04c_TCR_Chain_Pairing_UMAP.csv"))
bcr <- read.csv(file.path(.scf_shared, "Main_Figure_04_05/01_Chain_Pairing_Clonotypes/M05c_BCR_Chain_Pairing_UMAP.csv"))


plot_tcr(tcr)
plot_bcr(bcr)
message("M04a/b and M05a/b complete")
