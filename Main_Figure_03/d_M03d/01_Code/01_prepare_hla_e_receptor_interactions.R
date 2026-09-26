.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "tidyr"))

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

root <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
plot_ready <- file.path(.scf_project_root, "Main_Figure_03/00_CellChat_Analysis/02_Data/03_Plot_Ready")
output_data <- file.path(root, "03_PlotData")
dir.create(output_data, recursive = TRUE, showWarnings = FALSE)

signals <- c(
  "HLA-E_CD94:NKG2E",
  "HLA-E_CD94:NKG2C",
  "HLA-E_CD94:NKG2A",
  "HLA-E_KLRK1",
  "HLA-E_KLRC2",
  "HLA-E_KLRC1",
  "CLEC1B_KLRB1",
  "CLEC2C_KLRB1",
  "CLEC2B_KLRB1",
  "CLEC2D_KLRB1",
  "CD99_CD99"
)
targets <- c("CX3CR1+ NK", "CXCR6+ NK", "CD56highCD16low NK")
excluded_sources <- c("Megakaryocytes", "Early_ERY", "Mid_ERY", "Late_ERY")
order_table <- read.delim(file.path(plot_ready, "M03bcd_exact_28class_celltype_order.tsv"), check.names = FALSE)
source_table <- order_table %>%
  filter(!model_cell_type %in% excluded_sources) %>%
  transmute(source_order = order, source = model_cell_type, source_display = gsub("Naïve", "Naive", display_cell_type, fixed = TRUE))

read_stage <- function(stage) {
  read.csv(file.path(plot_ready, paste0("PBMC_", stage, "_exact28_communications.csv")), check.names = FALSE) %>%
    mutate(stage = stage)
}

selected <- bind_rows(read_stage("Early"), read_stage("Late")) %>%
  filter(interaction_name %in% signals, target %in% targets, source %in% source_table$source) %>%
  left_join(source_table, by = "source") %>%
  mutate(
    dataset = factor(stage, levels = c("Early", "Late")),
    target = factor(target, levels = targets),
    interaction_name = factor(interaction_name, levels = rev(signals)),
    source = factor(source_display, levels = source_table$source_display),
    prob.original = as.numeric(prob),
    prob = ifelse(prob.original > 0 & prob.original < 1, -1 / log(prob.original), NA_real_),
    p_value = case_when(
      pval < 0.01 ~ "p < 0.01",
      TRUE ~ "0.01 <= p < 0.05"
    ),
    p_value = factor(p_value, levels = c("p < 0.01", "0.01 <= p < 0.05")),
    selection_basis = "Full significant CellChat network; fixed reference mechanism list"
  ) %>%
  arrange(dataset, target, interaction_name, source_order) %>%
  select(source, target, ligand, receptor, prob, pval, interaction_name, interaction_name_2, pathway_name, annotation, evidence, dataset, p_value, prob.original, source_order, selection_basis)

if (!nrow(selected)) stop("No selected ligand-receptor interactions found")
missing_signals <- setdiff(signals, unique(as.character(selected$interaction_name)))
if (length(missing_signals)) stop(paste("missing requested signals:", paste(missing_signals, collapse = ", ")))
write.csv(selected, file.path(.scf_shared, "Main_Figure_03/d_M03d/03_PlotData/M03d_diffcell_pathwayNK3lateEarly_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

