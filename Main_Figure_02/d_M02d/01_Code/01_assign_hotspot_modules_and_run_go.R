.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("clusterProfiler", "DOSE", "dplyr", "ggplot2", "org.Hs.eg.db", "RColorBrewer"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_shared <- file.path(.scf_project_root, "00_Global_Data", "Shared_Inputs")
.scf_panel_root <- normalizePath(file.path(.scf_code_dir, ".."), winslash = "/", mustWork = FALSE)

options(stringsAsFactors = FALSE)

output_dir <- file.path(.scf_panel_root, "02_Data/02_Processed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

early_fresh_path <- file.path(.scf_shared, "Main_Figure_02/d_M02d/M02d_Early_Detected_Modules.csv")
late_fresh_path <- file.path(.scf_shared, "Main_Figure_02/d_M02d/M02d_Late_Detected_Modules.csv")
early_reference_path <- file.path(.scf_shared, "Main_Figure_02/d_M02d/M02d_Early_Reference_Modules.csv")
late_reference_path <- file.path(.scf_shared, "Main_Figure_02/d_M02d/M02d_Late_Reference_Modules.csv")

read_fresh_modules <- function(path) {
  data <- read.csv(path, check.names = FALSE)
  colnames(data)[1] <- "Gene"
  data$Gene <- as.character(data$Gene)
  data
}

read_reference_modules <- function(path) {
  data <- read.csv(path, row.names = 1, check.names = FALSE)
  data$Gene <- rownames(data)
  data
}

map_modules_by_gene_overlap <- function(fresh, reference, stage) {
  mapping_rows <- lapply(sort(unique(fresh$Module)), function(module_id) {
    genes <- fresh$Gene[fresh$Module == module_id]
    overlap <- reference[reference$Gene %in% genes, , drop = FALSE]
    counts <- sort(table(overlap$Module), decreasing = TRUE)
    if (length(counts) == 0) {
      stop("No overlap found for ", stage, " module ", module_id)
    }
    top_count <- max(counts)
    top_labels <- sort(names(counts)[counts == top_count])
    data.frame(
      Stage = stage,
      Fresh_Module = module_id,
      Canonical_Module = top_labels[[1]],
      Fresh_Genes = length(unique(genes)),
      Overlap_Genes = sum(counts),
      Winning_Overlap = top_count,
      Winning_Fraction = top_count / length(unique(genes)),
      Tie_Count = length(top_labels),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, mapping_rows)
}

early_fresh <- read_fresh_modules(early_fresh_path)
late_fresh <- read_fresh_modules(late_fresh_path)
early_reference <- read_reference_modules(early_reference_path)
late_reference <- read_reference_modules(late_reference_path)

early_mapping <- map_modules_by_gene_overlap(early_fresh, early_reference, "Early")
late_mapping <- map_modules_by_gene_overlap(late_fresh, late_reference, "Late")
module_mapping <- rbind(early_mapping, late_mapping)
write.csv(
  module_mapping,
  file.path(output_dir, "M02d_numeric_to_original_module_mapping.csv"),
  row.names = FALSE
)

early_fresh$Module <- early_mapping$Canonical_Module[
  match(early_fresh$Module, early_mapping$Fresh_Module)
]
late_fresh$Module <- late_mapping$Canonical_Module[
  match(late_fresh$Module, late_mapping$Fresh_Module)
]


query_genes <- function(genes, bk_grd, p_thresh = 0.05) {
  ego <- enrichGO(
    genes,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    keyType = "SYMBOL"
  )
  return(ego)
}

run_modules <- function(module_data) {
  results_by_module <- list()
  for (mod_itr in sort(unique(module_data$Module))) {
    message("enrichGO: ", mod_itr)
    range_genes <- module_data[module_data$Module == mod_itr, ]$Gene
    ego <- query_genes(range_genes)
    results_by_module[[paste("Module", mod_itr, sep = "-")]] <- ego
  }
  results_by_module
}

Earlyresults_by_module <- run_modules(early_fresh)
Lateresults_by_module <- run_modules(late_fresh)
saveRDS(
  Earlyresults_by_module,
  file.path(output_dir, "M02d_Early_GO_results_by_original_module.rds")
)
saveRDS(
  Lateresults_by_module,
  file.path(output_dir, "M02d_Late_GO_results_by_original_module.rds")
)


