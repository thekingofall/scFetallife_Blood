.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("CellChat", "circlize", "ComplexHeatmap", "dplyr", "ggplot2", "grid", "jsonlite", "tidyr"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(CellChat); library(dplyr); library(tidyr); library(jsonlite)})
root <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
processed <- file.path(.scf_shared, "Main_Figure_03", "00_CellChat_Models")
plot_ready <- file.path(root, "02_Data/03_Plot_Ready")
b_data <- file.path(.scf_project_root, "Main_Figure_03/b_M03b/03_PlotData")
c_data <- file.path(.scf_project_root, "Main_Figure_03/c_M03c/03_PlotData")
for (path in c(b_data, c_data)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
order_table <- read.delim(file.path(plot_ready, "M03bcd_exact_28class_celltype_order.tsv"), check.names = FALSE)
display_order <- order_table$display_cell_type
model_order <- order_table$model_cell_type
names(display_order) <- model_order

color_json <- file.path(.scf_shared, "Main_Figure_03/a_M03a/M03_Cell_Type_Colors.json")
cell_colors <- unlist(fromJSON(color_json))
if (!all(model_order %in% names(cell_colors))) stop("cell-type color dictionary is incomplete")
cell_colors <- cell_colors[model_order]
names(cell_colors) <- unname(display_order[model_order])

early <- readRDS(file.path(processed, "PBMC_Early_exact28_cellchat.rds"))
late <- readRDS(file.path(processed, "PBMC_Late_exact28_cellchat.rds"))


expand_network <- function(object, measure) {
  result <- matrix(0, nrow = length(model_order), ncol = length(model_order), dimnames = list(model_order, model_order))
  source <- object@net[[measure]]
  common_row <- intersect(rownames(source), model_order)
  common_col <- intersect(colnames(source), model_order)
  result[common_row, common_col] <- source[common_row, common_col, drop = FALSE]
  result
}

count_early <- expand_network(early, "count")
count_late <- expand_network(late, "count")
weight_early <- expand_network(early, "weight")
weight_late <- expand_network(late, "weight")
count_diff <- count_early - count_late
weight_diff <- weight_early - weight_late
dimnames(count_diff) <- list(unname(display_order[rownames(count_diff)]), unname(display_order[colnames(count_diff)]))
dimnames(weight_diff) <- list(unname(display_order[rownames(weight_diff)]), unname(display_order[colnames(weight_diff)]))

network_long <- bind_rows(
  as.data.frame(as.table(count_diff), stringsAsFactors = FALSE) %>% rename(source = Var1, target = Var2, difference = Freq) %>% mutate(measure = "count"),
  as.data.frame(as.table(weight_diff), stringsAsFactors = FALSE) %>% rename(source = Var1, target = Var2, difference = Freq) %>% mutate(measure = "weight")
)
write.csv(network_long, file.path(.scf_shared, "Main_Figure_03/b_M03b/03_PlotData/M03b_differential_interaction_matrix.csv"), row.names = FALSE)

pbmc <- mergeCellChat(list(late, early), add.names = c("PBMC_Late", "PBMC_Early"))
joint_levels <- levels(pbmc@idents$joint)

signaling_targets <- c("CX3CR1+ NK", "CXCR6+ NK", "CD56highCD16low NK", "Naïve CD8 T", "NK T")
signaling_data <- lapply(signaling_targets, function(target) {
  result <- netAnalysis_signalingChanges_scatter(pbmc, idents.use = target, color.use = c("#7A4E87", "#4F8F8A", "#C85A45"))
  data <- result$data
  data$Group <- target
  data
})
names(signaling_data) <- signaling_targets
write.csv(bind_rows(signaling_data), file.path(plot_ready, "M03c_all_signaling_changes_exact28.csv"), row.names = FALSE)

selected_signaling <- bind_rows(
  signaling_data[["CX3CR1+ NK"]] %>% filter(incoming > 0.3),
  signaling_data[["CXCR6+ NK"]] %>% filter(incoming > 0.3),
  signaling_data[["CD56highCD16low NK"]] %>% filter(incoming > 0.3),
  signaling_data[["Naïve CD8 T"]] %>% filter(incoming > 0.3),
  signaling_data[["NK T"]] %>% filter(incoming < -0.3)
)
selected_signaling$Group <- factor(selected_signaling$Group, levels = signaling_targets)
write.csv(selected_signaling, file.path(.scf_shared, "Main_Figure_03/c_M03c/03_PlotData/M03c_signaling_changes_data.csv"), row.names = FALSE)

