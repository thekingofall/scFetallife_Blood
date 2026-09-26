.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("CellChat", "future", "Matrix", "SeuratObject"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(CellChat)
  library(SeuratObject)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(file.path(.scf_code_dir, ".."), winslash = "/", mustWork = TRUE)
processed <- file.path(.scf_shared, "Main_Figure_03", "00_CellChat_Models")
plot_ready <- file.path(root, "02_Data", "03_Plot_Ready")
dir.create(processed, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_ready, recursive = TRUE, showWarnings = FALSE)

source_rds <- Sys.getenv("SCFETALLIFE_SEURAT_RDS", unset = file.path(.scf_global, "Fetal_Immune_Atlas_Seurat.rds"))

display_order <- c(
  "HSC_MPP", "MEP", "MEMP", "Pro-B", "Large pre-B", "Small pre-B",
  "CXCR5- Naïve B", "CXCR5+ Naïve B", "DN(Q) T", "DP(P) T", "DP(Q) T",
  "Treg", "Cycling Treg", "Naïve CD4 T", "Naïve CD8 T", "abT(entry)",
  "Tem", "Th17like_INNATE_T", "NK T", "Gamma Delta V1 T", "Gamma Delta V2 T",
  "GNG4 +CD8aa+T", "ILC2/3", "CX3CR1+ NK", "CXCR6+ NK",
  "CD56highCD16low NK", "Myeloid-CD177", "Classical Monocytes",
  "CD14+PPBP+ Monocytes", "Macrophages", "DC1", "DC2", "pDC",
  "Megakaryocytes", "Early_ERY", "Mid_ERY", "Late_ERY", "Endothelial cells", "Others"
)

contact_db <- function() {
  db <- subsetDB(CellChatDB.human, search = "Cell-Cell Contact", key = "annotation")

  db
}

run_one <- function(data, labels, group_name, db) {
  if (!identical(colnames(data), names(labels))) stop("expression and label order mismatch")
  set.seed(1L)
  meta <- data.frame(labels = as.character(labels), row.names = names(labels), stringsAsFactors = FALSE)
  cellchat <- createCellChat(object = data)
  cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels")
  cellchat@DB <- db
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    trim = 0.1,
    raw.use = TRUE,
    population.size = FALSE,
    nboot = 100L,
    seed.use = 1L
  )
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat@options$analysis_group <- group_name
  cellchat@options$display_class_count <- 28L
  cellchat@options$filterCommunication_applied <- FALSE
  cellchat
}

if (as.character(packageVersion("CellChat")) != "1.1.2") stop("CellChat version must be 1.1.2")
options(future.globals.maxSize = 64 * 1024^3)
future::plan("sequential")

message("Loading normalized Seurat input")
master <- readRDS(source_rds)

md <- master@meta.data
md$cell_id <- rownames(md)
md$Category <- ifelse(grepl("^B", as.character(md$MainID)), "PBMC", ifelse(grepl("^L", as.character(md$MainID)), "Liver", ifelse(grepl("^T", as.character(md$MainID)), "Thymus", ifelse(grepl("^S", as.character(md$MainID)), "Spleen", "Other"))))
md$NumericPart <- as.numeric(sub(".*?([0-9]+\\.[0-9]+).*", "\\1", as.character(md$MainID)))
md$Stage <- ifelse(md$NumericPart <= 26, "Early", "Late")
md$GroupStat <- paste0(md$Category, "_", md$Stage)

pbmc_original <- md[md$Category == "PBMC", , drop = FALSE]
original_counts <- table(as.character(pbmc_original$Last_cell_type))
selected_display <- display_order[display_order %in% names(original_counts)[original_counts > 100L] & display_order != "Others"]

selected_model <- gsub("[()]", "", selected_display)
pbmc_model_label <- gsub("[()]", "", as.character(md$Last_cell_type))
selected <- md$Category == "PBMC" & pbmc_model_label %in% selected_model




counts <- do.call(rbind, lapply(c("PBMC_Early", "PBMC_Late"), function(group_name) {
  data.frame(
    group = group_name,
    retained_cells_28class = sum(selected & md$GroupStat == group_name),
    stringsAsFactors = FALSE
  )
}))
write.table(counts, file.path(plot_ready, "M03bcd_exact_28class_input_census.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(order = seq_along(selected_display), display_cell_type = selected_display, model_cell_type = selected_model), file.path(plot_ready, "M03bcd_exact_28class_celltype_order.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

db <- contact_db()
expression <- master@assays$RNA@data
model_summary <- list()
for (group_name in c("PBMC_Early", "PBMC_Late")) {
  cells <- md$cell_id[selected & md$GroupStat == group_name]
  labels <- pbmc_model_label[match(cells, md$cell_id)]
  names(labels) <- cells
  if (any(!labels %in% selected_model)) stop("unexpected CellChat class")
  object_path <- file.path(processed, paste0(group_name, "_exact28_cellchat.rds"))
  message("Running ", group_name, " with ", length(cells), " cells and ", length(unique(labels)), " observed classes")
  object <- run_one(expression[, cells, drop = FALSE], labels, group_name, db)
  saveRDS(object, object_path, compress = TRUE)
  comm <- subsetCommunication(object)
  comm$dataset <- group_name
  write.csv(comm, file.path(plot_ready, paste0(group_name, "_exact28_communications.csv")), row.names = FALSE)
  count_table <- as.data.frame(as.table(object@net$count), stringsAsFactors = FALSE)
  colnames(count_table) <- c("source", "target", "count")
  count_table$dataset <- group_name
  write.csv(count_table, file.path(plot_ready, paste0(group_name, "_exact28_network_count.csv")), row.names = FALSE)
  weight_table <- as.data.frame(as.table(object@net$weight), stringsAsFactors = FALSE)
  colnames(weight_table) <- c("source", "target", "weight")
  weight_table$dataset <- group_name
  write.csv(weight_table, file.path(plot_ready, paste0(group_name, "_exact28_network_weight.csv")), row.names = FALSE)
  model_summary[[group_name]] <- data.frame(
    group = group_name,
    cells = ncol(object@data),
    genes = nrow(object@data),
    signaling_genes = nrow(object@data.signaling),
    observed_cell_types = nlevels(object@idents),
    significant_interactions = nrow(comm),
    pathways = length(object@netP$pathways),
    stringsAsFactors = FALSE
  )
  rm(object)
  gc()
}
write.table(do.call(rbind, model_summary), file.path(plot_ready, "M03bcd_model_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("CellChat complete: 28 cell classes
")
