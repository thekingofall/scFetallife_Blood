.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("Cairo", "circlize", "ComplexHeatmap", "grid", "ragg"))

.scf_code_dir <- .scf_start_dir
project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ragg)
  library(Cairo)
})

args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
code_dir <- if (length(args)) dirname(normalizePath(sub("^--file=", "", args[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
project_root <- normalizePath(file.path(code_dir, "../../.."), winslash = "/", mustWork = FALSE)
panel_root <- normalizePath(file.path(code_dir, ".."), winslash = "/", mustWork = FALSE)
input_dir <- file.path(project_root, "00_Global_Data", "Shared_Inputs", "Main_Figure_02/c_M02c")
figure_dir <- file.path(panel_root, "02_Figures")
minimal_dir <- file.path(panel_root, "03_PlotData")
source_dir <- file.path(panel_root, "03_Source_Data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(minimal_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

canonical_order <- c("M1-ERY", "M2-Mono", "M3-MK", "M4-STEM", "M5-Bcell", "M6-NK", "M7-Tcell")
canonical_labels <- c(
  "M1-ERY" = "M1 Erythroid",
  "M2-Mono" = "M2 Monocyte",
  "M3-MK" = "M3 Megakaryocyte",
  "M4-STEM" = "M4 Stem/progenitor",
  "M5-Bcell" = "M5 B cell",
  "M6-NK" = "M6 NK cell",
  "M7-Tcell" = "M7 T cell"
)
canonical_colors <- c(
  "M1-ERY" = "#C95752",
  "M2-Mono" = "#D17B49",
  "M3-MK" = "#D6A44A",
  "M4-STEM" = "#6C5B7B",
  "M5-Bcell" = "#4C956C",
  "M6-NK" = "#2A9D8F",
  "M7-Tcell" = "#3E6C8E"
)
heat_fun <- colorRamp2(c(-1, 0, 1), c("#343A9A", "#18A6B4", "#F4D943"))
mapping_path <- file.path(source_dir, "M02c_fresh_to_canonical_module_mapping.csv")
module_mapping <- read.csv(mapping_path, check.names = FALSE, stringsAsFactors = FALSE)
module_mapping$FreshModule <- as.integer(module_mapping$FreshModule)
module_mapping$CanonicalOrder <- as.integer(module_mapping$CanonicalOrder)
mapping_key <- paste(module_mapping$Stage, module_mapping$FreshModule, sep = "::")
if (anyDuplicated(mapping_key)) stop("Module-to-lineage mapping contains duplicate stage/module keys")

render_stage <- function(stage) {
  key <- tolower(stage)
  matrix_path <- file.path(input_dir, paste0(key, "_local_correlation_z.csv.gz"))
  module_path <- file.path(input_dir, paste0(key, "_modules.csv"))
  mat <- as.matrix(read.csv(gzfile(matrix_path), row.names = 1, check.names = FALSE))
  modules <- read.csv(module_path, check.names = FALSE)
  modules$gene <- as.character(modules$gene)
  modules$module <- as.integer(modules$module)
  keep <- intersect(rownames(mat), modules$gene)
  mat <- mat[keep, keep, drop = FALSE]
  module <- modules$module[match(keep, modules$gene)]
  mat_plot <- mat
  mat_plot[mat_plot < -1] <- -1
  mat_plot[mat_plot > 1] <- 1
  module_ids <- sort(unique(module))
  stage_mapping <- module_mapping[module_mapping$Stage == stage, , drop = FALSE]
  missing_mapping <- setdiff(module_ids, stage_mapping$FreshModule)
  if (length(missing_mapping)) stop(paste(stage, "modules missing lineage labels:", paste(missing_mapping, collapse = ", ")))
  fresh_to_canonical <- setNames(stage_mapping$CanonicalModule, stage_mapping$FreshModule)
  canonical <- unname(fresh_to_canonical[as.character(module)])
  if (anyNA(canonical)) stop(paste(stage, "contains genes without a lineage-module assignment"))
  gene_order <- unlist(lapply(canonical_order, function(canonical_id) {
    index <- which(canonical == canonical_id)
    if (length(index) < 3) return(index)
    index[hclust(as.dist(1 - mat_plot[index, index, drop = FALSE]), method = "average")$order]
  }), use.names = FALSE)
  mat_plot <- mat_plot[gene_order, gene_order, drop = FALSE]
  module <- module[gene_order]
  canonical <- unname(fresh_to_canonical[as.character(module)])
  canonical_factor <- factor(canonical, levels = canonical_order)
  canonical_counts <- table(factor(canonical, levels = canonical_order))
  display_lookup <- setNames(
    paste0(unname(canonical_labels[canonical_order]), " (n = ", as.integer(canonical_counts), ")"),
    canonical_order
  )
  canonical_display <- factor(
    unname(display_lookup[canonical]),
    levels = unname(display_lookup[canonical_order])
  )
  annotation <- rowAnnotation(
    Module = canonical_factor,
    col = list(Module = canonical_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(5, "mm")
  )
  heatmap <- Heatmap(
    mat_plot,
    name = "Z-score",
    col = heat_fun,
    row_split = canonical_display,
    row_order = seq_len(nrow(mat_plot)),
    column_order = seq_len(ncol(mat_plot)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title_side = "left",
    row_title_gp = gpar(fontfamily = "Arial", fontsize = 14, fontface = "bold", col = "#202020"),
    row_title_rot = 0,
    row_gap = unit(1.0, "mm"),
    border = "#202020",
    rect_gp = gpar(col = NA),
    use_raster = TRUE,
    raster_quality = 1.8,
    width = unit(132, "mm"),
    height = unit(132, "mm"),
    heatmap_legend_param = list(
      title = "Z-score",
      at = c(-1, -0.5, 0, 0.5, 1),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 13, col = "#202020"),
      title_gp = gpar(fontfamily = "Arial", fontsize = 14, fontface = "bold", col = "#202020"),
      legend_height = unit(34, "mm")
    )
  )
  draw_panel <- function() {
    draw(annotation + heatmap, heatmap_legend_side = "right", padding = unit(c(5, 9, 5, 6), "mm"))
  }
  pdf_path <- file.path(figure_dir, paste0("M02c_", stage, "_local_correlations.pdf"))
  png_path <- file.path(figure_dir, paste0("M02c_", stage, "_local_correlations.png"))
  CairoPDF(pdf_path, width = 9.4, height = 6.7, family = "Arial", bg = "white")
  draw_panel()
  dev.off()
  agg_png(png_path, width = 2820, height = 2010, units = "px", res = 300, background = "white")
  draw_panel()
  dev.off()
  count_table <- data.frame(
    Stage = stage,
    CanonicalModule = canonical_order,
    CanonicalLabel = unname(canonical_labels[canonical_order]),
    GeneCount = as.integer(canonical_counts),
    ContributingFreshModules = vapply(canonical_order, function(canonical_id) {
      paste(stage_mapping$FreshModule[stage_mapping$CanonicalModule == canonical_id], collapse = ",")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(count_table, file.path(minimal_dir, paste0("M02c_", stage, "_module_counts.csv")), row.names = FALSE)
  data.frame(
    Stage = stage,
    Gene = rownames(mat_plot),
    FreshModule = module,
    CanonicalModule = canonical,
    CanonicalLabel = unname(canonical_labels[canonical]),
    CanonicalOrder = match(canonical, canonical_order),
    PlotOrder = seq_len(nrow(mat_plot)),
    stringsAsFactors = FALSE
  )
}

early_order <- render_stage("Early")
late_order <- render_stage("Late")
write.csv(rbind(early_order, late_order), file.path(source_dir, "M02c_gene_module_plot_order.csv"), row.names = FALSE)
