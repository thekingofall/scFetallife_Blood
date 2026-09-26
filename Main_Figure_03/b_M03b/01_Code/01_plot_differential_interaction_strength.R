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
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize); library(grid); library(jsonlite)})
root <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
figures <- file.path(root, "02_Figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(file.path(.scf_shared, "Main_Figure_03/b_M03b/03_PlotData/M03b_differential_interaction_matrix.csv"), check.names = FALSE)
order_table <- read.delim(file.path(.scf_project_root, "Main_Figure_03/00_CellChat_Analysis/02_Data/03_Plot_Ready/M03bcd_exact_28class_celltype_order.tsv"))
display_order <- order_table$display_cell_type
model_order <- order_table$model_cell_type
cell_colors <- unlist(fromJSON(file.path(.scf_shared, "Main_Figure_03/a_M03a/M03_Cell_Type_Colors.json")))[model_order]
names(cell_colors) <- display_order
weight_diff <- matrix(0, length(display_order), length(display_order), dimnames = list(display_order, display_order))
d <- d[d$measure == "weight", ]
weight_diff[cbind(match(d$source, display_order), match(d$target, display_order))] <- d$difference
make_heatmap <- function(mat, title, legend_title) {
  limit <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) limit <- 1
  colors <- colorRamp2(c(-limit, 0, limit), c("#2B5C8A", "#F7F5EF", "#A6423A"))
  top_values <- colSums(abs(mat))
  right_values <- rowSums(abs(mat))
  top_colors <- cell_colors[colnames(mat)]
  right_colors <- cell_colors[rownames(mat)]
  top_annotation <- HeatmapAnnotation(
    Strength = anno_barplot(top_values, gp = gpar(fill = top_colors, col = "black", lwd = 0.35), border = TRUE, axis_param = list(gp = gpar(fontfamily = "Arial", fontsize = 9))),
    show_annotation_name = FALSE,
    height = unit(12, "mm")
  )
  right_annotation <- rowAnnotation(
    Strength = anno_barplot(right_values, gp = gpar(fill = right_colors, col = "black", lwd = 0.35), border = TRUE, axis_param = list(gp = gpar(fontfamily = "Arial", fontsize = 9))),
    show_annotation_name = FALSE,
    width = unit(12, "mm")
  )
  bottom_annotation <- HeatmapAnnotation(
    Cell_type = factor(colnames(mat), levels = names(cell_colors)),
    col = list(Cell_type = cell_colors),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm")
  )
  left_annotation <- rowAnnotation(
    Cell_type = factor(rownames(mat), levels = names(cell_colors)),
    col = list(Cell_type = cell_colors),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm")
  )
  Heatmap(
    mat,
    name = legend_title,
    col = colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.35),
    top_annotation = top_annotation,
    right_annotation = right_annotation,
    bottom_annotation = bottom_annotation,
    left_annotation = left_annotation,
    column_title = title,
    row_title = "Sources (sender)",
    column_names_rot = 90,
    row_names_side = "left",
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 10.5, fontface = "plain"),
    row_names_gp = gpar(fontfamily = "Arial", fontsize = 10.5, fontface = "plain"),
    column_title_gp = gpar(fontfamily = "Arial", fontsize = 14, fontface = "plain"),
    row_title_gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "plain"),
    heatmap_legend_param = list(
      title_gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "plain"),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 11),
      border = "black",
      legend_height = unit(35, "mm")
    )
  )
}

save_heatmap_pair <- function(plot, stem, width, height, dpi = 300) {
  pdf_path <- file.path(figures, paste0(stem, ".pdf"))
  png_path <- file.path(figures, paste0(stem, ".png"))
  cairo_pdf(pdf_path, width = width, height = height, family = "Arial", onefile = TRUE)
  draw(plot, heatmap_legend_side = "right")
  dev.off()
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi, background = "white")
  } else {
    png(png_path, width = width, height = height, units = "in", res = dpi, type = "cairo", family = "Arial")
  }
  draw(plot, heatmap_legend_side = "right")
  dev.off()
}

ht_weight <- make_heatmap(weight_diff, "Differential interaction strength", "Relative values")
save_heatmap_pair(ht_weight, "M03b_Differential_interaction_strength", 11.8, 10.2)
