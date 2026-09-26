.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ComplexHeatmap", "circlize", "dplyr", "ragg"))
suppressPackageStartupMessages({library(ComplexHeatmap);library(circlize);library(dplyr);library(grid);library(ragg)})
root <- normalizePath(file.path(.scf_start_dir, "../../../.."), winslash = "/")
panel <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
plot_data <- file.path(panel, "03_PlotData")
dir.create(plot_data, showWarnings = FALSE)
figures <- file.path(panel, "02_Figures")
dir.create(figures, showWarnings = FALSE)
stage_columns <- c("16-22PCW", "22-28PCW", "28-34PCW", "34-40PCW")
stage_labels <- c("16-22 pcw", "22-28 pcw", "28-34 pcw", "34-40 pcw")
stage_colors <- c("#BB3630", "#777EB3", "#2B75B4", "#499F50")
cluster_colors <- c("#47B1A5", "#B66AAA", "#D73931", "#28537A", "#8B9249", "#499F50")

cm <- readRDS(file.path(plot_data, "PBMC_Bulk_16_Modules.rds"))
programs <- read.delim(file.path(plot_data, "PBMC_Bulk_Temporal_Programs.tsv"))
mapping <- data.frame(representative_GO = programs$GO_label)
selected <- cm$wide.res[cm$wide.res$cluster %in% programs$cluster, ]
selected$manuscript_cluster <- programs$program[match(selected$cluster, programs$cluster)]
ordered_chunks <- lapply(1:6, function(cluster_id) {
  chunk <- selected[selected$manuscript_cluster == cluster_id, , drop = FALSE]
  values <- as.matrix(chunk[, stage_columns, drop = FALSE])
  if (nrow(values) > 1L) {
    row_index <- hclust(dist(values), method = "complete")$order
    chunk <- chunk[row_index, , drop = FALSE]
  }
  chunk$plot_row_order_within_cluster <- seq_len(nrow(chunk))
  chunk
})
selected_ordered <- bind_rows(ordered_chunks)
selected_ordered$plot_row_order_global <- seq_len(nrow(selected_ordered))

mat <- as.matrix(selected_ordered[, stage_columns, drop = FALSE])
rownames(mat) <- selected_ordered$gene
split_factor <- factor(selected_ordered$manuscript_cluster, levels = 1:6)
cluster_counts <- table(split_factor)

right_annotation <- rowAnnotation(
  program = anno_block(
    panel_fun = function(index, nm) {
      cluster_id <- as.integer(nm)
      grid.rect(
        x = unit(2.6, "mm"),
        width = unit(5.2, "mm"),
        gp = gpar(fill = cluster_colors[[cluster_id]], col = NA)
      )
      grid.text(
        paste0("n:", length(index)),
        x = unit(2.6, "mm"),
        rot = 90,
        gp = gpar(col = "white", fontsize = 7.2, fontface = "bold")
      )
      grid.text(
        mapping$representative_GO[[cluster_id]],
        x = unit(7.2, "mm"),
        just = "left",
        gp = gpar(col = cluster_colors[[cluster_id]], fontsize = 8.3)
      )
    },
    width = unit(70, "mm")
  ),
  width = unit(70, "mm")
)

top_annotation <- HeatmapAnnotation(
  pcw = stage_labels,
  col = list(pcw = stats::setNames(stage_colors, stage_labels)),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(4.5, "mm"),
  gp = gpar(col = "white", lwd = 0.8)
)

heatmap <- Heatmap(
  mat,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#0DA9CE", "white", "#E74A32")),
  row_split = split_factor,
  row_title = NULL,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_labels = stage_labels,
  column_names_side = "top",
  column_names_rot = 52,
  column_names_gp = gpar(fontsize = 8.2),
  row_gap = unit(1.5, "mm"),
  border = TRUE,
  rect_gp = gpar(col = NA),
  top_annotation = top_annotation,
  right_annotation = right_annotation,
  width = unit(39, "mm"),
  use_raster = TRUE,
  raster_quality = 4,
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, -1, 0, 1, 2),
    labels_gp = gpar(fontsize = 7.5),
    title_gp = gpar(fontsize = 8.2, fontface = "bold"),
    legend_height = unit(28, "mm")
  )
)

draw_panel <- function() {
  draw(
    heatmap,
    heatmap_legend_side = "right",
    merge_legends = TRUE,
    padding = unit(c(7, 5, 5, 7), "mm")
  )
}

pdf_path <- file.path(figures, "S17a_Bulk_Temporal_Expression_Programs.pdf")
png_path <- file.path(figures, "S17a_Bulk_Temporal_Expression_Programs.png")

cairo_pdf(pdf_path, width = 7.2, height = 6.6, family = "Arial")
draw_panel()
dev.off()

ragg::agg_png(png_path, width = 3600, height = 3300, units = "px", res = 500, background = "white")
draw_panel()
dev.off()

write.csv(selected_ordered, file.path(plot_data, "PBMC_Bulk_Temporal_Heatmap_Values.csv"), row.names = FALSE)
write.csv(data.frame(program = 1:6, genes = as.integer(cluster_counts)), file.path(plot_data, "PBMC_Bulk_Selected_Program_Counts.csv"), row.names = FALSE)
