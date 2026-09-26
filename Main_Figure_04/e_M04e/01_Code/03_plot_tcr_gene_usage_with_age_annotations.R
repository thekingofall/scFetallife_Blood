.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("circlize", "ComplexHeatmap", "grid"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_panel_root <- normalizePath(file.path(.scf_code_dir, ".."), winslash = "/", mustWork = FALSE)
.scf_min <- file.path(.scf_panel_root, "03_PlotData")
.scf_fig <- file.path(.scf_panel_root, "02_Figures")
dir.create(.scf_fig, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

if (!"Arial" %in% names(pdfFonts())) {
  pdfFonts(Arial = Type1Font("Arial", pdfFonts("Helvetica")[[1]]$metrics))
}

plot_data <- read.csv(file.path(.scf_min, "M04e_gene_usage_plot_data_with_gestational_week.csv"), check.names = FALSE, stringsAsFactors = FALSE)
age_stats <- read.csv(file.path(.scf_shared, "Main_Figure_04/e_M04e/03_PlotData/M04e_gene_usage_gestational_week_spearman_statistics.csv"), check.names = FALSE, stringsAsFactors = FALSE)
scale_ranges <- read.csv(file.path(.scf_shared, "Main_Figure_04/e_M04e/03_PlotData/M04e_segment_independent_scale_ranges_PBMC.csv"), check.names = FALSE, stringsAsFactors = FALSE)

segments <- c("TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ")
cell_types <- c("Naïve CD4 T", "Naïve CD8 T")
stopifnot(
  identical(sort(unique(plot_data$Organ)), "PBMC"),
  identical(sort(unique(plot_data$CellType)), sort(cell_types)),
  all(segments %in% unique(plot_data$segment))
)

gene_orders <- setNames(lapply(segments, function(seg) {
  x <- plot_data[plot_data$segment == seg, , drop = FALSE]
  order0 <- unique(x$variable)
  totals <- tapply(x$value, x$variable, sum, na.rm = TRUE)
  order0[is.finite(totals[order0]) & totals[order0] > 0]
}), segments)

heat_colors <- c("#313695", "#649AC7", "#FFFFBF", "#FDBE70", "#EA5839", "#A50026")
scale_max <- setNames(scale_ranges$scale_max_percent, scale_ranges$segment)
color_functions <- setNames(lapply(segments, function(seg) {
  colorRamp2(seq(0, scale_max[[seg]], length.out = length(heat_colors)), heat_colors)
}), segments)

sample_meta <- unique(plot_data[, c("MainID", "Gestational_week")])
sample_meta <- sample_meta[order(sample_meta$Gestational_week, sample_meta$MainID), , drop = FALSE]

sample_annotation <- function() {
  id_labels <- sub("_P[0-9]+$", "", sample_meta$MainID)
  rowAnnotation(
    ID = anno_text(
      id_labels,
      which = "row",
      location = unit(1, "npc"),
      just = "right",
      gp = gpar(fontfamily = "Arial", fontsize = 7, fontface = "bold", col = "black"),
      width = unit(13, "mm")
    ),
    show_annotation_name = FALSE
  )
}

age_annotation <- function(seg, cell_type) {
  ord <- gene_orders[[seg]]
  s <- age_stats[age_stats$segment == seg & age_stats$CellType == cell_type, , drop = FALSE]
  s <- s[match(ord, s$variable), , drop = FALSE]
  marker <- ifelse(
    !is.na(s$p_value) & s$p_value < 0.05 & s$direction == "positive", "\u2191",
    ifelse(!is.na(s$p_value) & s$p_value < 0.05 & s$direction == "negative", "\u2193", "")
  )
  HeatmapAnnotation(
    pcw = anno_text(
      marker,
      gp = gpar(fontfamily = "Arial", fontsize = 16, fontface = "bold", col = "black"),
      location = unit(0.8, "mm"),
      just = "bottom",
      rot = 0
    ),
    show_annotation_name = FALSE,
    height = unit(5.5, "mm")
  )
}

segment_matrix <- function(seg, cell_type) {
  ord <- gene_orders[[seg]]
  d <- plot_data[plot_data$segment == seg & plot_data$CellType == cell_type & plot_data$variable %in% ord, , drop = FALSE]
  mat <- matrix(0, nrow = nrow(sample_meta), ncol = length(ord), dimnames = list(sample_meta$MainID, ord))
  ri <- match(d$MainID, sample_meta$MainID)
  ci <- match(d$variable, ord)
  mat[cbind(ri, ci)] <- as.numeric(d$value)
  mat
}

segment_heatmap <- function(seg, cell_type, show_title, show_column_names) {
  Heatmap(
    segment_matrix(seg, cell_type),
    name = seg,
    col = color_functions[[seg]],
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = NULL,
    column_title = if (show_title) seg else NULL,
    column_title_gp = gpar(fontfamily = "Arial", fontsize = 19, fontface = "bold", col = "black"),
    top_annotation = age_annotation(seg, cell_type),
    show_row_names = FALSE,
    show_column_names = show_column_names,
    column_names_rot = 90,
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 11, fontface = "bold", col = "black"),
    column_names_max_height = unit(34, "mm"),
    rect_gp = gpar(col = NA),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 1.5),
    show_heatmap_legend = FALSE
  )
}

assemble_block <- function(cell_type, show_title, show_column_names) {
  ht <- sample_annotation() + segment_heatmap(segments[[1]], cell_type, show_title, show_column_names)
  for (seg in segments[-1]) ht <- ht + segment_heatmap(seg, cell_type, show_title, show_column_names)
  ht
}

compact_legends <- function() {
  legends <- setNames(lapply(segments, function(seg) {
    at <- pretty(c(0, scale_max[[seg]]), n = 3)
    at <- at[at >= 0 & at <= scale_max[[seg]]]
    Legend(
      title = seg,
      col_fun = color_functions[[seg]],
      at = at,
      labels = if (max(at) >= 10) formatC(at, format = "f", digits = 0) else formatC(at, format = "f", digits = 1),
      title_position = "topleft",
      title_gp = gpar(fontfamily = "Arial", fontsize = 10, fontface = "bold", col = "black"),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 9, fontface = "bold", col = "black"),
      grid_width = unit(2.3, "mm"),
      legend_height = unit(12, "mm")
    )
  }), segments)
  row_1 <- packLegend(legends$TRAV, legends$TRAJ, direction = "horizontal", gap = unit(2, "mm"))
  row_2 <- packLegend(legends$TRBV, legends$TRBD, legends$TRBJ, direction = "horizontal", gap = unit(2, "mm"))
  packLegend(row_1, row_2, direction = "vertical", gap = unit(2.4, "mm"))
}

draw_block <- function(ht, label) {
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.048, 0.952), "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(label, x = 0.26, y = 0.5, gp = gpar(fontfamily = "Arial", fontsize = 17, fontface = "bold", col = "black"))
  grid.lines(
    x = unit(c(0.78, 0.78), "npc"),
    y = unit(c(0.88, 0.12), "npc"),
    gp = gpar(col = "black", fill = "black", lwd = 2.8, lineend = "round"),
    arrow = arrow(type = "closed", length = unit(3.2, "mm"))
  )
  upViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw(
    ht,
    newpage = FALSE,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    ht_gap = unit(1.2, "mm"),
    padding = unit(c(0.8, 0.5, 0.8, 0.5), "mm")
  )
  upViewport(2)
}

draw_figure <- function(top_ht, bottom_ht, legends) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    2,
    2,
    heights = unit(c(0.424, 0.576), "npc"),
    widths = unit(c(0.925, 0.075), "npc")
  )))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw_block(top_ht, "CD4")
  upViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw_block(bottom_ht, "CD8")
  upViewport()
  pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 2))
  grid.draw(legends)
  upViewport(2)
  grid.text("pcw  ID", x = unit(0.064, "npc"), y = unit(0.965, "npc"), gp = gpar(fontfamily = "Arial", fontsize = 12.5, fontface = "bold", col = "black"))
}

top_ht <- assemble_block(cell_types[[1]], TRUE, FALSE)
bottom_ht <- assemble_block(cell_types[[2]], FALSE, TRUE)
legends <- compact_legends()
basename <- "M04e_TRAV_TRAJ_TRBV_TRBD_TRBJ_gene_usage"
pdf_path <- file.path(.scf_fig, paste0(basename, ".pdf"))
png_path <- file.path(.scf_fig, paste0(basename, ".png"))

cairo_pdf(pdf_path, width = 22.5, height = 5.625, family = "Arial")
draw_figure(top_ht, bottom_ht, legends)
dev.off()
png(png_path, width = 4000, height = 1000, units = "px", res = 177.777, type = "cairo")
draw_figure(top_ht, bottom_ht, legends)
dev.off()


cat("Saved M04e TCR gene-usage heatmap\n")
