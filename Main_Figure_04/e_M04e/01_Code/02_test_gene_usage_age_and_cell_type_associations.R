.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("circlize", "ComplexHeatmap", "grid", "tidyverse"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")
.scf_input_dir <- file.path(.scf_shared, "Main_Figure_04_05/02_Gene_Usage")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

if (!"Arial" %in% names(pdfFonts())) {
  pdfFonts(Arial = Type1Font("Arial", pdfFonts("Helvetica")[[1]]$metrics))
}

options(stringsAsFactors = FALSE)


prepare_long <- function(tab, segment_name, cell_type_mode = c("keep", "merge_b")) {
  cell_type_mode <- match.arg(cell_type_mode)
  x <- tab %>% filter(segment == segment_name)
  if (cell_type_mode == "merge_b") {
    counts <- x %>%
      group_by(MainID, Main_Organ, Gestational_Age_Weeks, gene, segment) %>%
      summarise(count = sum(as.numeric(count)), .groups = "drop")
    denoms <- x %>%
      distinct(Last_cell_type, MainID, Main_Organ, Gestational_Age_Weeks, segment, segment_total) %>%
      group_by(MainID, Main_Organ, Gestational_Age_Weeks, segment) %>%
      summarise(segment_total = sum(as.numeric(segment_total)), .groups = "drop")
    x <- left_join(counts, denoms, by = c("MainID", "Main_Organ", "Gestational_Age_Weeks", "segment")) %>%
      mutate(
        percentage = ifelse(segment_total > 0, 100 * count / segment_total, 0),
        Last_cell_type = "Naïve B"
      )
  }
  x %>%
    transmute(
      Sample = paste0(Last_cell_type, "_", gsub("_", ".", MainID)),
      MainID = MainID,
      Gestational_week = as.numeric(Gestational_Age_Weeks),
      Organ = Main_Organ,
      CellType = Last_cell_type,
      segment = segment,
      variable = gene,
      value = as.numeric(percentage)
    )
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(n = length(x), rho = NA_real_, p_value = NA_real_))
  }
  z <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(n = length(x), rho = unname(z$estimate), p_value = z$p.value)
}

calculate_age_statistics <- function(long_data) {
  out <- long_data %>%
    group_by(segment, variable, CellType) %>%
    group_modify(~ safe_spearman(.x$Gestational_week, .x$value)) %>%
    ungroup() %>%
    group_by(segment, CellType) %>%
    mutate(p_adjust_bh = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      significant_raw_p_0_05 = !is.na(p_value) & p_value < 0.05,
      direction = case_when(rho > 0 ~ "positive", rho < 0 ~ "negative", TRUE ~ NA_character_),
      marker = case_when(significant_raw_p_0_05 & rho > 0 ~ "↑", significant_raw_p_0_05 & rho < 0 ~ "↓", TRUE ~ "")
    )
  out
}

calculate_group_comparison <- function(long_data, group_1, group_2) {
  genes <- long_data %>% distinct(segment, variable)
  out <- pmap_dfr(genes, function(segment, variable) {
    d <- long_data %>%
      filter(.data$segment == .env$segment, .data$variable == .env$variable) %>%
      select(MainID, CellType, value) %>%
      pivot_wider(names_from = CellType, values_from = value)
    x <- d[[group_1]]
    y <- d[[group_2]]
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    paired_p <- if (length(x) > 0 && any(x != y)) suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value) else 1
    unpaired_p <- if (length(x) > 0 && any(x != y)) suppressWarnings(wilcox.test(x, y, paired = FALSE, exact = FALSE)$p.value) else 1
    tibble(
      segment = segment,
      variable = variable,
      group_1 = group_1,
      group_2 = group_2,
      n_pairs = length(x),
      median_group_1 = if (length(x)) median(x) else NA_real_,
      median_group_2 = if (length(y)) median(y) else NA_real_,
      median_paired_difference_group_2_minus_group_1 = if (length(x)) median(y - x) else NA_real_,
      p_value_paired_wilcoxon = paired_p,
      p_value_unpaired_wilcoxon_sensitivity = unpaired_p
    )
  }) %>%
    mutate(p_adjust_bh_paired = p.adjust(p_value_paired_wilcoxon, method = "BH")) %>%
    mutate(
      significant_bh_q_0_05 = p_adjust_bh_paired < 0.05,
      marker = ifelse(significant_bh_q_0_05, "*", "")
    )
  out
}

organ_levels <- "PBMC"
organ_colors <- c(PBMC = "#F6313E")
cell_colors <- c("Naïve CD4 T" = "#fa6e01", "Naïve CD8 T" = "#993333", "CXCR5+ Naïve B" = "#0081C9", "CXCR5- Naïve B" = "#E64B35", "Naïve B" = "#0081C9")
heat_colors <- c("#313695", "#649AC7", "#FFFFBF", "#FDBE70", "#EA5839", "#A50026")
segment_scale_max_overrides <- numeric()

nonzero_gene_orders <- function(long_data, orders) {
  map(orders, function(ord) {
    keep <- long_data %>%
      filter(variable %in% ord) %>%
      group_by(variable) %>%
      summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
      filter(total != 0) %>%
      pull(variable)
    ord[ord %in% keep]
  })
}

make_row_annotation <- function(meta, cell_type) {
  ann <- data.frame(
    CellType = factor(rep(cell_type, nrow(meta)), levels = names(cell_colors)),
    Organ = factor(meta$Organ, levels = organ_levels)
  )
  rowAnnotation(
    df = ann,
    col = list(CellType = cell_colors, Organ = organ_colors),
    show_legend = FALSE,
    annotation_name_gp = gpar(fontfamily = "Arial", fontsize = 11),
    simple_anno_size = unit(3.6, "mm")
  )
}

make_marker_annotation <- function(labels, location = c("top", "bottom")) {
  location <- match.arg(location)
  if (location == "top") {
    HeatmapAnnotation(
      Age = anno_text(labels, gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold"), just = "center", rot = 0),
      show_annotation_name = FALSE,
      height = unit(4.5, "mm"),
      which = "column"
    )
  } else {
    HeatmapAnnotation(
      Difference = anno_text(labels, gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold"), just = "center", rot = 0),
      show_annotation_name = FALSE,
      height = unit(4.5, "mm"),
      which = "column"
    )
  }
}

segment_color_fun <- function(long_data, segment_name) {
  max_value <- max(long_data$value[long_data$segment == segment_name], na.rm = TRUE)
  if (segment_name %in% names(segment_scale_max_overrides)) max_value <- unname(segment_scale_max_overrides[[segment_name]])
  if (!is.finite(max_value) || max_value <= 0) max_value <- 1
  colorRamp2(seq(0, max_value, length.out = length(heat_colors)), heat_colors)
}

segment_scale_table <- function(long_data, orders) {
  imap_dfr(orders, function(ord, seg) {
    d <- long_data %>% filter(segment == seg, variable %in% ord)
    data_max <- max(d$value, na.rm = TRUE)
    raw_data_max <- if ("Raw_value" %in% names(d)) max(d$Raw_value, na.rm = TRUE) else data_max
    scale_max <- if (seg %in% names(segment_scale_max_overrides)) unname(segment_scale_max_overrides[[seg]]) else data_max
    tibble(
      segment = seg,
      organ = paste(sort(unique(d$Organ)), collapse = ";"),
      n_cell_types = n_distinct(d$CellType),
      n_samples = n_distinct(d$MainID),
      n_displayed_genes = length(ord),
      scale_min_percent = 0,
      raw_data_max_percent = raw_data_max,
      data_max_percent = data_max,
      scale_max_percent = scale_max
    )
  })
}

make_external_legends <- function(long_data, orders, displayed_cell_types, max_height_in, include_context = TRUE) {
  heat_legends <- imap(orders, function(ord, seg) {
    col_fun <- segment_color_fun(long_data, seg)
    max_value <- max(long_data$value[long_data$segment == seg], na.rm = TRUE)
    if (seg %in% names(segment_scale_max_overrides)) max_value <- unname(segment_scale_max_overrides[[seg]])
    at <- unique(pretty(c(0, max_value), n = 3))
    labels <- if (max(at) >= 10) formatC(at, format = "f", digits = 0) else formatC(at, format = "f", digits = 1)
    Legend(
      title = seg,
      col_fun = col_fun,
      at = at,
      labels = labels,
      title_gp = gpar(fontfamily = "Arial", fontsize = 9.5, fontface = "bold"),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 8.5),
      grid_width = unit(2.6, "mm"),
      grid_height = unit(2, "mm")
    )
  })
  if (length(heat_legends) == 5 && all(c("TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ") %in% names(heat_legends))) {
    legend_col_1 <- packLegend(heat_legends[["TRAV"]], heat_legends[["TRBD"]], direction = "vertical", gap = unit(1.2, "mm"))
    legend_col_2 <- packLegend(heat_legends[["TRAJ"]], heat_legends[["TRBJ"]], direction = "vertical", gap = unit(1.2, "mm"))
    legend_col_3 <- packLegend(heat_legends[["TRBV"]], direction = "vertical")
    segment_legends <- packLegend(legend_col_1, legend_col_2, legend_col_3, direction = "horizontal", gap = unit(2.2, "mm"))
  } else if (length(heat_legends) == 7 && all(c("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ") %in% names(heat_legends))) {
    legend_col_1 <- packLegend(heat_legends[["IGHV"]], heat_legends[["IGKV"]], heat_legends[["IGLV"]], direction = "vertical", gap = unit(1.2, "mm"))
    legend_col_2 <- packLegend(heat_legends[["IGHD"]], heat_legends[["IGKJ"]], heat_legends[["IGLJ"]], direction = "vertical", gap = unit(1.2, "mm"))
    legend_col_3 <- packLegend(heat_legends[["IGHJ"]], direction = "vertical")
    segment_legends <- packLegend(legend_col_1, legend_col_2, legend_col_3, direction = "horizontal", gap = unit(2.2, "mm"))
  } else {
    segment_legends <- do.call(packLegend, c(heat_legends, list(direction = "vertical", max_height = unit(max_height_in, "in"), gap = unit(2, "mm"))))
  }
  cell_legend <- Legend(
    title = "Cell type",
    at = displayed_cell_types,
    legend_gp = gpar(fill = cell_colors[displayed_cell_types], col = NA),
    title_gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontfamily = "Arial", fontsize = 11),
    grid_width = unit(3.8, "mm"),
    grid_height = unit(3.8, "mm")
  )
  organ_legend <- Legend(
    title = "Organ",
    at = organ_levels,
    legend_gp = gpar(fill = organ_colors[organ_levels], col = NA),
    title_gp = gpar(fontfamily = "Arial", fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontfamily = "Arial", fontsize = 11),
    grid_width = unit(3.8, "mm"),
    grid_height = unit(3.8, "mm")
  )
  context_legends <- if (include_context) list(cell_legend, organ_legend) else list()
  if (length(context_legends) == 0) return(segment_legends)
  do.call(packLegend, c(context_legends, list(segment_legends), list(direction = "vertical", max_height = unit(max_height_in, "in"), gap = unit(2, "mm"))))
}

build_segment_heatmap <- function(long_data, segment_name, gene_order, cell_type, age_stats, comparison_stats = NULL, top_title = TRUE, bottom_names = FALSE, add_difference = FALSE, show_row_names = FALSE) {
  d <- long_data %>%
    filter(segment == segment_name, CellType == cell_type, variable %in% gene_order) %>%
    mutate(variable = factor(variable, levels = gene_order))
  meta <- d %>%
    distinct(Sample, MainID, Gestational_week, Organ, CellType) %>%
    mutate(Organ = factor(Organ, levels = organ_levels)) %>%
    arrange(Organ, Gestational_week, MainID)
  wide <- d %>%
    select(Sample, variable, value) %>%
    pivot_wider(names_from = variable, values_from = value, values_fill = 0)
  wide <- left_join(meta %>% select(Sample), wide, by = "Sample")
  mat <- as.matrix(wide[, gene_order, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- sub("\\.P[0-9]+$", "", gsub("_", ".", meta$MainID))
  age_labels <- age_stats %>%
    filter(segment == segment_name, CellType == cell_type) %>%
    right_join(tibble(variable = gene_order), by = "variable") %>%
    mutate(variable = factor(variable, levels = gene_order)) %>%
    arrange(variable) %>%
    pull(marker)
  age_labels[is.na(age_labels)] <- ""
  top_ann <- make_marker_annotation(age_labels, "top")
  bottom_ann <- NULL
  if (add_difference) {
    difference_labels <- comparison_stats %>%
      filter(segment == segment_name) %>%
      right_join(tibble(variable = gene_order), by = "variable") %>%
      mutate(variable = factor(variable, levels = gene_order)) %>%
      arrange(variable) %>%
      pull(marker)
    difference_labels[is.na(difference_labels)] <- ""
    bottom_ann <- make_marker_annotation(difference_labels, "bottom")
  }
  ht <- Heatmap(
    mat,
    name = segment_name,
    col = segment_color_fun(long_data, segment_name),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = NULL,
    column_title = if (top_title) segment_name else NULL,
    column_title_gp = gpar(fontfamily = "Arial", fontsize = 14, fontface = "plain"),
    top_annotation = top_ann,
    bottom_annotation = bottom_ann,
    show_row_names = show_row_names,
    row_names_side = "left",
    row_names_gp = gpar(fontfamily = "Arial", fontsize = 10),
    show_column_names = bottom_names,
    column_names_rot = 90,
    column_names_gp = gpar(fontfamily = "Arial", fontsize = 10),
    rect_gp = gpar(col = NA),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 1.4),
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(
      title = paste0(segment_name, " (%)"),
      title_gp = gpar(fontfamily = "Arial", fontsize = 11),
      labels_gp = gpar(fontfamily = "Arial", fontsize = 10),
      grid_width = unit(3, "mm"),
      grid_height = unit(3, "mm")
    )
  )
  list(annotation = make_row_annotation(meta, cell_type), heatmap = ht)
}

assemble_tcr_list <- function(long_data, orders, cell_type, age_stats, comparison_stats, top_title, bottom_names, add_difference) {
  parts <- imap(orders, function(ord, seg) build_segment_heatmap(
    long_data, seg, ord, cell_type, age_stats, comparison_stats,
    top_title = top_title,
    bottom_names = bottom_names,
    add_difference = add_difference,
    show_row_names = FALSE
  ))
  ht <- parts[[1]]$annotation + parts[[1]]$heatmap
  if (length(parts) > 1) {
    for (i in 2:length(parts)) ht <- ht + parts[[i]]$heatmap
  }
  ht
}

draw_tcr_pair <- function(top_ht, bottom_ht, external_legends) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.46, 0.54), "npc"), widths = unit(c(0.84, 0.16), "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(top_ht, newpage = FALSE, show_heatmap_legend = FALSE, show_annotation_legend = FALSE, ht_gap = unit(7, "mm"), padding = unit(c(1, 1, 0.5, 1), "mm"))
  upViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw(bottom_ht, newpage = FALSE, show_heatmap_legend = FALSE, show_annotation_legend = FALSE, ht_gap = unit(7, "mm"), padding = unit(c(0.5, 1, 1, 1), "mm"))
  upViewport()
  pushViewport(viewport(layout.pos.row = 1:2, layout.pos.col = 2))
  grid.draw(external_legends)
  upViewport(2)
}

save_tcr_pair <- function(top_ht, bottom_ht, external_legends, pdf_path, png_path, width, height) {
  cairo_pdf(pdf_path, width = width, height = height, family = "Arial")
  draw_tcr_pair(top_ht, bottom_ht, external_legends)
  dev.off()
  png(png_path, width = width, height = height, units = "in", res = 200, type = "cairo")
  draw_tcr_pair(top_ht, bottom_ht, external_legends)
  dev.off()
}

draw_single_heatmap <- function(ht, external_legends) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(0.84, 0.16), "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(ht, newpage = FALSE, show_heatmap_legend = FALSE, show_annotation_legend = FALSE, ht_gap = unit(3.2, "mm"), padding = unit(c(1, 1, 1, 1), "mm"))
  upViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid.draw(external_legends)
  upViewport(2)
}

save_single_heatmap <- function(ht, external_legends, pdf_path, png_path, width, height) {
  cairo_pdf(pdf_path, width = width, height = height, family = "Arial")
  draw_single_heatmap(ht, external_legends)
  dev.off()
  png(png_path, width = width, height = height, units = "in", res = 200, type = "cairo")
  draw_single_heatmap(ht, external_legends)
  dev.off()
}

tcr_raw <- read_csv(file.path(.scf_input_dir, "M04e_complete_per_sample_segment_gene_usage.csv"), show_col_types = FALSE) %>%
  filter(Main_Organ == "PBMC")

stopifnot(nrow(tcr_raw) > 0, identical(sort(unique(tcr_raw$Main_Organ)), "PBMC"))
tcr_gene_root <- file.path(.scf_shared, "00_Common/03_Receptor_Gene_Order", "TCR")
tra_locus <- read_csv(file.path(tcr_gene_root, "TRA_Gene_Locus.csv"), show_col_types = FALSE)
trb_locus <- read_csv(file.path(tcr_gene_root, "TRB_Gene_Locus.csv"), show_col_types = FALSE)
tcr_orders <- list(
  TRAV = tra_locus$IMGT_gene_name[grepl("^TRAV", tra_locus$IMGT_gene_name)],
  TRAJ = tra_locus$IMGT_gene_name[grepl("^TRAJ", tra_locus$IMGT_gene_name)],
  TRBV = trb_locus$IMGT_gene_name[grepl("^TRBV", trb_locus$IMGT_gene_name)],
  TRBD = trb_locus$IMGT_gene_name[grepl("^TRBD", trb_locus$IMGT_gene_name)],
  TRBJ = trb_locus$IMGT_gene_name[grepl("^TRBJ", trb_locus$IMGT_gene_name)]
)
tcr_long <- map_dfr(names(tcr_orders), ~ prepare_long(tcr_raw, .x, "keep"))
tcr_orders <- nonzero_gene_orders(tcr_long, tcr_orders)
tcr_age_stats <- calculate_age_statistics(tcr_long)
tcr_comparison_stats <- calculate_group_comparison(tcr_long, "Naïve CD4 T", "Naïve CD8 T")
tcr_scale_ranges <- segment_scale_table(tcr_long, tcr_orders)

m04e_root <- file.path(.scf_project_root, "Main_Figure_04", "e_M04e")
m04e_fig <- file.path(m04e_root, "02_Figures")
m04e_min <- file.path(m04e_root, "03_PlotData")
m04e_src <- file.path(m04e_root, "03_Source_Data")
dir.create(m04e_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(m04e_min, recursive = TRUE, showWarnings = FALSE)
dir.create(m04e_src, recursive = TRUE, showWarnings = FALSE)

tcr_top <- assemble_tcr_list(tcr_long, tcr_orders, "Naïve CD4 T", tcr_age_stats, tcr_comparison_stats, TRUE, FALSE, FALSE)
tcr_bottom <- assemble_tcr_list(tcr_long, tcr_orders, "Naïve CD8 T", tcr_age_stats, tcr_comparison_stats, FALSE, TRUE, TRUE)
tcr_legends <- make_external_legends(tcr_long, tcr_orders, c("Naïve CD4 T", "Naïve CD8 T"), 3.25, include_context = FALSE)
save_tcr_pair(
  tcr_top,
  tcr_bottom,
  tcr_legends,
  file.path(m04e_fig, "M04e_TRAV_TRAJ_TRBV_TRBD_TRBJ_gene_usage.pdf"),
  file.path(m04e_fig, "M04e_TRAV_TRAJ_TRBV_TRBD_TRBJ_gene_usage.png"),
  20,
  7.5
)

tcr_subset_specs <- list(
  TRAVJ = c("TRAV", "TRAJ"),
  TRBVDJ = c("TRBV", "TRBD", "TRBJ")
)
for (nm in names(tcr_subset_specs)) {
  sub_orders <- tcr_orders[tcr_subset_specs[[nm]]]
  top_ht <- assemble_tcr_list(tcr_long, sub_orders, "Naïve CD4 T", tcr_age_stats, tcr_comparison_stats, TRUE, FALSE, FALSE)
  bottom_ht <- assemble_tcr_list(tcr_long, sub_orders, "Naïve CD8 T", tcr_age_stats, tcr_comparison_stats, FALSE, TRUE, TRUE)
  legends <- make_external_legends(tcr_long, sub_orders, c("Naïve CD4 T", "Naïve CD8 T"), 3.25, include_context = FALSE)
  width <- min(13.2, max(7.5, sum(lengths(sub_orders)) / 11))
  save_tcr_pair(top_ht, bottom_ht, legends, file.path(m04e_fig, paste0("M04e_", nm, "_gene_usage.pdf")), file.path(m04e_fig, paste0("M04e_", nm, "_gene_usage.png")), width, 5.2)
}
for (seg in names(tcr_orders)) {
  sub_orders <- tcr_orders[seg]
  top_ht <- assemble_tcr_list(tcr_long, sub_orders, "Naïve CD4 T", tcr_age_stats, tcr_comparison_stats, TRUE, FALSE, FALSE)
  bottom_ht <- assemble_tcr_list(tcr_long, sub_orders, "Naïve CD8 T", tcr_age_stats, tcr_comparison_stats, FALSE, TRUE, TRUE)
  legends <- make_external_legends(tcr_long, sub_orders, c("Naïve CD4 T", "Naïve CD8 T"), 3.25, include_context = FALSE)
  width <- min(13.2, max(6.5, length(tcr_orders[[seg]]) / 11))
  save_tcr_pair(top_ht, bottom_ht, legends, file.path(m04e_fig, paste0("M04e_", seg, "_gene_usage.pdf")), file.path(m04e_fig, paste0("M04e_", seg, "_gene_usage.png")), width, 5.2)
}

write_csv(tcr_long, file.path(m04e_min, "M04e_gene_usage_plot_data_with_gestational_week.csv"))
write_csv(tcr_long %>% select(Sample, variable, value), file.path(.scf_shared, "Main_Figure_04/e_M04e/03_PlotData/M04e_gene_usage_plot_data.csv"))
write_csv(tcr_scale_ranges, file.path(.scf_shared, "Main_Figure_04/e_M04e/03_PlotData/M04e_segment_independent_scale_ranges_PBMC.csv"))
write_csv(tcr_age_stats, file.path(.scf_shared, "Main_Figure_04/e_M04e/03_PlotData/M04e_gene_usage_gestational_week_spearman_statistics.csv"))
write_csv(tcr_comparison_stats, file.path(.scf_shared, "Main_Figure_04/e_M04e/03_PlotData/M04e_CD4_vs_CD8_paired_wilcoxon_BH_statistics.csv"))

bcr_raw <- read_csv(file.path(.scf_input_dir, "M05e_complete_per_sample_segment_gene_usage.csv"), show_col_types = FALSE) %>%
  filter(Main_Organ == "PBMC")

stopifnot(nrow(bcr_raw) > 0, identical(sort(unique(bcr_raw$Main_Organ)), "PBMC"))
bcr_gene_root <- file.path(.scf_shared, "00_Common/03_Receptor_Gene_Order", "BCR")
read_gene_order <- function(file, sep = "\t") {
  tab <- read.table(file.path(bcr_gene_root, file), header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  as.character(tab[[1]])
}
bcr_orders <- list(
  IGHV = read_gene_order("IGHV.txt"),
  IGHD = read_gene_order("IGHD.csv", ","),
  IGHJ = paste0("IGHJ", 1:6),
  IGKV = read_gene_order("IGKV.txt"),
  IGKJ = read_gene_order("IGKJ.txt"),
  IGLV = read_gene_order("IGLV.txt"),
  IGLJ = read_gene_order("IGLJ.txt")
)
bcr_subtype_long <- map_dfr(names(bcr_orders), ~ prepare_long(bcr_raw, .x, "keep"))
bcr_long <- map_dfr(names(bcr_orders), ~ prepare_long(bcr_raw, .x, "merge_b"))
bcr_orders <- nonzero_gene_orders(bcr_long, bcr_orders)
bcr_age_stats <- calculate_age_statistics(bcr_long)
bcr_comparison_stats <- calculate_group_comparison(bcr_subtype_long, "CXCR5+ Naïve B", "CXCR5- Naïve B")
segment_scale_max_overrides <- c(IGHV = 15, IGHD = 30, IGHJ = 60, IGKV = 30, IGKJ = 40, IGLV = 20, IGLJ = 80)
bcr_segment_medians <- bcr_long %>% group_by(segment) %>% summarise(segment_median = median(value, na.rm = TRUE), .groups = "drop")
bcr_plot_long <- bcr_long %>%
  left_join(bcr_segment_medians, by = "segment") %>%
  mutate(
    Raw_value = value,
    value = case_when(
      segment == "IGHV" & value >= 25 ~ segment_median,
      segment == "IGHD" & value > 51 ~ 30,
      segment == "IGLV" & (value == 100 | value > 50) ~ segment_median,
      segment == "IGKV" & (value == 100 | value > 80) ~ segment_median,
      segment == "IGKJ" & value == 100 ~ segment_median,
      TRUE ~ value
    )
  ) %>%
  select(-segment_median)
bcr_scale_ranges <- segment_scale_table(bcr_plot_long, bcr_orders)

assemble_bcr_list <- function(orders, cell_type, top_title, bottom_names, add_difference) {
  parts <- imap(orders, function(ord, seg) build_segment_heatmap(
    bcr_plot_long, seg, ord, cell_type, bcr_age_stats, bcr_comparison_stats,
    top_title = top_title,
    bottom_names = bottom_names,
    add_difference = add_difference,
    show_row_names = identical(seg, names(orders)[[1]])
  ))
  ht <- parts[[1]]$heatmap
  if (length(parts) > 1) {
    for (i in 2:length(parts)) ht <- ht + parts[[i]]$heatmap
  }
  ht
}

m05e_root <- file.path(.scf_project_root, "Main_Figure_05", "e_M05e")
m05e_fig <- file.path(m05e_root, "02_Figures")
m05e_min <- file.path(m05e_root, "03_PlotData")
m05e_src <- file.path(m05e_root, "03_Source_Data")
dir.create(m05e_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(m05e_min, recursive = TRUE, showWarnings = FALSE)
dir.create(m05e_src, recursive = TRUE, showWarnings = FALSE)

bcr_ht <- assemble_bcr_list(bcr_orders, "Naïve B", TRUE, TRUE, FALSE)
bcr_legends <- make_external_legends(bcr_plot_long, bcr_orders, "Naïve B", 4.75, include_context = FALSE)
save_single_heatmap(bcr_ht, bcr_legends, file.path(m05e_fig, "M05e_BCR_all_gene_usage.pdf"), file.path(m05e_fig, "M05e_BCR_all_gene_usage.png"), 20, 6.5)
for (seg in names(bcr_orders)) {
  ht <- assemble_bcr_list(bcr_orders[seg], "Naïve B", TRUE, TRUE, FALSE)
  legends <- make_external_legends(bcr_plot_long, bcr_orders[seg], "Naïve B", 3.25, include_context = FALSE)
  width <- min(13.2, max(6.5, length(bcr_orders[[seg]]) / 11))
  save_single_heatmap(ht, legends, file.path(m05e_fig, paste0("M05e_", seg, "_gene_usage.pdf")), file.path(m05e_fig, paste0("M05e_", seg, "_gene_usage.png")), width, 3.8)
}

write_csv(bcr_plot_long %>% mutate(Display_value = value, value = Raw_value) %>% select(Sample, MainID, Gestational_week, Organ, CellType, segment, variable, value, Display_value), file.path(m05e_min, "M05e_gene_usage_plot_data_with_gestational_week.csv"))
write_csv(bcr_plot_long %>% mutate(Display_value = value, value = Raw_value) %>% select(Sample, variable, value, Display_value), file.path(m05e_min, "M05e_gene_usage_plot_data.csv"))
write_csv(bcr_long %>% select(Sample, variable, value), file.path(m05e_src, "M05e_gene_usage_plot_data.csv"))
write_csv(bcr_scale_ranges, file.path(.scf_shared, "Main_Figure_05/e_M05e/03_PlotData/M05e_segment_independent_scale_ranges_PBMC.csv"))
write_csv(bcr_age_stats, file.path(.scf_shared, "Main_Figure_05/e_M05e/03_PlotData/M05e_gene_usage_gestational_week_spearman_statistics.csv"))
write_csv(bcr_comparison_stats, file.path(.scf_shared, "Main_Figure_05/e_M05e/03_PlotData/M05e_CXCR5pos_vs_CXCR5neg_paired_wilcoxon_BH_statistics.csv"))

cat("M04e and M05e statistics and annotations complete\n")
