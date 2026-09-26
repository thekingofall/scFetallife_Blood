.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "ggrastr", "grid", "jsonlite", "ragg", "readr", "stringr"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(
  Sys.getenv("SCF_PROJECT_ROOT", unset = file.path(.scf_code_dir, "../../../..")),
  winslash = "/",
  mustWork = FALSE
)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(readr)
  library(stringr)
  library(grid)
})

repo <- .scf_project_root
runner_dir <- Sys.getenv(
  "SCF_RASTER_RUNNER_DIR",
  unset = file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_02_03/03_PlotData")
)
output_root <- Sys.getenv(
  "SCF_RASTER_OUTPUT_ROOT",
  unset = repo
)
original_umap <- file.path(.scf_shared, "Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv")
original_colors <- file.path(.scf_shared, "Main_Figure_01/01_Cell_Annotation/M01_Cell_Type_Colors.json")
source(file.path(.scf_code_dir, "01_define_umap_plot_functions.R"))

all_data <- read_csv(original_umap, show_col_types = FALSE)
lineage_columns <- grep("^Cell_lineage", names(all_data), value = TRUE)
cellname_columns <- grep("^Cellname", names(all_data), value = TRUE)
if (!"Cell_lineage" %in% names(all_data)) {
  if (!length(lineage_columns)) stop("UMAP table is missing the Cell_lineage column")
  all_data$Cell_lineage <- all_data[[lineage_columns[[1L]]]]
}
if (!"Cellname" %in% names(all_data)) {
  if (!length(cellname_columns)) stop("UMAP table is missing the Cellname column")
  all_data$Cellname <- all_data[[cellname_columns[[1L]]]]
}
if (length(lineage_columns) > 1L &&
    !identical(as.character(all_data[[lineage_columns[[1L]]]]),
               as.character(all_data[[lineage_columns[[2L]]]]))) {
  stop("Duplicate Cell_lineage columns contain different values")
}

colors_dict <- unlist(fromJSON(original_colors))

mark_arrow_layers <- function(data) {
  xr <- range(data$UMAP1, finite = TRUE)
  yr <- range(data$UMAP2, finite = TRUE)
  x0 <- xr[1] + 0.04 * diff(xr)
  y0 <- yr[1] + 0.04 * diff(yr)
  list(
    annotate("segment", x = x0, xend = x0 + 0.10 * diff(xr), y = y0, yend = y0,
             arrow = arrow(length = unit(0.10, "cm"), type = "closed"), linewidth = 0.35),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + 0.10 * diff(yr),
             arrow = arrow(length = unit(0.10, "cm"), type = "closed"), linewidth = 0.35)
  )
}

save_pair <- function(plot, directory, stem, width = 8, height = 8) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(directory, paste0(stem, ".png")), plot, width = width, height = height,
         device = ragg::agg_png, dpi = 320, bg = "white", limitsize = FALSE)
}

make_grey_panel <- function(selected, title, legendcol) {
  if (!nrow(selected)) stop("No cells selected for ", title)
  create_umap_plot_grey(
    umap_data = selected, alldata = all_data, position = "bottom",
    legendcol = legendcol, colors_dict = colors_dict
  ) +
    mark_arrow_layers(all_data) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.text = element_blank(), axis.ticks = element_blank()
    )
}

write_panel_source <- function(selected, directory, stem) {
  keep <- intersect(
    c("Cellname", "UMAP1", "UMAP2", "MainID", "Gestational_Age_Weeks", "Post_Conception_Age_Weeks", "organ", "Cell_lineage",
      "Last_cell_type", "Last_cell_type_num", "TCRBCR"),
    names(selected)
  )
  data_dir <- file.path(.scf_shared, "Supplementary_Figures",
                        basename(dirname(dirname(directory))), basename(dirname(directory)))
  dir.create(data_dir, recursive=TRUE, showWarnings=FALSE)
  write_csv(selected[, keep, drop = FALSE], file.path(data_dir, paste0(stem, "_source.csv")))
}

organ_specs <- list(
  S02a = list(value = "PBMC", title = "Blood", legendcol = 3L),
  S02b = list(value = "Liver", title = "Liver", legendcol = 3L),
  S02c = list(value = "Thymus", title = "Thymus", legendcol = 3L),
  S02d = list(value = "Spleen", title = "Spleen", legendcol = 3L)
)
for (panel_id in names(organ_specs)) {
  spec <- organ_specs[[panel_id]]
  selected <- all_data[all_data$organ == spec$value, , drop = FALSE]
  directory <- file.path(output_root, "Supplementary_Figures/Supplementary_Figure_02", paste0(substr(panel_id,4,4),"_",panel_id), "02_Figures")
  plot <- make_grey_panel(selected, spec$title, spec$legendcol)
  save_pair(plot, directory, paste0(panel_id, "_", spec$value, "_UMAP"))
  write_panel_source(selected, directory, panel_id)
}

lineage_specs <- list(
  S03a = list(title = "Precursor", legendcol = 1L, keep = all_data$Cell_lineage == "PRECURSOR"),
  S03b = list(title = "B cells", legendcol = 1L, keep = all_data$Cell_lineage == "B_CELL"),
  S03c = list(title = "T/ILC", legendcol = 2L, keep = all_data$Cell_lineage == "T/ILC"),
  S03d = list(title = "NK", legendcol = 1L, keep = grepl(" NK", all_data$Last_cell_type)),
  S03e = list(title = "MYELOID", legendcol = 2L, keep = grepl("Mono|My|Mac", all_data$Last_cell_type)),
  S03f = list(title = "DC", legendcol = 2L, keep = all_data$Cell_lineage == "DC"),
  S03g = list(title = "MK/ERY", legendcol = 3L, keep = grepl("MK", all_data$Cell_lineage)),
  S03h = list(title = "Others", legendcol = 2L, keep = grepl("En|Others", all_data$Last_cell_type))
)
for (panel_id in names(lineage_specs)) {
  spec <- lineage_specs[[panel_id]]
  selected <- all_data[spec$keep, , drop = FALSE]
  safe_title <- gsub("[^A-Za-z0-9]+", "_", spec$title)
  directory <- file.path(output_root, "Supplementary_Figures/Supplementary_Figure_03", paste0(substr(panel_id,4,4),"_",panel_id), "02_Figures")
  plot <- make_grey_panel(selected, spec$title, spec$legendcol)
  save_pair(plot, directory, paste0(panel_id, "_", safe_title, "_UMAP"))
  write_panel_source(selected, directory, panel_id)
}

tcrbcr <- all_data
tcrbcr$TCRBCR <- factor(tcrbcr$TCRBCR, levels = c("BCR", "None", "TCR"))
tcrbcr_plot <- ggplot(tcrbcr, aes(x = UMAP1, y = UMAP2, color = TCRBCR)) +
  ggrastr::geom_point_rast(size = 0.01, alpha = 0.6, raster.dpi = 300) +
  scale_color_manual(values = c("BCR" = "#C71000FF", "None" = "#008EA0FF", "TCR" = "#FF6F00FF")) +
  coord_fixed() +
  mark_arrow_layers(tcrbcr) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) +
  labs(color = "") + ggtitle("TCRBCR") +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
tcrbcr_dir <- file.path(.scf_project_root, "Supplementary_Figures/Supplementary_Figure_03", "i_S03i", "02_Figures")
save_pair(tcrbcr_plot, tcrbcr_dir, "S03i_TCRBCR_UMAP")
write_panel_source(tcrbcr, tcrbcr_dir, "S03i")

