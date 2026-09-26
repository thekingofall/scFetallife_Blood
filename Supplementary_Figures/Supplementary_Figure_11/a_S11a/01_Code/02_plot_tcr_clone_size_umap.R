.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "grDevices", "grid", "scales"))

library(ggplot2)
s11a_panel <- normalizePath(file.path(.scf_start_dir,".."),winslash="/")
dir.create(file.path(s11a_panel,"02_Figures"),showWarnings=FALSE,recursive=TRUE)
umap <- read.csv(file.path(s11a_panel,"03_PlotData/S11a_TCR_clone_annotations.csv.gz"))
save_plot_r <- function(plot, pdf_path, png_path, width_in, height_in) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()
  grDevices::png(
    png_path,
    width = width_in,
    height = height_in,
    units = "in",
    res = 300,
    type = "cairo",
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}

umap$clonal_expansion <- factor(as.character(umap$clonal_expansion), levels = c("1", "2", ">= 3"))
plot_source_umap <- umap[, c("cell_id", "UMAP1", "UMAP2", "clonal_expansion", "MainID", "Main_Organ")]
write.csv(
  plot_source_umap,
  gzfile(file.path(s11a_panel, "03_PlotData", "S11a_TCR_UMAP_plot_source_data.csv.gz")),
  row.names = FALSE,
  na = ""
)


umap <- umap[order(umap$clonal_expansion), , drop = FALSE]
umap_plot <- ggplot(umap, aes(x = UMAP1, y = UMAP2, colour = clonal_expansion)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  scale_colour_manual(values = c("1" = "#868686", "2" = "#145096", ">= 3" = "#E9412F"), drop = FALSE) +
  coord_equal() +
  labs(colour = "Cells per exact clonotype") +
  theme_void(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    legend.position = "right",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    legend.key.height = grid::unit(8, "pt"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.7, alpha = 1)))
save_plot_r(
  umap_plot,
  file.path(s11a_panel, "02_Figures", "S11a_TCR_categorical_clone_size_UMAP.pdf"),
  file.path(s11a_panel, "02_Figures", "S11a_TCR_categorical_clone_size_UMAP.png"),
  4.6,
  4.0
)
