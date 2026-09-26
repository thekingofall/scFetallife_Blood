args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
script <- if (length(files)) tail(files, 1)[[1]] else sub("^--file=", "", args[[1]])
code_dir <- dirname(normalizePath(script))
root <- normalizePath(file.path(code_dir, "../../.."))
source(file.path(root, "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("grid", "png", "patchwork", "ragg"))
library(grid)

panels <- list(
  b = c("b_M08b", "01_plot_variance_explained.R", "p8b"),
  c = c("c_M08c", "01_plot_factor_age_associations.R", "p8c"),
  d = c("d_M08d", "01_plot_factor3_age_association.R", "F3scatter_plot"),
  e = c("e_M08e", "01_plot_factor3_feature_expression.R", "combined_plot3"),
  f = c("f_M08f", "01_plot_mofa_factor_go_enrichment.R", "combined_plot")
)
plots <- lapply(panels, function(panel) {
  environment <- new.env(parent = globalenv())
  source(file.path(root, "Main_Figure_08", panel[[1]], "01_Code", panel[[2]]), local = environment)
  environment[[panel[[3]]]]
})
schematic <- png::readPNG(file.path(root, "Main_Figure_08/a_M08a/02_Figures/M08a_multiomics_study_design.png"))
output <- file.path(code_dir, "../02_Figures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

layout <- list(
  a = c(0.035, 0.535, 0.325, 0.425),
  b = c(0.370, 0.705, 0.620, 0.265),
  c = c(0.380, 0.535, 0.340, 0.155),
  d = c(0.755, 0.535, 0.225, 0.155),
  e = c(0.005, 0.010, 0.525, 0.490),
  f = c(0.535, 0.010, 0.460, 0.490)
)

draw_figure <- function() {
  grid.newpage()
  for (letter in names(layout)) {
    box <- layout[[letter]]
    pushViewport(viewport(x = box[[1]], y = box[[2]], width = box[[3]], height = box[[4]], just = c("left", "bottom")))
    if (letter == "a") {
      width <- min(convertWidth(unit(1, "npc"), "in", valueOnly = TRUE),
                   convertHeight(unit(1, "npc"), "in", valueOnly = TRUE) * ncol(schematic) / nrow(schematic))
      grid.raster(schematic, width = unit(width, "in"), height = unit(width * nrow(schematic) / ncol(schematic), "in"))
    } else {
      plot <- plots[[letter]]
      grob <- if (inherits(plot, "patchwork")) patchwork::patchworkGrob(plot) else ggplot2::ggplotGrob(plot)
      grid.draw(grob)
    }
    popViewport()
    grid.text(letter, x = box[[1]], y = box[[2]] + box[[4]] + 0.006,
              just = c("left", "bottom"), gp = gpar(fontfamily = "Arial", fontsize = 22, fontface = "bold"))
  }
}

grDevices::cairo_pdf(file.path(output, "Figure_08.pdf"), width = 20, height = 21.5, family = "Arial")
draw_figure()
dev.off()
ragg::agg_png(file.path(output, "Figure_08.png"), width = 20, height = 21.5, units = "in", res = 160, background = "white")
draw_figure()
dev.off()
