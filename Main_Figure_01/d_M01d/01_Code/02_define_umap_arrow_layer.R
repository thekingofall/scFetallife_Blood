.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "grid"))

geom_markArrow <- function(mapping = NULL, data = NULL, stat = "identity",
                           position = "identity", na.rm = FALSE,
                           show.legend = NA, inherit.aes = TRUE,
                           lineend = "butt",linejoin = "round",
                           rel.pos = 0.1,rel.len = 0.3,
                           arrow = grid::arrow(type = "closed",length = grid::unit(0.4,"cm")),
                           arrow.fill = NULL,colour = "black",
                           label = c("UMAP1","UMAP2"),label.size = 3,label.shift = c(0.025,0.025),
                           fontface = "bold.italic",tail.shift = 0,corner.pos = "left_b",
                           facet.vars = NULL,
                           ...) {

  if(!is.null(facet.vars)){
    data = data.frame(facet.vars)
  }


  ggplot2::layer(
    geom = GeomMarkArrow, mapping = mapping,
    data = data,
    stat = stat, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm,
                  lineend = lineend,
                  linejoin = linejoin,
                  rel.pos = rel.pos,
                  rel.len = rel.len,
                  arrow.fill = arrow.fill,
                  arrow = arrow,
                  label = label,
                  label.size = label.size,
                  label.shift = label.shift,
                  fontface = fontface,
                  tail.shift = tail.shift,
                  corner.pos = corner.pos,
                  colour = colour,
                  ...)
  )
}


GeomMarkArrow <- ggplot2::ggproto("GeomMarkArrow", ggplot2::Geom,

                                  non_missing_aes = c("linetype", "linewidth", "shape"),
                                  default_aes = aes(linewidth = 0.5, linetype = 1,
                                                    alpha = NA),
                                  draw_key = draw_key_path,
                                  draw_panel = function(data, panel_scales, coord,
                                                        rel.pos = 0.1,rel.len = 0.3,
                                                        arrow.fill = NULL,colour = "black",
                                                        arrow = grid::arrow(type = "closed",length = grid::unit(0.4,"cm")),
                                                        lineend = "butt", linejoin = "round",
                                                        label = c("UMAP1","UMAP2"),label.size = 3,
                                                        label.shift = c(0.025,0.025),
                                                        fontface = "bold.italic",tail.shift = 0,
                                                        corner.pos = "left_b") {


                                    coords <- coord$transform(data, panel_scales)


                                    if(corner.pos == "left_b"){
                                      x = c(rel.pos,rel.pos)
                                      xend = c(rel.pos + rel.len,rel.pos)
                                      y = c(rel.pos,rel.pos)
                                      yend = c(rel.pos,rel.pos + rel.len)
                                    }else if(corner.pos == "left_u"){
                                      x = c(rel.pos,rel.pos)
                                      xend = c(rel.pos + rel.len,rel.pos)
                                      y = 1 - c(rel.pos,rel.pos)
                                      yend = 1 - c(rel.pos,rel.pos + rel.len)
                                    }else if(corner.pos == "right_b"){
                                      x = 1 - c(rel.pos,rel.pos)
                                      xend = 1 - c(rel.pos + rel.len,rel.pos)
                                      y = c(rel.pos,rel.pos)
                                      yend = c(rel.pos,rel.pos + rel.len)
                                    }else if(corner.pos == "right_u"){
                                      x = 1 - c(rel.pos,rel.pos)
                                      xend = 1 - c(rel.pos + rel.len,rel.pos)
                                      y = 1 - c(rel.pos,rel.pos)
                                      yend = 1 - c(rel.pos,rel.pos + rel.len)
                                    }


                                    coords <- data.frame(
                                      x = x,xend = xend,
                                      y = y,yend = yend,
                                      PANEL = 1,group = c(1:2),
                                      linewidth = unique(coords$linewidth),
                                      colour = rep(colour,2),
                                      linetype = unique(coords$linetype),
                                      alpha = unique(coords$alpha),
                                      label = label)


                                    arrow.fill <- arrow.fill %||% coords$colour

                                    seg_grob <-
                                      grid::segmentsGrob(
                                        x0 = c(coords$x[1] - tail.shift,coords$x[2]),
                                        x1 = coords$xend,
                                        y0 = c(coords$y[1],coords$y[2] - tail.shift),
                                        y1 = coords$yend,
                                        gp = grid::gpar(
                                          col = ggplot2::alpha(coords$colour, coords$alpha),
                                          fill = ggplot2::alpha(arrow.fill, coords$alpha),
                                          lwd = coords$linewidth * .pt,
                                          lty = coords$linetype,
                                          lineend = lineend,
                                          linejoin = linejoin
                                        ),
                                        arrow = arrow)

                                    text_grob <-
                                      grid::textGrob(label = coords$label,
                                                     x = c(((coords$x + coords$xend)/2)[1],
                                                           ((coords$x + coords$xend)/2)[2] - label.shift[2]),
                                                     y = c(((coords$y + coords$yend)/2)[1] - label.shift[1],
                                                           ((coords$y + coords$yend)/2)[2]),
                                                     hjust = 0.5,
                                                     rot = c(0,90),
                                                     check.overlap = TRUE,
                                                     gp = grid::gpar(col = ggplot2::alpha(coords$colour, coords$alpha),
                                                                     fontsize = label.size*.pt,
                                                                     fontface = fontface))


                                    grobs <- grid::grobTree(seg_grob,text_grob)
                                    grid::gTree(children = grid::gList(grobs))
                                  })
