.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "ggrastr", "ggrepel", "grid", "stringr"))

library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(grid)

draw_number_circle <- function(data, params, size) {
  grid::grobTree(
    pointsGrob(
      x = 0.5, y = 0.5,
      size = unit(2.8, "char"),
      pch = 16,
      gp = gpar(
        col = alpha(data$colour %||% "grey50", data$alpha),
        fill = alpha(data$fill %||% "grey50", data$alpha),
        lwd = (data$linewidth %||% 0.5) * .pt,
        lty = data$linetype %||% 1
      )
    ),
    textGrob(
      label = data$label,
      x = rep(0.5, 3), y = rep(0.5, 3),
      gp = gpar(col = "black", fontfamily = "Arial", fontsize = 16)
    )
  )
}

create_umap_plot <- function(umap_data, position = 'right', legendcol = 2, cornerpos = "right_b", colors_dict) {

  umap_data$Cell_Type_Code <- factor(umap_data$Cell_Type_Code )
  umap_data$Cell_Type2 <- paste0(umap_data$Cell_Type_Code, "_", umap_data$Cell_Type)
  umap_data$Cell_Type2 <- factor(umap_data$Cell_Type2)

  sorted_levels <- umap_data$Cell_Type2 %>%
    levels() %>%
    data.frame(Cell_Type2 = .) %>%
    mutate(Num = as.numeric(str_extract(Cell_Type2, "\\d+"))) %>%
    arrange(Num) %>%
    pull(Cell_Type2)


  umap_data$Cell_Type2 <- factor(umap_data$Cell_Type2, levels = sorted_levels)


library(dplyr)

set.seed(123)


umap_data_center <- umap_data %>%
  group_by(Cell_Type_Code) %>%
  summarise(x_m = median(UMAP1), y_m = median(UMAP2)) %>%
  mutate(Cell_Type_Code = factor(Cell_Type_Code))


threshold <- 0.5


distance_matrix <- as.matrix(dist(umap_data_center[, c('x_m', 'y_m')]))


for (i in 1:(nrow(umap_data_center) - 1)) {
  for (j in (i + 1):nrow(umap_data_center)) {
    if (distance_matrix[i, j] < threshold) {

      direction_x <- ifelse(runif(1) > 0.5, 1, -1)
      direction_y <- ifelse(runif(1) > 0.5, 1, -1)


      displacement_x <- runif(1, min = 0.2, max = 1) * direction_x
      displacement_y <- runif(1, min = 0.2, max = 1) * direction_y


      umap_data_center$x_m[i] <- umap_data_center$x_m[i] + displacement_x
      umap_data_center$x_m[j] <- umap_data_center$x_m[j] - displacement_x

      umap_data_center$y_m[i] <- umap_data_center$y_m[i] + displacement_y
      umap_data_center$y_m[j] <- umap_data_center$y_m[j] - displacement_y
    }
  }
}
  umap_plot <- ggplot() +
    ggrastr::geom_point_rast(data = umap_data, aes(x = UMAP1, y = UMAP2, color = Cell_Type2), size = 0.9, alpha = 0.50, key_glyph = draw_number_circle, raster.dpi = 300) +
    ggrepel::geom_label_repel(
      data = umap_data_center,
      aes(x = x_m, y = y_m, label = Cell_Type_Code),
      size = 7,
      family = "Arial",
      fill = "grey90",
      label.size = NA,
      label.padding = unit(0.28, "lines"),
      label.r = unit(0.32, "cm"),
      box.padding = unit(0.24, "lines"),
      point.padding = unit(0.02, "lines"),
      force = 4,
      force_pull = 0.12,
      max.overlaps = Inf,
      max.time = 10,
      seed = 123,
      segment.color = NA
    ) +
    scale_color_manual(values = colors_dict) +
    guides(color = guide_legend(override.aes = list(label = levels(factor(umap_data$Cell_Type_Code)), size = 7),
                                ncol = legendcol,
                                title.theme = element_text(size = 16),
                                label.theme = element_text(size = 16))) +

    theme_minimal() +
    theme(panel.spacing.y = unit(0, "mm"),
          axis.text = element_text(color = "
                                   black"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          legend.spacing.y = unit(0.65, "cm"),
          legend.key.height = unit(0.95, "cm"),
          legend.key.width = unit(0.9, "cm"),
          legend.key.spacing.y = unit(0.10, "cm"),
          legend.key.spacing.x = unit(0.75, "cm"),
          strip.text = element_text(face = "bold"),
          legend.position = position,
          aspect.ratio = 1
    ) + labs(color = "")

  return(umap_plot)
}


create_umap_plot_grey <- function(umap_data,alldata, position = 'right', legendcol = 2, cornerpos = "right_b", colors_dict) {

  umap_data$Cell_Type_Code <- factor(umap_data$Cell_Type_Code )
  umap_data$Cell_Type2 <- paste0(umap_data$Cell_Type_Code, "_", umap_data$Cell_Type)
  umap_data$Cell_Type2 <- factor(umap_data$Cell_Type2)

  sorted_levels <- umap_data$Cell_Type2 %>%
    levels() %>%
    data.frame(Cell_Type2 = .) %>%
    mutate(Num = as.numeric(str_extract(Cell_Type2, "\\d+"))) %>%
    arrange(Num) %>%
    pull(Cell_Type2)


  umap_data$Cell_Type2 <- factor(umap_data$Cell_Type2, levels = sorted_levels)


set.seed(123)


library(dplyr)


umap_data_center <- umap_data %>%
  group_by(Cell_Type_Code) %>%
  summarise(x_m = median(UMAP1), y_m = median(UMAP2)) %>%
  mutate(Cell_Type_Code = factor(Cell_Type_Code))

filtered_data <- umap_data_center %>%
  filter(Cell_Type_Code %in% c(5,6,19, 20, 21))
print(filtered_data)
threshold <- 1

distance_matrix <- as.matrix(dist(umap_data_center[, c('x_m', 'y_m')]))

for (i in 1:(nrow(umap_data_center) - 1)) {
  for (j in (i + 1):nrow(umap_data_center)) {
    if (distance_matrix[i, j] < threshold) {
      if (umap_data_center$Cell_Type_Code[i] == 5 || umap_data_center$Cell_Type_Code[j] == 5) {
        displacement_x <- 0
        displacement_y <- 0
      } else if (umap_data_center$Cell_Type_Code[i] == 6 || umap_data_center$Cell_Type_Code[j] == 6) {
        displacement_x <- 1
        displacement_y <- -1
      } else if (umap_data_center$Cell_Type_Code[i] == 19 || umap_data_center$Cell_Type_Code[j] == 19) {
        displacement_x <- -0.8
        displacement_y <- 0.4
      } else if (umap_data_center$Cell_Type_Code[i] == 21 || umap_data_center$Cell_Type_Code[j] == 21) {
        displacement_x <- -1
        displacement_y <- 0.1
      } else {
        direction_x <- ifelse(runif(1) > 0.5, 1, -1)
        direction_y <- ifelse(runif(1) > 0.5, 1, -1)
        displacement_x <- runif(1, min = 0.2, max = 1) * direction_x
        displacement_y <- runif(1, min = 0.2, max = 1) * direction_y
      }
      umap_data_center$x_m[i] <- umap_data_center$x_m[i] + displacement_x
      umap_data_center$x_m[j] <- umap_data_center$x_m[j] - displacement_x
      umap_data_center$y_m[i] <- umap_data_center$y_m[i] + displacement_y
      umap_data_center$y_m[j] <- umap_data_center$y_m[j] - displacement_y
    }
  }
}

filtered_data <- umap_data_center %>%
  filter(Cell_Type_Code %in% c(5,6,19, 20, 21))

print(filtered_data)
  umap_plot <- ggplot() + ggrastr::geom_point_rast(data = alldata, aes(x = UMAP1, y = UMAP2), size = 0.5, alpha = 0.8, color = "grey", raster.dpi = 300)+
    ggrastr::geom_point_rast(data = umap_data, aes(x = UMAP1, y = UMAP2, color = Cell_Type2), size = 0.5, alpha = 0.8, key_glyph = draw_number_circle, raster.dpi = 300) +
    geom_label(
      data = umap_data_center,
      aes(x = x_m, y = y_m, label = Cell_Type_Code),
      size = 5.2,
      family = "Arial",
      fill = "grey90",
      label.size = NA,
      label.padding = unit(0.22, "lines"),
      label.r = unit(0.25, "cm")
    ) +
    scale_color_manual(values = colors_dict) +
    guides(color = guide_legend(override.aes = list(label = levels(factor(umap_data$Cell_Type_Code)), size = 5),
                                ncol = legendcol,
                                title.theme = element_text(size = 14),
                                label.theme = element_text(size = 12))) +

    theme_minimal() +
    theme(panel.spacing.y = unit(0, "mm"),
          axis.text = element_text(color = "
                                   black"),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          legend.spacing.y = unit(0.5, "cm"),
          strip.text = element_text(face = "bold"),
          legend.position = position,
          aspect.ratio = 1
    ) + labs(color = "")

  return(umap_plot)
}
