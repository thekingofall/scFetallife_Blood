.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
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
      size = unit(2, "char"),
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
      gp = gpar(col = "black", fontfamily = "Arial", fontsize = 12)
    )
  )
}

repel_circle_centers <- function(center_data, reference_data, min_sep = 0.045, iterations = 800) {
  x_limits <- range(reference_data$UMAP1, na.rm = TRUE)
  y_limits <- range(reference_data$UMAP2, na.rm = TRUE)
  x_span <- diff(x_limits)
  y_span <- diff(y_limits)
  x_norm <- (center_data$x_m - x_limits[1]) / x_span
  y_norm <- (center_data$y_m - y_limits[1]) / y_span
  n <- nrow(center_data)
  if (n > 1) {
    for (iteration in seq_len(iterations)) {
      moved <- FALSE
      for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
          dx <- x_norm[j] - x_norm[i]
          dy <- y_norm[j] - y_norm[i]
          distance <- sqrt(dx^2 + dy^2)
          if (distance < min_sep) {
            if (distance < 1e-8) {
              angle <- ((i * 37 + j * 53) %% 360) * pi / 180
              dx <- cos(angle)
              dy <- sin(angle)
              distance <- 1
            }
            shift <- (min_sep - distance) / 2 + 0.0005
            x_norm[i] <- x_norm[i] - shift * dx / distance
            x_norm[j] <- x_norm[j] + shift * dx / distance
            y_norm[i] <- y_norm[i] - shift * dy / distance
            y_norm[j] <- y_norm[j] + shift * dy / distance
            moved <- TRUE
          }
        }
      }
      if (!moved) break
    }
  }
  center_data$x_m <- x_limits[1] + x_norm * x_span
  center_data$y_m <- y_limits[1] + y_norm * y_span
  center_data
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
umap_data_center <- repel_circle_centers(umap_data_center, umap_data)
  umap_plot <- ggplot() +
    ggrastr::geom_point_rast(data = umap_data, aes(x = UMAP1, y = UMAP2, color = Cell_Type2), size = 0.5, alpha = 0.8, key_glyph = draw_number_circle, raster.dpi = 300) +
    geom_point(
      data = umap_data_center,
      aes(x = x_m, y = y_m),
      shape = 21,
      size = 6.2,
      stroke = 0.45,
      color = "#303030",
      fill = "#F2F2F2"
    ) +
    geom_text(
      data = umap_data_center,
      aes(x = x_m, y = y_m, label = Cell_Type_Code),
      size = 3.8,
      family = "Arial",
      color = "black"
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

umap_data_center <- repel_circle_centers(umap_data_center, alldata)

filtered_data <- umap_data_center %>%
  filter(Cell_Type_Code %in% c(5,6,19, 20, 21))

print(filtered_data)
  umap_plot <- ggplot() + ggrastr::geom_point_rast(data = alldata, aes(x = UMAP1, y = UMAP2), size = 0.5, alpha = 0.8, color = "grey", raster.dpi = 300)+
    ggrastr::geom_point_rast(data = umap_data, aes(x = UMAP1, y = UMAP2, color = Cell_Type2), size = 0.5, alpha = 0.8, key_glyph = draw_number_circle, raster.dpi = 300) +
    geom_point(
      data = umap_data_center,
      aes(x = x_m, y = y_m),
      shape = 21,
      size = 6.2,
      stroke = 0.45,
      color = "#303030",
      fill = "#F2F2F2"
    ) +
    geom_text(
      data = umap_data_center,
      aes(x = x_m, y = y_m, label = Cell_Type_Code),
      size = 3.8,
      family = "Arial",
      color = "black"
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
