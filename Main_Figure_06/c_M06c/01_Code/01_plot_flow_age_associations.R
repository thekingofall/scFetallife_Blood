.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "grid", "ragg", "readr", "scales"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(readr)
  library(scales)
})

repo <- .scf_project_root
input_dir <- file.path(.scf_shared, "Main_Figure_06/bc_M06bc")
panel_root <- normalizePath(file.path(.scf_start_dir, ".."), winslash="/")
plot_data <- file.path(.scf_shared, "Main_Figure_06/c_M06c/03_PlotData")
dir.create(plot_data,recursive=TRUE,showWarnings=FALSE)

save_plot_pair <- function(plot, directory, stem, width, height, dpi = 320) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(directory, paste0(stem, ".png")), plot, width = width, height = height, device = ragg::agg_png, dpi = dpi, bg = "white")
}


plot_data_bar <- function(df, title = NULL) {
  cutpoints <- c(0.001, 0.005, 0.01, 0.05, 0.1)
  cutlabels <- c("≤0.001", "0.005", "0.01", "0.05", "≥0.1")
  df <- df %>%
    mutate(
      pvalue2 = case_when(
        p >= max(cutpoints) ~ max(cutpoints),
        p < min(cutpoints) ~ min(cutpoints),
        TRUE ~ p
      ),
      group = case_when(
        pvalue2 <= 0.005 ~ 1,
        pvalue2 > 0.005 & pvalue2 <= 0.01 ~ 2,
        pvalue2 > 0.01 & pvalue2 <= 0.05 ~ 3,
        TRUE ~ 4
      )
    ) %>%
    group_by(group) %>%
    mutate(
      pvalue3 = scales::rescale(
        pvalue2,
        to = c(group[1] - 1, group[1]),
        from = c(cutpoints[group[1]], cutpoints[group[1] + 1])
      )
    ) %>%
    ungroup()

  ggplot(df, aes(x = feature, y = r, fill = pvalue3, label = feature)) +
    geom_bar(stat = "identity") +
    scale_fill_gradientn(
      colours = c("#e22b2b", "#ed7c24", "#f8bd19", "#71c8dc", "#6a73cf"),
      limits = c(0, 4), breaks = 0:4, labels = cutlabels
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      legend.key.height = unit(1, "cm"),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(x = "", y = "Spearman r", fill = "p-value", title = title) +
    scale_y_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_text(
      aes(y = ifelse(r > 0, -0.03, 0.03), angle = 90),
      hjust = ifelse(df$r > 0, 1, 0), vjust = 0.5,
      size = 4, fontface = "bold", color = "black"
    ) +
    geom_text(
      aes(y = ifelse(r > 0, r + 0.02, r - 0.04), label = ifelse(p < 0.05, "*", "")),
      hjust = 0.5
    )
}

render_association_bar <- function(stats_path, output_dir, stem, title, feature_col, rho_col, p_col) {
  x <- read_tsv(stats_path, show_col_types = FALSE) %>%
    transmute(
      feature = .data[[feature_col]],
      r = as.numeric(.data[[rho_col]]),
      p = as.numeric(.data[[p_col]])
    ) %>%
    filter(is.finite(r), is.finite(p)) %>%
    arrange(desc(r))
  x$feature <- factor(x$feature, levels = x$feature)
  p <- plot_data_bar(x, title)
  save_plot_pair(p, output_dir, stem, 15, 5.5)
  write_csv(x, file.path(plot_data, paste0(stem, "_plot_data.csv")))
}

out6c <- file.path(.scf_project_root, "Main_Figure_06", "c_M06c", "02_Figures")
render_association_bar(
  file.path(input_dir, "M06c_Flow_Age_Association_Statistics.tsv"),
  out6c, "M06c_SFC_44type_spearman", NULL,
  "feature_display", "rho", "p_value"
)

