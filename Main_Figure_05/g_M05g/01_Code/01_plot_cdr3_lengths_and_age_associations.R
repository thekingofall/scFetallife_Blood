.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("gridExtra", "tidyverse"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

clean_root <- .scf_project_root
input_dir <- file.path(.scf_shared, "Main_Figure_04_05/04_CDR3_Length")
output_root <- file.path(clean_root, "Output")

suppressPackageStartupMessages({
  library(tidyverse)
  library(gridExtra)
})

original_distribution_plot <- function(df, plot_title, smooth_method = NULL,
                                       smooth_size = 0.3) {
  df <- df %>%
    transmute(
      Length = factor(as.numeric(length), levels = sort(unique(as.numeric(length)))),
      Sample = MainID,
      ratio = as.numeric(percentage),
      Post_Conception_Age_Weeks = as.numeric(sub("^[A-Za-z]+([0-9]+([.][0-9]+)?).*", "\\1", MainID))
    )
  

  p <- ggplot(df, aes(Length, ratio, color = Post_Conception_Age_Weeks, group = Sample)) +
    geom_point(shape = 1, size = 1)
  if (is.null(smooth_method)) {
    p <- p + geom_smooth(se = FALSE, size = smooth_size)
  } else {
    p <- p + geom_smooth(se = FALSE, method = smooth_method, size = smooth_size)
  }
  p +
    ggtitle(plot_title) +
    theme_linedraw() +
    xlab("Length") +
    ylab("ratio(%)") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      panel.border = element_rect(linetype = "solid", colour = "black", size = 1.5)
    ) +
    scale_color_gradientn(
      colours = rev(colorRampPalette(
        c("#C71000B2", "#FF6F00B2", "#6a73cf", "#00AF99")
      )(100))
    )
}

original_pvalue_plot <- function(df, plot_title, y_step = 0.25) {
  stat <- df %>%
    transmute(
      lengths = as.numeric(length),
      R_value = as.numeric(rho),
      p_value = as.numeric(p_value)
    ) %>%
    distinct()
  stat$pvalue_cat <- cut(
    stat$p_value,
    breaks = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
    labels = c("< 0.0001", "0.0001-0.001", "0.001-0.01", "0.01-0.05", "> 0.05"),
    include.lowest = TRUE
  )
  colors <- c(
    "< 0.0001" = "#C71000B2",
    "0.0001-0.001" = "#FF6348B2",
    "0.001-0.01" = "#FF95A8B2",
    "0.01-0.05" = "#8A4198B2",
    "> 0.05" = "#008EA0B2"
  )
  ggplot(stat, aes(x = factor(lengths), y = R_value, fill = pvalue_cat)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    labs(x = "Lengths", y = "R Value(spearman) ", fill = "P Value") +
    theme_bw() +
    ggtitle(plot_title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      panel.border = element_rect(linetype = "solid", colour = "black", size = 1.5)
    ) +
    scale_y_continuous(
      breaks = seq(-1, 1, by = y_step),
      labels = seq(-1, 1, by = y_step)
    )
}

save_panel_pair <- function(dat, out_dir, stem, dist_title, stat_title,
                            smooth_method = NULL, smooth_size = 0.3,
                            dist_height = 4, y_step = 0.25) {
  
  p_dist <- original_distribution_plot(dat, dist_title, smooth_method, smooth_size)
  p_stat <- original_pvalue_plot(dat, stat_title, y_step)
  p_combined <- arrangeGrob(p_dist, p_stat, ncol = 1, heights = c(4, 3))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(stem, "_distribution.pdf")), p_dist,
         width = 8, height = dist_height)
  ggsave(file.path(out_dir, paste0(stem, "_distribution.png")), p_dist,
         width = 8, height = dist_height, dpi = 300)
  ggsave(file.path(out_dir, paste0(stem, "_age_correlation.pdf")), p_stat,
         width = 8, height = 3)
  ggsave(file.path(out_dir, paste0(stem, "_age_correlation.png")), p_stat,
         width = 8, height = 3, dpi = 300)
  ggsave(file.path(out_dir, paste0(stem, "_combined.pdf")), p_combined,
         width = 8, height = 7)
  ggsave(file.path(out_dir, paste0(stem, "_combined.png")), p_combined,
         width = 8, height = 7, dpi = 300)
  table_dir <- file.path(.scf_shared, basename(dirname(dirname(out_dir))),
                         basename(dirname(out_dir)), "03_PlotData")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(dat, file.path(table_dir, paste0(stem, "_plot_data.csv")))
}

m04g_dir <- file.path(.scf_project_root, "Main_Figure_04", "g_M04g", "02_Figures")
tra <- read_csv(file.path(input_dir, "M04g_TRA_CDR3_Length_and_Age_Associations.csv"),
                show_col_types = FALSE)
trb <- read_csv(file.path(input_dir, "M04g_TRB_CDR3_Length_and_Age_Associations.csv"),
                show_col_types = FALSE)
save_panel_pair(
  tra, m04g_dir, "M04g_TRA_CDR3",
  "TRA CDR3(aa) in PBMC", "Percentage of TRA Lengths Over Time",
  smooth_method = "gam", smooth_size = 0.3, dist_height = 4, y_step = 0.25
)
save_panel_pair(
  trb, m04g_dir, "M04g_TRB_CDR3",
  "TRB CDR3(aa) in PBMC", "Percentage of TRB Lengths Over Time",
  smooth_method = "gam", smooth_size = 0.3, dist_height = 4, y_step = 0.25
)

m05g_dir <- file.path(.scf_project_root, "Main_Figure_05", "g_M05g", "02_Figures")
igh <- read_csv(file.path(input_dir, "M05g_IGH_CDR3_Length_and_Age_Associations.csv"),
                show_col_types = FALSE)
iglk <- read_csv(file.path(input_dir, "M05g_IGLK_CDR3_Length_and_Age_Associations.csv"),
                 show_col_types = FALSE)
save_panel_pair(
  igh, m05g_dir, "M05g_BCRH_CDR3",
  "BCRH CDR3(aa) in PBMC", "Percentage of BCRH Lengths Over Time",
  smooth_method = NULL, smooth_size = 0.3, dist_height = 5, y_step = 0.25
)
save_panel_pair(
  iglk, m05g_dir, "M05g_BCRLK_CDR3",
  "BCRL/K CDR3(aa) in PBMC", "Percentage of BCRK/L Lengths Over Time",
  smooth_method = NULL, smooth_size = 0.5, dist_height = 5, y_step = 0.2
)

message("M04g and M05g complete")
