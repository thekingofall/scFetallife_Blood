.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "grid", "patchwork", "ragg", "readr", "tidyr"))

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
  library(patchwork)
  library(readr)
  library(tidyr)
})

repo <- .scf_project_root
input_dir <- file.path(.scf_start_dir, "../03_PlotData")

save_pair <- function(plot, directory, stem, width, height, dpi = 320) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, limitsize = FALSE)
  ggsave(file.path(directory, paste0(stem, ".png")), plot, width = width, height = height, device = ragg::agg_png, dpi = dpi, bg = "white", limitsize = FALSE)
}


result1 <- read_csv(file.path(input_dir, "M08c_factor_age_spearman_statistics.csv"), show_col_types = FALSE)
create_qplot <- function(df) {
  df$qvalue <- cut(
    df$q,
    breaks = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
    labels = c("< 0.0001", "0.0001-0.001", "0.001-0.01", "0.01-0.05", "> 0.05"),
    include.lowest = TRUE
  )
  colors <- c(
    "< 0.0001" = "#F6313E", "0.0001-0.001" = "#FF7F50",
    "0.001-0.01" = "#f6c619", "0.01-0.05" = "#999933", "> 0.05" = "#00AF99"
  )
  ggplot(df, aes(x = reorder(factor, rho), y = rho, fill = qvalue, label = sub("Factor", "Factor ", factor))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    theme_classic(base_family = "Arial") +
    theme(
      legend.position = "inside", legend.position.inside = c(0.15, 0.86),
      legend.text = element_text(size = 9), legend.title = element_text(size = 9),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.line.x = element_blank(), axis.line.y = element_line(color = "black", linewidth = 0.5)
    ) +
    labs(x = "Gene", y = "Spearman's R correlation with pcw", fill = "q value") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_text(
      aes(y = ifelse(rho > 0, -0.01, 0.01), angle = 90),
      hjust = ifelse(df$rho > 0, 1, 0), vjust = 0.5, family = "Arial"
    ) +
    xlab("") +
    scale_y_continuous(breaks = seq(-1, 1, 0.2), labels = seq(-1, 1, by = 0.2))
}

p8c <- create_qplot(result1)
out8c <- file.path(.scf_project_root, "Main_Figure_08", "c_M08c", "02_Figures")
save_pair(p8c, out8c, "M08c_factor_age_spearman", 5.2, 3.4)


