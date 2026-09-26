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

#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

panel_dir <- normalizePath(file.path(.scf_start_dir, ".."), winslash="/")
input_path <- file.path(.scf_shared, "Main_Figure_07/a_M07a/03_PlotData/M07a_Plasma_Olink_spearman_plot_data.csv")
output_dir <- file.path(panel_dir, "02_Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(input_path, check.names = FALSE)
stopifnot(all(c("feature", "r", "p") %in% names(df)))
stopifnot(all(is.finite(df$r)), all(is.finite(df$p)))

df <- df[order(df$r, decreasing = TRUE), , drop = FALSE]
df$feature <- factor(df$feature, levels = df$feature)

cutpoints <- c(0.001, 0.005, 0.01, 0.05, 0.1)
cutlabels <- c("<=0.001", "0.005", "0.01", "0.05", ">=0.1")
df$pvalue2 <- pmin(pmax(df$p, min(cutpoints)), max(cutpoints))
df$group <- ifelse(
  df$pvalue2 <= 0.005,
  1,
  ifelse(df$pvalue2 <= 0.01, 2, ifelse(df$pvalue2 <= 0.05, 3, 4))
)
df$pvalue3 <- mapply(
  function(value, group) {
    scales::rescale(
      value,
      to = c(group - 1, group),
      from = c(cutpoints[group], cutpoints[group + 1])
    )
  },
  df$pvalue2,
  df$group
)

reference_palette <- c(
  "#C15040",
  "#D6884A",
  "#D8AF64",
  "#7CB6C6",
  "#6873B4"
)

plot <- ggplot(df, aes(x = feature, y = r, fill = pvalue3, label = feature)) +
  geom_col(width = 0.9) +
  scale_fill_gradientn(
    colours = reference_palette,
    limits = c(0, 4),
    breaks = 0:4,
    labels = cutlabels
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +
  geom_text(
    aes(y = ifelse(r > 0, -0.03, 0.03), angle = 90),
    hjust = ifelse(df$r > 0, 1, 0),
    vjust = 0.5,
    size = 4,
    fontface = "bold",
    family = "Arial",
    color = "black"
  ) +
  geom_text(
    aes(
      y = ifelse(r > 0, r + 0.02, r - 0.04),
      label = ifelse(p < 0.05, "*", "")
    ),
    hjust = 0.5,
    size = 4,
    family = "Arial",
    color = "black"
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5),
    limits = c(-1.08, 1.08),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Spearman r", fill = "p-value") +
  theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.8),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.height = grid::unit(0.9, "cm"),
    plot.margin = margin(8, 8, 8, 8)
  )

pdf_path <- file.path(output_dir, "M07a_Plasma_Olink_spearman.pdf")
png_path <- file.path(output_dir, "M07a_Plasma_Olink_spearman.png")

cairo_pdf(pdf_path, width = 15, height = 5.1, family = "Arial")
print(plot)
dev.off()

png(
  png_path,
  width = 15,
  height = 5.1,
  units = "in",
  res = 300,
  type = "cairo",
  family = "Arial",
  bg = "white"
)
print(plot)
dev.off()

write.csv(
  transform(
    df,
    pvalue_mapped = pvalue2,
    pvalue_colour_coordinate = pvalue3
  ),
  file.path(panel_dir, "03_PlotData", "M07a_reference_palette_plot_data.csv"),
  row.names = FALSE
)

message("Saved: ", pdf_path)
message("Saved: ", png_path)
