.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "grid", "scales"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

panel_dir <- normalizePath(file.path(.scf_start_dir,".."),winslash="/")
data_dir <- file.path(panel_dir, "03_PlotData")
figure_dir <- file.path(panel_dir, "02_Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.delim(
  file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_16/a_S16a/03_PlotData/S16_spearman_results.tsv"),
  check.names = FALSE
)
stopifnot(all(c("assay", "spearman_r", "p_value") %in% names(df)))
stopifnot(all(is.finite(df$spearman_r)), all(is.finite(df$p_value)))

df <- df[order(df$spearman_r, decreasing = TRUE), , drop = FALSE]
df$assay <- factor(df$assay, levels = df$assay)


cutpoints <- c(0.001, 0.005, 0.01, 0.05, 0.1)
cutlabels <- c("<=0.001", "0.005", "0.01", "0.05", ">=0.1")
df$pvalue2 <- pmin(pmax(df$p_value, min(cutpoints)), max(cutpoints))
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

wangflow_palette <- c(
  "#C15040",
  "#D6884A",
  "#D8AF64",
  "#7CB6C6",
  "#6873B4"
)

plot <- ggplot(df, aes(x = assay, y = spearman_r, fill = pvalue3, label = assay)) +
  geom_col(width = 0.9) +
  scale_fill_gradientn(
    colours = wangflow_palette,
    limits = c(0, 4),
    breaks = 0:4,
    labels = cutlabels
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.75) +
  geom_text(
    aes(y = ifelse(spearman_r > 0, -0.025, 0.025), angle = 90),
    hjust = ifelse(df$spearman_r > 0, 1, 0),
    vjust = 0.5,
    size = 4,
    fontface = "bold",
    family = "Arial",
    color = "black"
  ) +
  geom_text(
    aes(
      y = ifelse(spearman_r > 0, spearman_r + 0.025, spearman_r - 0.04),
      label = ifelse(p_value < 0.05, "*", "")
    ),
    hjust = 0.5,
    size = 4,
    family = "Arial",
    color = "black"
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5),
    limits = c(-1.05, 1.05),
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
    plot.margin = margin(6, 6, 6, 6)
  )

pdf_path <- file.path(figure_dir, "FigureS16a_Bulk_mRNA_Associations.pdf")
png_path <- file.path(figure_dir, "FigureS16a_Bulk_mRNA_Associations.png")

cairo_pdf(pdf_path, width = 12, height = 4.2, family = "Arial")
print(plot)
dev.off()

png(
  png_path,
  width = 12,
  height = 4.2,
  units = "in",
  res = 300,
  type = "cairo",
  family = "Arial",
  bg = "white"
)
print(plot)
dev.off()

write.table(
  transform(
    df,
    pvalue_mapped = pvalue2,
    pvalue_colour_coordinate = pvalue3
  ),
  file.path(data_dir, "S16a_Wangflow_gradient_plot_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

