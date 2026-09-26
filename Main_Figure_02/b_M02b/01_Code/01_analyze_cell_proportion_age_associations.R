.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "grid", "ragg", "readr"))

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
})

repo <- .scf_project_root
input_pbmc_file <- file.path(.scf_shared, "Main_Figure_02/ab_M02ab/M02b_PBMC_CD45_Cell_Proportions.csv")
input_ery_file <- file.path(.scf_shared, "Main_Figure_02/ab_M02ab/M02b_Erythroid_Cell_Proportions.csv")
output_dir <- file.path(.scf_project_root, "Main_Figure_02", "b_M02b", "02_Figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plot_data_dir <- file.path(.scf_project_root, "Main_Figure_02", "b_M02b", "03_PlotData")
dir.create(plot_data_dir, recursive = TRUE, showWarnings = FALSE)

numeric_order <- function(values) {
  values <- unique(as.character(values))
  values[order(as.numeric(sub("_.*$", "", values)), values)]
}

correlation_table <- function(frame, cohort, method) {
  rows <- lapply(numeric_order(frame$Last_cell_type_num), function(code) {
    part <- frame[frame$Last_cell_type_num == code & is.finite(frame$Post_Conception_Age_Weeks) & is.finite(frame$percentage), ]
    status <- "ok"
    estimate <- NA_real_
    p_value <- NA_real_
    if (nrow(part) < 3) {
      status <- "insufficient"
    } else if (length(unique(part$Post_Conception_Age_Weeks)) < 2 || length(unique(part$percentage)) < 2) {
      status <- "constant"
    } else {
      test <- suppressWarnings(cor.test(part$Post_Conception_Age_Weeks, part$percentage, method = method))
      estimate <- unname(test$estimate)
      p_value <- test$p.value
      if (!is.finite(estimate) || !is.finite(p_value)) status <- "nonfinite"
    }
    significant <- status == "ok" && p_value < 0.05
    data.frame(
      cohort = cohort, organ = cohort, Last_cell_type_num = code,
      Last_cell_type = as.character(part$Last_cell_type[1]), method = method,
      sample_count = nrow(part), estimate = estimate, p_value = p_value,
      threshold = 0.05, status = status, significant = significant,
      plotted = significant && as.character(part$Last_cell_type[1]) != "Others",
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

pbmc_input <- read_csv(input_pbmc_file, show_col_types = FALSE)
ery_input <- read_csv(input_ery_file, show_col_types = FALSE)
input <- bind_rows(pbmc_input, ery_input)

stats <- correlation_table(input, "PBMC", "spearman")
valid_tests <- stats$status == "ok" & is.finite(stats$p_value)
stats$q_value <- NA_real_
stats$q_value[valid_tests] <- p.adjust(stats$p_value[valid_tests], method = "BH")
stats$significant_q <- valid_tests & stats$q_value < 0.05
selected <- stats[stats$status == "ok" & stats$p_value < 0.05, ]
manuscript_types <- c(
  "1_HSC_MPP", "2_MEP", "3_MEMP", "4_Pro-B", "7_CXCR5- Naïve B", "12_Treg",
  "20_Gamma Delta V1 T", "21_Gamma Delta V2 T", "28_Classical Monocytes", "37_Late_ERY"
)
plotted_stats <- stats[match(manuscript_types, stats$Last_cell_type_num), ]
if (anyNA(plotted_stats$Last_cell_type_num)) stop("A manuscript Figure 2b cell type is missing from the complete PBMC composition table")

type_levels <- manuscript_types
plot_data <- input[input$Last_cell_type_num %in% type_levels, ]
plot_data$Last_cell_type_num <- factor(plot_data$Last_cell_type_num, levels = type_levels)
display_labels <- c(
  "HSC/MPP", "MEP", "MEMP", "Pro-B", "CXCR5- Naïve B", "Treg",
  "Vδ1", "Vδ2", "Classical Monocytes", "Late ERY"
)
names(display_labels) <- type_levels
plot_data$DisplayType <- factor(display_labels[as.character(plot_data$Last_cell_type_num)], levels = display_labels)
plotted_stats$DisplayType <- factor(display_labels[plotted_stats$Last_cell_type_num], levels = display_labels)
plot_data$ColorKey <- as.character(plot_data$Last_cell_type_num)

colors <- c(
  "1_HSC_MPP" = "#FF6F00FF", "2_MEP" = "#008EA0FF", "3_MEMP" = "#84D7E1FF",
  "4_Pro-B" = "#FF95A8FF", "7_CXCR5- Naïve B" = "#3D3B25FF", "12_Treg" = "#C71000FF",
  "20_Gamma Delta V1 T" = "#8A4198FF", "21_Gamma Delta V2 T" = "#5A9599FF",
  "28_Classical Monocytes" = "#FF6348FF", "37_Late_ERY" = "#999933"
)
missing_colors <- setdiff(unique(plot_data$ColorKey), names(colors))
if (length(missing_colors)) stop("Missing cell-type colors: ", paste(missing_colors, collapse = ", "))

format_stat <- function(value) {
  ifelse(value < 0.001, formatC(value, format = "e", digits = 2), formatC(value, format = "f", digits = 3))
}
plotted_stats$stat_label <- paste0(
  "r = ", round(plotted_stats$estimate, 2),
  "\np = ", format_stat(plotted_stats$p_value),
  "\nq = ", format_stat(plotted_stats$q_value)
)

plot <- ggplot(plot_data, aes(x = Post_Conception_Age_Weeks, y = percentage)) +
  geom_point(aes(color = ColorKey), size = 2.2) +
  geom_smooth(
    aes(color = ColorKey), method = "lm", se = TRUE,
    show.legend = FALSE, fill = "#aaddcc"
  ) +
  facet_wrap(~ DisplayType, scales = "free", labeller = label_value, ncol = 5) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
    strip.text.x = element_text(size = 9), strip.text.y = element_text(size = 9),
    aspect.ratio = 0.8, panel.spacing = unit(0.65, "lines")
  ) +
  scale_x_continuous(limits = c(10, 40), breaks = seq(10, 40, 5), labels = seq(10, 40, 5)) +
  scale_color_manual(values = colors, drop = FALSE) +
  geom_text(
    data = plotted_stats, inherit.aes = FALSE,
    aes(
      x = Inf, y = Inf,
      label = stat_label
    ),
    hjust = 1.05, vjust = 1.1, size = 3, lineheight = 0.9
  ) +
  xlab("PCW") + ylab("Composition(%)")

ggsave(file.path(output_dir, "M02b_PBMC_age_correlations.pdf"), plot, width = 12, height = 6.5, device = cairo_pdf)
ggsave(file.path(output_dir, "M02b_PBMC_age_correlations.png"), plot, width = 12, height = 6.5, device = ragg::agg_png, dpi = 320, bg = "white")
write_csv(plot_data, file.path(plot_data_dir, "M02b_PBMC_age_correlations_plot_data.csv"))
write_csv(stats, file.path(.scf_shared, "Main_Figure_02/b_M02b/03_PlotData/M02b_PBMC_age_correlations_statistics.csv"))
write_csv(selected, file.path(.scf_shared, "Main_Figure_02/b_M02b/03_PlotData/M02b_PBMC_age_correlations_significant.csv"))
cat("M02b complete; samples=", length(unique(input$MainID)), "; manuscript panels=", nrow(plotted_stats), "\n", sep = "")
