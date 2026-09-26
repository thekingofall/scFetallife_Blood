.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2", "readr"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

clean_root <- .scf_project_root
input_dir <- file.path(.scf_shared, "Main_Figure_04_05/03_VJ_DJ_Pairing")
output_root <- file.path(clean_root, "Output")

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

is_true <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

natural_gene_order <- function(x) {
  x[order(gsub("[0-9]+", "", x),
          suppressWarnings(as.numeric(sub(".*?([0-9]+).*", "\\1", x))),
          x,
          na.last = TRUE)]
}

read_order_column <- function(path, sep = "\t") {
  tab <- read.table(path, header = TRUE, sep = sep, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "")
  as.character(tab[[1]])
}

render_original_vj <- function(input_name, out_dir, stem, plot_title,
                               x_name, x_label, width, height,
                               x_order = NULL, y_order = NULL,
                               text_sizes = c("0-0.5" = 4, "0.6-0.7" = 5,
                                              "0.7-0.8" = 6, "0.8-0.9" = 7,
                                              "0.8-1.0" = 8),
                               fill_missing_zero = FALSE) {
  dat <- read_csv(file.path(input_dir, input_name), show_col_types = FALSE)
  
  dat$Per <- as.numeric(dat$percentage)

  if (is.null(x_order)) x_order <- natural_gene_order(unique(dat[[x_name]]))
  if (is.null(y_order)) y_order <- natural_gene_order(unique(dat$Jgene))
  x_order <- x_order[x_order %in% unique(dat[[x_name]])]
  y_order <- y_order[y_order %in% unique(dat$Jgene)]
  dat[[x_name]] <- factor(dat[[x_name]], levels = x_order)
  dat$Jgene <- factor(dat$Jgene, levels = rev(y_order))

  plot_dat <- dat
  if (fill_missing_zero) {
    plot_dat <- expand.grid(
      x_gene = x_order,
      Jgene = y_order,
      stringsAsFactors = FALSE
    )
    names(plot_dat)[1] <- x_name
    dat_join <- dat
    dat_join[[x_name]] <- as.character(dat_join[[x_name]])
    dat_join$Jgene <- as.character(dat_join$Jgene)
    plot_dat <- left_join(plot_dat, dat_join, by = c(x_name, "Jgene"))
    missing_combo <- is.na(plot_dat$Per)
    plot_dat$Per[missing_combo] <- 0
    if ("percentage" %in% names(plot_dat)) {
      plot_dat$percentage[missing_combo] <- 0
    }
    if ("frequency" %in% names(plot_dat)) {
      plot_dat$frequency[missing_combo] <- 0
    }
    if ("legacy_display" %in% names(plot_dat)) {
      plot_dat$legacy_display[missing_combo] <- FALSE
    }
    if ("component" %in% names(plot_dat)) {
      component_value <- unique(na.omit(as.character(dat$component)))
      if (length(component_value) == 1) {
        plot_dat$component[missing_combo] <- component_value
      }
    }
    if ("combo" %in% names(plot_dat)) {
      plot_dat$combo[missing_combo] <- paste0(
        plot_dat[[x_name]][missing_combo], "_", plot_dat$Jgene[missing_combo]
      )
    }
    plot_dat[[x_name]] <- factor(plot_dat[[x_name]], levels = x_order)
    plot_dat$Jgene <- factor(plot_dat$Jgene, levels = rev(y_order))
  }

  sig <- dat[is_true(dat$legacy_display), , drop = FALSE]
  sig$cor <- as.numeric(sig$rho)
  sig$VDJcor <- as.character(sig$direction)
  sig$size_factor <- cut(
    abs(sig$cor),
    breaks = c(-Inf, 0.5, 0.6, 0.7, 0.8, Inf),
    labels = c("0-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-1.0")
  )

  p <- ggplot() +
    geom_tile(
      data = plot_dat,
      aes(x = .data[[x_name]], y = Jgene, fill = Per),
      color = NA,
      linewidth = 0
    ) +
    scale_fill_gradientn(
      colors = rev(c("#C71000FF", "#edae11", "white", "#246BAE")),
      limits = c(0, max(plot_dat$Per, na.rm = TRUE)),
      na.value = "#246BAE"
    ) +
    theme_minimal() +
    labs(
      title = plot_title,
      x = x_label,
      y = "J gene",
      fill = "Percentage(%)"
    ) +
    theme(
      text = element_text(family = "Arial"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    geom_text(
      data = sig,
      aes(
        x = .data[[x_name]], y = Jgene, label = VDJcor,
        color = VDJcor, size = size_factor
      )
    ) +
    scale_color_manual(values = c("#3D3B25FF", "#8F1336")) +
    scale_size_manual(values = text_sizes)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p,
         width = width, height = height, dpi = 300)
  table_dir <- file.path(.scf_shared, basename(dirname(dirname(out_dir))),
                         basename(dirname(out_dir)), "03_PlotData")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(plot_dat, file.path(table_dir, paste0(stem, "_plot_data.csv")))
  write_csv(sig, file.path(table_dir, paste0(stem, "_significant_labels.csv")))
  invisible(p)
}

tcr_gene_root <- file.path(.scf_shared, "00_Common/03_Receptor_Gene_Order", "TCR")
tra_locus <- read_csv(file.path(tcr_gene_root, "TRA_Gene_Locus.csv"), show_col_types = FALSE)
trb_locus <- read_csv(file.path(tcr_gene_root, "TRB_Gene_Locus.csv"), show_col_types = FALSE)
tra_v_order <- tra_locus$IMGT_gene_name[grepl("^TRAV", tra_locus$IMGT_gene_name)]
tra_j_order <- tra_locus$IMGT_gene_name[grepl("^TRAJ", tra_locus$IMGT_gene_name)]
trb_v_order <- trb_locus$IMGT_gene_name[grepl("^TRBV", trb_locus$IMGT_gene_name)]
trb_j_order <- trb_locus$IMGT_gene_name[grepl("^TRBJ", trb_locus$IMGT_gene_name)]

m04f_dir <- file.path(.scf_project_root, "Main_Figure_04", "f_M04f", "02_Figures")
if (!identical(Sys.getenv("M05F_ONLY"), "1")) {
  render_original_vj(
    "M04f_TRA_VJ_Frequencies_and_Age_Associations.csv", m04f_dir,
    "M04f_TRA_all_VJcombined", "PBMC All TRA",
    "Vgene", "V gene", 7, 7, tra_v_order, tra_j_order,
    c("0-0.5" = 5, "0.6-0.7" = 5, "0.7-0.8" = 6,
      "0.8-0.9" = 7, "0.8-1.0" = 8),
    fill_missing_zero = TRUE
  )
  render_original_vj(
    "M04f_TRB_VJ_Frequencies_and_Age_Associations.csv", m04f_dir,
    "M04f_TRB_all_VJcombined", "PBMC All TRB",
    "Vgene", "V gene", 8, 4, trb_v_order, trb_j_order,
    fill_missing_zero = TRUE
  )
}

if (identical(Sys.getenv("M04F_ONLY"), "1")) {
  message("M04f complete")
  quit(save = "no", status = 0)
}

bcr_gene_root <- file.path(.scf_shared, "00_Common/03_Receptor_Gene_Order", "BCR")
igh_v_order <- read_order_column(file.path(bcr_gene_root, "IGHV.txt"))
igk_v_order <- read_order_column(file.path(bcr_gene_root, "IGKV.txt"))
igl_v_order <- read_order_column(file.path(bcr_gene_root, "IGLV.txt"))
igk_j_order <- read_order_column(file.path(bcr_gene_root, "IGKJ.txt"))
igl_j_order <- read_order_column(file.path(bcr_gene_root, "IGLJ.txt"))
igh_j_order <- paste0("IGHJ", 1:6)

m05f_dir <- file.path(.scf_project_root, "Main_Figure_05", "f_M05f", "02_Figures")
render_original_vj(
  "M05f_IGH_VJ_Frequencies_and_Age_Associations.csv", m05f_dir,
  "M05f_BCRH_all_VJcombined", "PBMC All BCRH",
  "Vgene", "V gene", 8, 4, igh_v_order, igh_j_order,
  fill_missing_zero = TRUE
)
render_original_vj(
  "M05f_IGK_VJ_Frequencies_and_Age_Associations.csv", m05f_dir,
  "M05f_BCRK_all_VJcombined", "PBMC All BCRK",
  "Vgene", "V gene", 8, 4, igk_v_order, igk_j_order,
  fill_missing_zero = TRUE
)
render_original_vj(
  "M05f_IGL_VJ_Frequencies_and_Age_Associations.csv", m05f_dir,
  "M05f_BCRL_all_VJcombined", "PBMC All BCRL",
  "Vgene", "V gene", 8, 4, igl_v_order, igl_j_order,
  fill_missing_zero = TRUE
)

if (identical(Sys.getenv("M05F_ONLY"), "1")) {
  message("M05f complete")
  quit(save = "no", status = 0)
}

igh_d_order <- read_order_column(file.path(bcr_gene_root, "IGHD.csv"), sep = ",")
message("M04f and M05f complete")
