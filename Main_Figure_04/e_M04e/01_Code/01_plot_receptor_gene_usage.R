.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ComplexHeatmap", "RColorBrewer", "scales", "tidyverse"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

clean_root <- .scf_project_root
input_dir <- file.path(.scf_shared, "Main_Figure_04_05/02_Gene_Usage")
output_root <- file.path(clean_root, "Output")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(RColorBrewer)
  library(ComplexHeatmap)
})


prepare_long <- function(tab, segment_name, cell_type_mode = c("keep", "merge_b")) {
  cell_type_mode <- match.arg(cell_type_mode)
  x <- tab %>% filter(segment == segment_name)
  if (cell_type_mode == "merge_b") {
    counts <- x %>%
      group_by(MainID, Main_Organ, gene, segment) %>%
      summarise(count = sum(as.numeric(count)), .groups = "drop")
    denoms <- x %>%
      distinct(Last_cell_type, MainID, Main_Organ, segment, segment_total) %>%
      group_by(MainID, Main_Organ, segment) %>%
      summarise(segment_total = sum(as.numeric(segment_total)), .groups = "drop")
    x <- left_join(counts, denoms, by = c("MainID", "Main_Organ", "segment")) %>%
      mutate(
        percentage = ifelse(segment_total > 0, 100 * count / segment_total, 0),
        Last_cell_type = "Naïve B"
      )
  }
  x %>%
    transmute(
      Sample = paste0(Last_cell_type, "_", gsub("_", ".", MainID)),
      variable = gene,
      value = as.numeric(percentage)
    )
}

TCRgene_heatmap <- function(TCRH_chainlong, gene_order, titlename) {
  TCRH <- TCRH_chainlong
  TCRH$variable <- factor(TCRH$variable, levels = gene_order)
  split_TCRH <- strsplit(TCRH$Sample, "_")
  TCRH$CellType <- sapply(split_TCRH, "[", 1)
  TCRH$AdjusteID <- sapply(split_TCRH, "[", 2)
  TCRH$Organ <- substr(TCRH$AdjusteID, 1, 1)
  TCRH$Sample <- factor(TCRH$Sample, levels = rev(sort(unique(TCRH$Sample))))
  TCRH2 <- TCRH %>%
    group_by(variable) %>%
    filter(sum(as.numeric(value)) != 0) %>%
    ungroup()
  TCRH2$CellType_Organ <- paste(TCRH2$CellType, TCRH2$Organ, sep = "_")

  wide_data <- spread(TCRH2, variable, value)
  rownames(wide_data) <- wide_data$Sample
  wide_data$Sample <- NULL
  wide_data <- wide_data %>%
    mutate(MainID_num = as.numeric(str_extract(AdjusteID, "\\d+\\.\\d")))
  wide_data$Organ <- factor(
    wide_data$Organ, levels = c("B", "L", "T", "S"), ordered = TRUE
  )
  wide_data <- wide_data %>%
    arrange(CellType, Organ, MainID_num) %>%
    select(-MainID_num)
  mat <- sapply(wide_data[, -c(1:4)], as.numeric)
  rownames(mat) <- wide_data$AdjusteID

  annotations2 <- wide_data %>% select(CellType, Organ)
  annotations2$Organ <- factor(
    annotations2$Organ,
    levels = c("B", "L", "T", "S"),
    labels = c("PBMC", "Liver", "Thymus", "Spleen"),
    ordered = TRUE
  )
  annotation_colors2 <- list(
    CellType = c("Naïve CD4 T" = "#fa6e01", "Naïve CD8 T" = "#993333"),
    Organ = c(
      PBMC = "#F6313E", Liver = "#fb862b", Thymus = "#0eb0c8",
      Spleen = "#6a73cf"
    )
  )
  row_anno3 <- rowAnnotation(df = annotations2, col = annotation_colors2)
  row_split_vector <- interaction(
    wide_data$CellType, wide_data$Organ, drop = TRUE, lex.order = TRUE
  )
  ht <- Heatmap(
    mat,
    name = "value",
    col = c("#313695", "#649AC7", "#FFFFBF", "#FDBE70", "#EA5839", "#A50026"),
    cluster_rows = FALSE,
    row_title = NULL,
    cluster_columns = FALSE,
    column_title = titlename,
    row_split = row_split_vector,
    show_row_names = FALSE
  )
  list(row_anno3, ht)
}

save_heatmap <- function(ht, pdf_path, png_path, width, height) {
  pdf(pdf_path, width = width, height = height)
  draw(ht)
  dev.off()
  png(png_path, width = width, height = height, units = "in", res = 300)
  draw(ht)
  dev.off()
}

tcr <- read_csv(
  file.path(input_dir, "M04e_complete_per_sample_segment_gene_usage.csv"),
  show_col_types = FALSE
)

tcr_gene_root <- file.path(.scf_shared, "00_Common/03_Receptor_Gene_Order", "TCR")
tra_locus <- read_csv(file.path(tcr_gene_root, "TRA_Gene_Locus.csv"), show_col_types = FALSE)
trb_locus <- read_csv(file.path(tcr_gene_root, "TRB_Gene_Locus.csv"), show_col_types = FALSE)
orders_tcr <- list(
  TRAV = tra_locus$IMGT_gene_name[grepl("^TRAV", tra_locus$IMGT_gene_name)],
  TRAJ = tra_locus$IMGT_gene_name[grepl("^TRAJ", tra_locus$IMGT_gene_name)],
  TRBV = trb_locus$IMGT_gene_name[grepl("^TRBV", trb_locus$IMGT_gene_name)],
  TRBD = trb_locus$IMGT_gene_name[grepl("^TRBD", trb_locus$IMGT_gene_name)],
  TRBJ = trb_locus$IMGT_gene_name[grepl("^TRBJ", trb_locus$IMGT_gene_name)]
)
tcr_long <- lapply(names(orders_tcr), function(seg) prepare_long(tcr, seg, "keep"))
names(tcr_long) <- names(orders_tcr)
tcr_ht <- Map(
  function(dat, ord, nm) TCRgene_heatmap(dat, ord, nm),
  tcr_long, orders_tcr, names(orders_tcr)
)

m04e_dir <- file.path(.scf_project_root, "Main_Figure_04", "e_M04e", "02_Figures")
dir.create(m04e_dir, recursive = TRUE, showWarnings = FALSE)
ht_tra <- tcr_ht$TRAV[[1]] + tcr_ht$TRAV[[2]] + tcr_ht$TRAJ[[2]]
ht_trb <- tcr_ht$TRBV[[1]] + tcr_ht$TRBV[[2]] + tcr_ht$TRBD[[2]] + tcr_ht$TRBJ[[2]]
ht_tcr <- tcr_ht$TRAV[[1]] + tcr_ht$TRAV[[2]] + tcr_ht$TRAJ[[2]] +
  tcr_ht$TRBV[[2]] + tcr_ht$TRBD[[2]] + tcr_ht$TRBJ[[2]]
save_heatmap(
  ht_tcr,
  file.path(m04e_dir, "M04e_TRAV_TRAJ_TRBV_TRBD_TRBJ_gene_usage.pdf"),
  file.path(m04e_dir, "M04e_TRAV_TRAJ_TRBV_TRBD_TRBJ_gene_usage.png"),
  30, 6
)
if (!file.exists(file.path(m04e_dir, "M04e_TRBVDJ_gene_usage.png"))) {
  save_heatmap(
    ht_tra,
    file.path(m04e_dir, "M04e_TRAVJ_gene_usage.pdf"),
    file.path(m04e_dir, "M04e_TRAVJ_gene_usage.png"),
    16, 6
  )
  save_heatmap(
    ht_trb,
    file.path(m04e_dir, "M04e_TRBVDJ_gene_usage.pdf"),
    file.path(m04e_dir, "M04e_TRBVDJ_gene_usage.png"),
    14, 6
  )
  for (seg in names(tcr_ht)) {
    save_heatmap(
      tcr_ht[[seg]][[1]] + tcr_ht[[seg]][[2]],
      file.path(m04e_dir, paste0("M04e_", seg, "_gene_usage.pdf")),
      file.path(m04e_dir, paste0("M04e_", seg, "_gene_usage.png")),
      max(4, length(orders_tcr[[seg]]) / 5), 6
    )
  }
}
write_csv(bind_rows(tcr_long), file.path(m04e_dir, "M04e_gene_usage_plot_data.csv"))

if (identical(Sys.getenv("M04E_ONLY"), "1")) {
  message("M04e complete")
  quit(save = "no", status = 0)
}

bcr <- read_csv(
  file.path(input_dir, "M05e_complete_per_sample_segment_gene_usage.csv"),
  show_col_types = FALSE
)

bcr_gene_root <- file.path(.scf_shared, "00_Common/03_Receptor_Gene_Order", "BCR")
read_gene_order <- function(file, sep = "\t") {
  tab <- read.table(file.path(bcr_gene_root, file), header = TRUE, sep = sep,
                    check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  as.character(tab[[1]])
}
orders_bcr <- list(
  IGHV = read_gene_order("IGHV.txt"),
  IGHJ = paste0("IGHJ", 1:6),
  IGHD = read_gene_order("IGHD.csv", ","),
  IGKV = read_gene_order("IGKV.txt"),
  IGKJ = read_gene_order("IGKJ.txt"),
  IGLV = read_gene_order("IGLV.txt"),
  IGLJ = read_gene_order("IGLJ.txt")
)
bcr_long <- lapply(names(orders_bcr), function(seg) prepare_long(bcr, seg, "merge_b"))
names(bcr_long) <- names(orders_bcr)

BCRTR_heatmap <- function(BCRH_chainlong, gene_order, titlename) {
  BCRH <- BCRH_chainlong
  BCRH$variable <- factor(BCRH$variable, levels = gene_order)
  split_BCRH <- strsplit(BCRH$Sample, "_")
  BCRH$CellType <- sapply(split_BCRH, "[", 1)
  BCRH$AdjusteID <- sapply(split_BCRH, "[", 2)
  BCRH$Organ <- substr(BCRH$AdjusteID, 1, 1)
  BCRH$Sample <- factor(BCRH$Sample, levels = rev(sort(unique(BCRH$Sample))))
  BCRH2 <- BCRH %>%
    group_by(variable) %>%
    filter(sum(as.numeric(value)) != 0) %>%
    ungroup()
  BCRH2$CellType_Organ <- paste(BCRH2$CellType, BCRH2$Organ, sep = "_")

  wide_data <- spread(BCRH2, variable, value)
  rownames(wide_data) <- wide_data$Sample
  wide_data$Sample <- NULL
  wide_data <- wide_data %>%
    mutate(MainID_num = as.numeric(str_extract(AdjusteID, "\\d+\\.\\d")))
  wide_data$Organ <- factor(
    wide_data$Organ, levels = c("B", "L", "T", "S"), ordered = TRUE
  )
  wide_data <- wide_data %>%
    arrange(CellType, Organ, MainID_num) %>%
    select(-MainID_num)
  matrix_data <- data.frame(
    lapply(wide_data[, -c(1:4), drop = FALSE], as.numeric),
    check.names = FALSE
  )
  mat <- as.matrix(matrix_data)
  rownames(mat) <- wide_data$AdjusteID

  annotations2 <- wide_data %>% select(CellType, Organ)
  annotations2$Organ <- factor(
    annotations2$Organ,
    levels = c("B", "L", "T", "S"),
    labels = c("PBMC", "Liver", "Thymus", "Spleen"),
    ordered = TRUE
  )
  annotation_colors2 <- list(
    Organ = c(
      PBMC = "#F6313E", Liver = "#fb862b", Thymus = "#0eb0c8",
      Spleen = "#6a73cf"
    ),
    CellType = c("Naïve B" = "#0081C9")
  )
  row_anno3 <- rowAnnotation(df = annotations2, col = annotation_colors2)
  ht <- Heatmap(
    mat,
    name = "value",
    col = c("#313695", "#649AC7", "#FFFFBF", "#FDBE70", "#EA5839", "#A50026"),
    cluster_rows = FALSE,
    row_title = NULL,
    cluster_columns = FALSE,
    column_title = titlename,
    show_row_names = TRUE
  )
  list(row_anno3, ht)
}

bcr_ht <- Map(
  function(dat, ord, nm) {
    message("Building BCR gene-usage segment: ", nm)
    BCRTR_heatmap(dat, ord, titlename = nm)
  },
  bcr_long, orders_bcr, names(orders_bcr)
)

m05e_dir <- file.path(.scf_project_root, "Main_Figure_05", "e_M05e", "02_Figures")
dir.create(m05e_dir, recursive = TRUE, showWarnings = FALSE)
ht_bcr <- bcr_ht$IGHV[[1]] + bcr_ht$IGHV[[2]] + bcr_ht$IGHJ[[2]] +
  bcr_ht$IGHD[[2]] + bcr_ht$IGKV[[2]] + bcr_ht$IGKJ[[2]] +
  bcr_ht$IGLV[[2]] + bcr_ht$IGLJ[[2]]
save_heatmap(
  ht_bcr,
  file.path(m05e_dir, "M05e_BCR_all_gene_usage.pdf"),
  file.path(m05e_dir, "M05e_BCR_all_gene_usage.png"),
  30, 6
)
for (seg in names(bcr_ht)) {
  save_heatmap(
    bcr_ht[[seg]][[1]] + bcr_ht[[seg]][[2]],
    file.path(m05e_dir, paste0("M05e_", seg, "_gene_usage.pdf")),
    file.path(m05e_dir, paste0("M05e_", seg, "_gene_usage.png")),
    max(4, length(orders_bcr[[seg]]) / 5), 6
  )
}
write_csv(bind_rows(bcr_long), file.path(m05e_dir, "M05e_gene_usage_plot_data.csv"))

message("M04e and M05e complete")
