.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "grid", "readxl"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")
.scf_panel_root <- normalizePath(file.path(.scf_code_dir, ".."), winslash = "/", mustWork = FALSE)

options(stringsAsFactors = FALSE)
library(ggplot2)
table_dir <- file.path(.scf_panel_root, "03_PlotData")
plot_dir <- file.path(.scf_panel_root, "02_Figures")
dir.create(table_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
}

stop_unless <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

tissue_of <- function(ids) {
  prefix <- substr(ids, 1L, 1L)
  unname(c(B = "PBMC", L = "Liver", T = "Thymus", S = "Spleen")[prefix])
}

denominator_of <- function(tissues) {
  ifelse(tissues == "PBMC",
         "total lymphocytes (FlowJo parent: Single live cells)",
         "total CD45+ cells (FlowJo parent: CD45+cells)")
}

processed_blood_path <- file.path(.scf_shared, "Main_Figure_06/00_Flow_Summaries/Blood/Flow_Blood_CD45Freq.csv")
processed_organ_path <- file.path(.scf_shared, "Main_Figure_06/00_Flow_Summaries/Organ/Flow_organ_FreqofCD45.csv")
manuscript_source_xlsx <- file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_15/a_S15a/manuscript_source_data.xlsx")
blood_all <- c(
  "B12.0_P2", "B17.4_P3", "B18.0_P4", "B18.6_P5", "B20.9_P7",
  "B21.7_P8", "B22.4_P9", "B23.4_P10", "B24.6_P11", "B26.9_P12",
  "B29.1_P13", "B31.3_P15", "B32.4_P16", "B33.3_P17",
  "B34.1_P18", "B36.1_P19", "B37.9_P20", "B38.1_P21", "B39.1_P22"
)

organ_all <- c(
  "S18.6_P5", "S22.0_P26", "S24.6_P27", "S24.6_P11",
  "L10.0_P1", "L18.6_P5", "L22.0_P26", "L24.6_P27", "L24.6_P11",
  "T10.0_P1", "T18.6_P5", "T22.0_P26", "T24.6_P27", "T24.6_P11"
)

selected_samples <- c(
  blood_all,
  "L10.0_P1", "L18.6_P5", "L24.6_P11",
  "T10.0_P1", "T18.6_P5", "T24.6_P11",
  "S18.6_P5", "S24.6_P11"
)


source17 <- c(
  "CXCR5-B", "CXCR5+B", "Immature NK", "Mature NK",
  "Central memory-like CD8+T", "Naive CD8+T",
  "Effective memory-like CD8+T", "Terminal differentiated CD8+T",
  "Central memory-like CD4+T", "Naive CD4+T",
  "Effective memory-like CD4+T", "Terminal differentiated CD4+T",
  "DNT", "LDP", "EDP", "CD8+ISP", "CD4+ISP"
)

class11 <- c(
  "CXCR5-B", "CXCR5+B", "Immature NK", "Mature NK", "CD8+ T",
  "CD4+ T", "DNT", "CD4+ ISP", "CD8+ ISP", "DPT", "Ungate"
)

mapping11 <- list(
  "CXCR5-B" = "CXCR5-B",
  "CXCR5+B" = "CXCR5+B",
  "Immature NK" = "Immature NK",
  "Mature NK" = "Mature NK",
  "CD8+ T" = c(
    "Central memory-like CD8+T", "Naive CD8+T",
    "Effective memory-like CD8+T", "Terminal differentiated CD8+T"
  ),
  "CD4+ T" = c(
    "Central memory-like CD4+T", "Naive CD4+T",
    "Effective memory-like CD4+T", "Terminal differentiated CD4+T"
  ),
  "DNT" = "DNT",
  "CD4+ ISP" = "CD4+ISP",
  "CD8+ ISP" = "CD8+ISP",
  "DPT" = c("LDP", "EDP")
)

stop_unless(identical(sort(unlist(mapping11, use.names = FALSE)), sort(source17)),
            "The ten measured classes must partition all 17 source populations exactly once")

palette11 <- c(
  "CXCR5-B" = "#FF0000",
  "CXCR5+B" = "#FF9B9B",
  "Immature NK" = "#13C0DF",
  "Mature NK" = "#077E97",
  "CD8+ T" = "#FF6000",
  "CD4+ T" = "#B856D7",
  "DNT" = "#E8E800",
  "CD4+ ISP" = "#CE0666",
  "CD8+ ISP" = "#FA6494",
  "DPT" = "#A2A200",
  "Ungate" = "#8C8ED2"
)

read_frequency_matrix <- function(path) {
  x <- read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  names(x)[1] <- "source_feature"
  x$source_feature <- trimws(as.character(x$source_feature))
  x
}

matrix_to_long <- function(x) {
  sample_columns <- names(x)[-1]
  do.call(rbind, lapply(sample_columns, function(sample_id) {
    data.frame(
      MainID = sample_id,
      Tissue = tissue_of(sample_id),
      source_feature = x$source_feature,
      Percentage = as.numeric(x[[sample_id]]),
      stringsAsFactors = FALSE
    )
  }))
}

blood_frequency <- read_frequency_matrix(processed_blood_path)
organ_frequency <- read_frequency_matrix(processed_organ_path)
stop_unless(identical(names(blood_frequency)[-1], blood_all),
            "Processed blood sample order differs from the sample order")
stop_unless(identical(names(organ_frequency)[-1], organ_all),
            "Processed organ sample order differs from the sample order")
stop_unless(all(source17 %in% blood_frequency$source_feature),
            "One or more 17-class source populations are missing from the blood matrix")
stop_unless(all(source17 %in% organ_frequency$source_feature),
            "One or more 17-class source populations are missing from the organ matrix")

frequency_all_long <- rbind(matrix_to_long(blood_frequency), matrix_to_long(organ_frequency))


manuscript_fill_source <- as.data.frame(
  readxl::read_excel(manuscript_source_xlsx, sheet = "Fig. 6b"),
  check.names = FALSE, stringsAsFactors = FALSE
)
names(manuscript_fill_source) <- gsub(
  "Naïve", "Naive", names(manuscript_fill_source), fixed = TRUE
)
manuscript_fill_source <- manuscript_fill_source[
  grepl("^[BLTS][0-9]", as.character(manuscript_fill_source$ID)),
]
manuscript_fill_source$ID[manuscript_fill_source$ID == "T18.6_P6"] <- "T18.6_P5"

missing_source_idx <- which(
  frequency_all_long$source_feature %in% source17 & is.na(frequency_all_long$Percentage)
)

# Read the DNT percentage from the unrounded flow-cytometry table.
for (i in missing_source_idx) {
  j <- match(frequency_all_long$MainID[i], manuscript_fill_source$ID)
  feature <- frequency_all_long$source_feature[i]
  stop_unless(!is.na(j), "Sample is absent from manuscript source data")
  frequency_all_long$Percentage[i] <- as.numeric(manuscript_fill_source[j, feature])
}
source17_all_long <- frequency_all_long[frequency_all_long$source_feature %in% source17, ]
source17_all_long <- source17_all_long[
  order(match(source17_all_long$MainID, c(blood_all, organ_all)),
        match(source17_all_long$source_feature, source17)),
]
stop_unless(nrow(source17_all_long) == length(c(blood_all, organ_all)) * 17L, "Expected 33 samples x 17 source populations")
stop_unless(!anyNA(source17_all_long$Percentage), "Source 17-class matrix contains missing values")
stop_unless(all(source17_all_long$Percentage >= 0), "Source 17-class matrix contains negative values")

build_11class <- function(source_long, ordered_samples) {
  measured <- do.call(rbind, lapply(ordered_samples, function(sample_id) {
    sample_data <- source_long[source_long$MainID == sample_id, ]
    stop_unless(nrow(sample_data) == 17L, paste("Expected 17 source rows for", sample_id))
    do.call(rbind, lapply(names(mapping11), function(cell_class) {
      features <- mapping11[[cell_class]]
      data.frame(
        MainID = sample_id,
        Tissue = tissue_of(sample_id),
        Cell_type = cell_class,
        Percentage = sum(sample_data$Percentage[sample_data$source_feature %in% features]),
        stringsAsFactors = FALSE
      )
    }))
  }))

  residual <- aggregate(Percentage ~ MainID + Tissue, measured, sum)
  residual$Percentage <- 100 - residual$Percentage
  residual$Cell_type <- "Ungate"
  stop_unless(all(residual$Percentage >= -1e-10), "A computed Ungate residual is negative")
  residual$Percentage[abs(residual$Percentage) < 1e-12] <- 0

  result <- rbind(measured, residual[, names(measured)])
  result <- result[
    order(match(result$MainID, ordered_samples), match(result$Cell_type, class11)),
  ]
  rownames(result) <- NULL
  result
}

source17_selected <- source17_all_long[source17_all_long$MainID %in% selected_samples, ]
composition <- build_11class(source17_all_long, selected_samples)
stopifnot(nrow(composition) == length(selected_samples) * length(class11))
stopifnot(all(abs(aggregate(Percentage ~ MainID, composition, sum)$Percentage - 100) < 1e-8))
wide_from_long <- function(x, ordered_samples) {
  out <- data.frame(
    MainID = ordered_samples,
    Tissue = tissue_of(ordered_samples),
    Denominator = denominator_of(tissue_of(ordered_samples)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (cell_class in class11) {
    values <- x$Percentage[x$Cell_type == cell_class]
    ids <- x$MainID[x$Cell_type == cell_class]
    out[[cell_class]] <- values[match(ordered_samples, ids)]
  }
  out
}

write_tsv(source17_selected, file.path(table_dir, "S15a_source17_long.tsv"))
write_tsv(composition, file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_15/a_S15a/03_PlotData/S15a_11class_long.tsv"))
write_tsv(wide_from_long(composition, selected_samples), file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_15/a_S15a/03_PlotData/S15a_11class_wide.tsv"))
plot_data <- composition
plot_data$MainID <- factor(plot_data$MainID, levels = selected_samples)
plot_data$Tissue <- factor(plot_data$Tissue,
                           levels = c("PBMC", "Liver", "Thymus", "Spleen"))
plot_data$Cell_type <- factor(plot_data$Cell_type, levels = class11)

p <- ggplot(plot_data, aes(x = MainID, y = Percentage, fill = Cell_type)) +
  geom_col(width = 0.92, colour = NA) +
  facet_grid(. ~ Tissue, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette11, breaks = class11, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 100, 25), limits = c(0, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Spectral flow cytometry",
    x = NULL, y = "Composition (%)", fill = "Cell_type"
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = FALSE)) +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 21, face = "plain", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "#4D4D4D"),
    axis.text.y = element_text(colour = "#4D4D4D"),
    axis.title.y = element_text(size = 15),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    panel.spacing.x = grid::unit(0.35, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.key.height = grid::unit(0.55, "lines"),
    plot.margin = margin(8, 10, 8, 8)
  )

pdf_path <- file.path(plot_dir, "FigureS15a_SFC_11Class_Composition.pdf")
png_path <- file.path(plot_dir, "FigureS15a_SFC_11Class_Composition.png")
ggsave(pdf_path, p, width = 16, height = 6.4, units = "in", device = cairo_pdf)
ggsave(png_path, p, width = 16, height = 6.4, units = "in", dpi = 300, bg = "white")
