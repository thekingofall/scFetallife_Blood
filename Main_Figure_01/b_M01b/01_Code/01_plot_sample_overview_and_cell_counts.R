.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2"))

suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
if (length(script_arg) != 1) stop("Unable to locate the rendering script")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork=TRUE)
panel_dir <- dirname(dirname(script_path))
data_dir <- file.path(panel_dir, "03_PlotData")
figure_dir <- file.path(panel_dir, "02_Figures")
dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)

plot_data <- read.csv(
  file.path(data_dir, "M01b_multimodal_plot_data.csv"),
  stringsAsFactors=FALSE,
  check.names=FALSE
)
count_data <- read.delim(
  file.path(data_dir, "M01b_cell_counts.tsv"),
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_plot_columns <- c(
  "MainID", "AlmostWeek", "Organ", "has_SCRNAseq", "has_Spectral",
  "has_Plasma_Olink", "has_Olink_stimulation", "hasRNA"
)
required_count_columns <- c("Organ", "display_label", "cell_count")
stopifnot(all(required_plot_columns %in% names(plot_data)))
stopifnot(all(required_count_columns %in% names(count_data)))

sample_metadata <- read.csv(
  file.path(panel_dir, "../../00_Global_Data/Fetal_Immune_Atlas_Sample_Metadata.csv"),
  stringsAsFactors=FALSE,
  check.names=FALSE
)
assay_columns <- c(
  has_SCRNAseq="scRNAseq_Available",
  has_Spectral="Flow_Cytometry_Available",
  has_Plasma_Olink="Plasma_Olink_Available",
  has_Olink_stimulation="Stimulated_Olink_Available",
  hasRNA="PHA_Bulk_RNA_Available"
)
stopifnot(all(c("MainID", "Main_Organ", unname(assay_columns)) %in% names(sample_metadata)))
sample_index <- match(plot_data$MainID, sample_metadata$MainID)
stopifnot(!anyNA(sample_index), !anyDuplicated(sample_metadata$MainID))
stopifnot(all(plot_data$Organ == sample_metadata$Main_Organ[sample_index]))
for (column in names(assay_columns)) {
  availability <- sample_metadata[[assay_columns[[column]]]][sample_index]
  stopifnot(!anyNA(availability), all(availability %in% c(0, 1)))
  plot_data[[column]] <- availability == 1
}
plot_data$n_techs <- rowSums(plot_data[names(assay_columns)])
write.csv(plot_data, file.path(data_dir, "M01b_multimodal_plot_data.csv"), row.names=FALSE)

organ_levels <- c("PBMC", "Liver", "Thymus", "Spleen")
organ_centers <- c(PBMC=3.70, Liver=2.70, Thymus=1.70, Spleen=0.70)
count_data <- count_data[match(organ_levels, count_data$Organ), ]
count_data$center <- unname(organ_centers[count_data$Organ])
count_data$count_label <- format(
  count_data$cell_count,
  big.mark=",",
  scientific=FALSE,
  trim=TRUE
)

palette <- c(
  scrna="#ED3734",
  spectral="#F3A140",
  plasma="#7268AF",
  stimulation_proteomics="#7268AF",
  stimulation_bulk_rna="#43AE78",
  organ_strip="#57B8BF",
  frame="#202020",
  grid="#C9C9C9"
)

modality_specs <- data.frame(
  column=c(
    "has_SCRNAseq", "has_Spectral", "has_Plasma_Olink",
    "has_Olink_stimulation", "hasRNA"
  ),
  label=c(
    "scRNA-seq\nsc αβTCR/BCR-seq", "Spectral flow\ncytometry",
    "Plasma\nproteomics", "Stimulation\nproteomics",
    "Stimulation\nbulk RNA-seq"
  ),
  colour=unname(palette[c(
    "scrna", "spectral", "plasma", "stimulation_proteomics",
    "stimulation_bulk_rna"
  )]),
  shape=c(15, 15, 15, 16, 16),
  y_offset=c(-0.32, -0.16, 0.00, 0.16, 0.32),
  legend_x=c(11.1, 17.8, 24.5, 31.3, 37.6),
  stringsAsFactors=FALSE
)

point_data <- do.call(rbind, lapply(seq_len(nrow(modality_specs)), function(i) {
  spec <- modality_specs[i, ]
  keep <- as.logical(plot_data[[spec$column]])
  d <- plot_data[keep, c("MainID", "AlmostWeek", "Organ")]
  d$colour <- spec$colour
  d$shape <- spec$shape
  d$y <- unname(organ_centers[d$Organ]) + spec$y_offset
  d
}))

row_half_height <- 0.45
row_rectangles <- data.frame(
  Organ=organ_levels,
  center=unname(organ_centers[organ_levels]),
  stringsAsFactors=FALSE
)
row_rectangles$ymin <- row_rectangles$center - row_half_height
row_rectangles$ymax <- row_rectangles$center + row_half_height

grid_data <- expand.grid(
  x=9:40,
  Organ=organ_levels,
  stringsAsFactors=FALSE
)
grid_data$center <- unname(organ_centers[grid_data$Organ])
grid_data$ymin <- grid_data$center - row_half_height
grid_data$ymax <- grid_data$center + row_half_height

font_size_pt <- 11
font_size_mm <- font_size_pt / ggplot2::.pt

p <- ggplot() +
  geom_rect(
    data=row_rectangles,
    aes(xmin=8.60, xmax=40.55, ymin=ymin, ymax=ymax),
    fill="white", colour=palette[["frame"]], linewidth=0.30
  ) +
  geom_rect(
    data=row_rectangles,
    aes(xmin=6.00, xmax=8.60, ymin=ymin, ymax=ymax),
    fill=palette[["organ_strip"]], colour=NA
  ) +
  geom_rect(
    data=row_rectangles,
    aes(xmin=40.55, xmax=43.15, ymin=ymin, ymax=ymax),
    fill=palette[["organ_strip"]], colour=NA
  ) +
  geom_segment(
    data=grid_data,
    aes(x=x, xend=x, y=ymin, yend=ymax),
    colour=palette[["grid"]], linewidth=0.35, linetype="dotted"
  ) +
  geom_point(
    data=point_data,
    aes(x=AlmostWeek, y=y, colour=colour, shape=shape),
    size=2.65, stroke=0
  ) +
  scale_colour_identity() +
  scale_shape_identity() +
  geom_text(
    data=count_data,
    aes(x=7.30, y=center, label=display_label),
    family="Arial", size=font_size_mm, colour=palette[["frame"]]
  ) +
  geom_text(
    data=count_data,
    aes(x=41.85, y=center, label=count_label),
    family="Arial", size=font_size_mm, colour=palette[["frame"]]
  ) +
  geom_point(
    data=modality_specs,
    aes(x=legend_x, y=4.80, colour=colour, shape=shape),
    size=4.60, stroke=0
  ) +
  geom_text(
    data=modality_specs,
    aes(x=legend_x, y=4.41, label=label),
    family="Arial", size=font_size_mm, lineheight=0.88,
    colour=palette[["frame"]]
  ) +
  annotate(
    "text", x=41.85, y=4.35, label="No. of cells",
    family="Arial", size=font_size_mm, colour=palette[["frame"]]
  ) +
  geom_segment(
    data=data.frame(x=9:40),
    aes(x=x, xend=x, y=0.25, yend=0.20),
    colour=palette[["frame"]], linewidth=0.30
  ) +
  geom_text(
    data=data.frame(x=9:40),
    aes(x=x, y=0.08, label=x),
    family="Arial", size=font_size_mm, colour=palette[["frame"]]
  ) +
  annotate(
    "text", x=24.55, y=-0.17, label="pcw",
    family="Arial", size=font_size_mm,
    colour=palette[["frame"]]
  ) +
  coord_cartesian(xlim=c(5.40, 44.65), ylim=c(-0.745, 5.115), expand=FALSE) +
  theme_void(base_family="Arial", base_size=font_size_pt) +
  theme(plot.margin=margin(0, 0, 0, 0, unit="mm"))

pdf_path <- file.path(figure_dir, "M01b_sample_multiomics_overview.pdf")
png_path <- file.path(figure_dir, "M01b_sample_multiomics_overview.png")

cairo_pdf(pdf_path, width=10.6, height=5.22, family="Arial", onefile=TRUE)
print(p)
invisible(dev.off())

png(
  png_path,
  width=10.6,
  height=5.22,
  units="in",
  res=300,
  type="cairo-png",
  family="Arial",
  antialias="subpixel"
)
print(p)
invisible(dev.off())

cat(sprintf(
  "Blood cells: %s; total cells: %s\nSaved: %s\nSaved: %s\n",
  format(count_data$cell_count[count_data$Organ == "PBMC"], big.mark=","),
  format(sum(count_data$cell_count), big.mark=","),
  pdf_path,
  png_path
))
