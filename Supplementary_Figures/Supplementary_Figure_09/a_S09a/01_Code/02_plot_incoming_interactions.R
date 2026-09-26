.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "ragg"))
suppressPackageStartupMessages(library(ggplot2))
panel <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
root <- normalizePath(file.path(.scf_start_dir, "../../../.."), winslash = "/")
plot_data <- file.path(panel, "03_PlotData")
dir.create(plot_data, showWarnings = FALSE)
d <- readRDS(file.path(plot_data, "Incoming_Interaction_Plot.rds"))
source_order <- readRDS(file.path(root, "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_09/01_CellChat_Models/S09_CellChat_Communication_Networks.rds"))$sources
d$source <- factor(d$source, levels = source_order)
d$stage <- factor(sub("PBMC_(Early|Late).*", "\\1", d$dataset), levels = c("Early", "Late"))
d$target_label <- "Naïve CD8⁺ T"
colors <- colorRampPalette(c("#2166ac", "#518bbd", "#F0E68C", "#d25e55", "#b2182b"))(100)
p <- ggplot(d, aes(source, interaction_name, size = significance, color = prob)) +
  geom_point() + theme_bw(base_family = "Arial") +
  theme(axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.key = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 12),
        strip.background.y = element_blank(),
        strip.background.x = element_rect(fill = "white", color = "black"),
        strip.text.x = element_text(size = 12),
        axis.title.y = element_blank()) +
  facet_grid(stage ~ target_label, space = "free", scales = "free", switch = "y") +
  scale_color_gradientn(colours = colors, name = "Communication\nprobability") + labs(x = NULL, size = "P value")
figures <- file.path(panel, "02_Figures")
dir.create(figures, showWarnings = FALSE)
ggsave(file.path(figures, "S09a_Incoming_Ligand_Receptor_Interactions.pdf"), p, width = 12, height = 12, device = cairo_pdf)
ggsave(file.path(figures, "S09a_Incoming_Ligand_Receptor_Interactions.png"), p, width = 12, height = 12, device = ragg::agg_png, dpi = 300, bg = "white")
