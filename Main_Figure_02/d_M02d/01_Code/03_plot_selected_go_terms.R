.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "ragg"))
suppressPackageStartupMessages(library(ggplot2))
panel <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
figures <- file.path(panel, "02_Figures")
dir.create(figures, showWarnings = FALSE)
palette_table <- read.delim(file.path(.scf_project_root, "00_Settings/palette_registry.tsv"))
palette_table <- palette_table[palette_table$scope == "M02d", ]
palette <- setNames(palette_table$hex, palette_table$element)
d <- read.csv(file.path(panel, "03_PlotData/M02d_selected_pathways.csv"), check.names = FALSE)
layout_capture <- ragg::agg_capture(width = 3900, height = 3780, res = 300)
modules <- c("M1-ERY", "M2-Mono", "M3-MK", "M5-Bcell", "M6-NK", "M7-Tcell")
module_labels <- c("M1 Ery", "M2 Mono", "M3 MK", "M5 B", "M6 NK", "M7 T")
stopifnot(all(is.finite(d$BH_adjusted_P) & d$BH_adjusted_P > 0 & d$BH_adjusted_P < 0.05))
d <- d[order(d$DisplayOrder), ]
d$Stage <- factor(d$Stage, c("Early", "Late"))
d$ModuleLabel <- factor(module_labels[match(d$Module, modules)], module_labels)
# Module-qualified labels retain repeated biological terms in their own panels.
d$TermKey <- paste(d$Module, d$GO_ID, sep = "::")
stopifnot(!anyDuplicated(d$TermKey))
labels <- setNames(vapply(d$Description, function(x) paste(strwrap(x, width = 47), collapse = "\n"), character(1)), d$TermKey)
d$TermKey <- factor(d$TermKey, rev(d$TermKey))
p <- ggplot(d, aes(Stage, TermKey)) +
  geom_point(aes(size = GeneCount, colour = PlotNegLog10BH), alpha = 0.94) +
  facet_grid(ModuleLabel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = labels, expand = expansion(add = 0.6)) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_colour_gradientn(colours = unname(palette[c("Lower significance", "Intermediate significance", "Higher significance")]),
                         limits = c(0, 12), breaks = c(0, 4, 8, 12), labels = c("0", "4", "8", "\u226512"), name = "-log10(BH P)") +
  scale_size_continuous(range = c(5.0, 10.0), breaks = c(10, 20, 30, 40), name = "Gene count") +
  guides(colour = guide_colourbar(order = 1, barheight = grid::unit(34, "mm"), barwidth = grid::unit(5, "mm")),
         size = guide_legend(order = 2)) +
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL) +
  theme_classic(base_size = 18, base_family = "Arial") +
  theme(text = element_text(family = "Arial", colour = "#202124"),
        axis.line = element_blank(), axis.ticks.length = grid::unit(1.9, "mm"),
        legend.key = element_blank(), legend.box.margin = margin(0, 0, 0, 4),
        plot.background = element_rect(fill = "white", colour = NA)) +
  theme(panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.35),
        panel.grid.major.x = element_line(colour = "#E5E5E5", linewidth = 0.35),
        panel.grid.major.y = element_line(colour = "#E5E5E5", linewidth = 0.35),
        panel.border = element_rect(colour = palette[["Frames and text"]], fill = NA, linewidth = 0.8),
        panel.spacing.y = grid::unit(2.2, "mm"),
        strip.background = element_rect(fill = "#E7E7E7", colour = palette[["Frames and text"]], linewidth = 0.8),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 18, margin = margin(8, 12, 8, 12)),
        axis.text.y = element_text(size = 16, colour = "black", lineheight = 0.95, margin = margin(r = 5)),
        axis.text.x = element_text(size = 18, face = "bold", colour = "black"),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
        legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 15),
        legend.box.spacing = grid::unit(4, "mm"), legend.spacing.y = grid::unit(5, "mm"),
        plot.margin = margin(9, 12, 9, 9))
ggsave(file.path(figures, "M02d_Hotspot_GO.png"), p, width = 13, height = 12.6, dpi = 300, device = ragg::agg_png, bg = "white")
ggsave(file.path(figures, "M02d_Hotspot_GO.pdf"), p, width = 13, height = 12.6, device = cairo_pdf, bg = "white")
invisible(dev.off())
