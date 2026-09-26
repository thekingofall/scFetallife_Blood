.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "grid"))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- c(file.path(.scf_start_dir, "../03_PlotData/M05i_IGH_DJ_plot_data.csv"), file.path(.scf_start_dir, "../02_Figures"))
stopifnot(length(args) == 2L, file.exists(args[1]))
suppressPackageStartupMessages(library(ggplot2))
d <- read.csv(args[1], stringsAsFactors = FALSE)
d_order <- unique(d$Dgene)
d$x <- match(d$Dgene, d_order)
d$y <- 7 - as.integer(sub("IGHJ", "", d$Jgene))
s <- d[!is.na(d$marker) & d$marker != "", ]
s$symbol_size <- c(4, 5, 6, 7, 9)[findInterval(abs(s$rho), c(-Inf, .5, .6, .7, .8, Inf), rightmost.closed = TRUE)]
p <- ggplot(d, aes(x, y)) +
  geom_tile(aes(fill = percentage), color = "white", linewidth = .12) +
  geom_text(data = s, aes(label = marker, color = marker, size = symbol_size),
            family = "Arial", fontface = "bold", show.legend = TRUE) +
  scale_size_identity(guide = "none") +
  scale_fill_gradientn(colors = c("#246BAE", "#FFFFFF", "#EDAE11", "#C71000"),
                       limits = c(0, max(d$percentage)), name = "Percentage (%)") +
  scale_color_manual(values = c("-" = "#3D3B25", "+" = "#8F1336"),
                     breaks = c("-", "+"), labels = c("Negative", "Positive"),
                     name = "Correlation between\nD-J pair and pcw") +
  scale_x_continuous(breaks = seq_along(d_order), labels = d_order, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:6, labels = paste0("IGHJ", 6:1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(.5, 27.5), ylim = c(.5, 6.5), clip = "off") +
  annotate("segment", x = .2, xend = .2, y = 6.1, yend = .9, linewidth = .7,
           arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed")) +
  annotate("segment", x = 1, xend = 27, y = .29, yend = .29, linewidth = .7,
           arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed")) +
  annotate("text", x = c(.2,.2,.6,27.5), y = c(6.5,.5,.23,.23),
           label = c("5′", "3′", "5′", "3′"), family = "Arial", size = 4) +
  labs(x = "D gene", y = "J gene", title = "IGH") +
  guides(fill = guide_colorbar(order = 1, barheight = grid::unit(28, "mm"),
                               barwidth = grid::unit(4, "mm")),
         color = guide_legend(order = 2, override.aes = list(label = c("−", "+"), size = 5))) +
  theme_classic(base_size = 15, base_family = "Arial") +
  theme(plot.title = element_text(hjust = .5, size = 19),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5,
                                   color = "black", margin = margin(t = 14)),
        axis.text.y = element_text(size = 13, color = "black", margin = margin(r = 13)),
        axis.title = element_text(size = 16), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.title = element_text(size = 12),
        legend.text = element_text(size = 12), legend.position = "right",
        legend.box.spacing = grid::unit(5, "mm"), plot.margin = margin(10, 10, 8, 10))
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
stem <- file.path(args[2], "M05i_BCRH_all_DJcombined")
ggsave(paste0(stem, ".pdf"), p, width = 9.4, height = 5.0, device = cairo_pdf, bg = "white")
ggsave(paste0(stem, ".png"), p, width = 9.4, height = 5.0, dpi = 300, bg = "white")

