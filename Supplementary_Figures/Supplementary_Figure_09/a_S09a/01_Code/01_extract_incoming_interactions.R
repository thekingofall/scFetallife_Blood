.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("CellChat"))
suppressPackageStartupMessages(library(CellChat))
panel <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
root <- normalizePath(file.path(.scf_start_dir, "../../../.."), winslash = "/")
plot_data <- file.path(panel, "03_PlotData")
dir.create(plot_data, showWarnings = FALSE)
saved <- readRDS(file.path(root, "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_09/01_CellChat_Models/S09_CellChat_Communication_Networks.rds"))
incoming <- netVisual_bubble(saved$network, sources.use = saved$sources,
               targets.use = "Naïve CD8 T", comparison = c(2, 1), remove.isolate = TRUE)$data
d <- incoming[grepl("MHC|LCK|SELPLG|ITGB2", incoming$pathway_name), ]
d <- d[!grepl("Megakaryocytes|ERY",d$source), ]
d$dataset <- paste0(d$dataset, ":Naïve CD8 T")
d$significance <- ifelse(d$pval == 3, "P < 0.01", "0.01 <= P < 0.05")
d$significance <- factor(d$significance, levels = c("0.01 <= P < 0.05", "P < 0.01"))
saveRDS(d, file.path(plot_data, "Incoming_Interaction_Plot.rds"))
write.csv(d, file.path(plot_data, "Incoming_Interaction_Plot.csv"), row.names = FALSE)
