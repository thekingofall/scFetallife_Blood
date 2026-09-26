.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages("readxl")
panel <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
d <- as.data.frame(readxl::read_excel(file.path(panel, "03_PlotData/M02d_Selected_GO_Terms.xlsx")))
d <- d[!is.na(d$GO_ID), ]
modules <- c("M1-ERY", "M2-Mono", "M3-MK", "M5-Bcell", "M6-NK", "M7-Tcell")
for (field in c("Raw_P", "BH_adjusted_P", "GeneCount", "InputGeneCount", "UniverseGeneCount")) d[[field]] <- as.numeric(d[[field]])
stopifnot(all(d$BH_adjusted_P > 0 & d$BH_adjusted_P < 0.05))
d$SourceRow <- seq_len(nrow(d)) + 1L
d$NegLog10BH <- -log10(d$BH_adjusted_P)
d$PlotNegLog10BH <- pmin(d$NegLog10BH, 12)
d <- d[order(match(d$Module, modules), match(d$Stage, c("Early", "Late")), d$BH_adjusted_P, -d$GeneCount, d$Description), ]
d$DisplayOrder <- seq_len(nrow(d))
dir.create(file.path(panel, "03_PlotData"), showWarnings = FALSE)
write.csv(d, file.path(panel, "03_PlotData/M02d_selected_pathways.csv"), row.names = FALSE, na = "")
