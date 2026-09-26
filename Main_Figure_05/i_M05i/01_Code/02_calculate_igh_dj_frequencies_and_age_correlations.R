.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(character())

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- c(file.path(.scf_shared, "Main_Figure_04_05/05_IGH_Junction_Annotations/M05i_IGH_Junction_Annotations.txt"),
                             file.path(.scf_start_dir, "../03_Source_Data/BCR_cell_metadata.csv"),
                             file.path(.scf_start_dir, "../03_PlotData"))
if (length(args) != 3L) stop("Usage: Rscript 02_calculate_igh_dj_frequencies_and_age_correlations.R M05i_IGH_Junction_Annotations.txt BCRobs.csv output_directory")
stopifnot(all(file.exists(args[1:2])))
junction <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
obs <- read.csv(args[2], stringsAsFactors = FALSE)
stopifnot(all(c("Sequence ID", "D-GENE and allele", "J-GENE and allele") %in% names(junction)))
stopifnot(all(c("Cellname", "MainID") %in% names(obs)))
junction$Clean_IGHD <- gsub("^Homsap\\s+|\\*.*$", "", junction[["D-GENE and allele"]])
junction$Clean_IGHJ <- gsub("^Homsap\\s+|\\*.*$", "", junction[["J-GENE and allele"]])
junction <- junction[!is.na(junction$Clean_IGHD) & junction$Clean_IGHD != "", ]
joined <- merge(junction[, c("Sequence ID", "Clean_IGHD", "Clean_IGHJ")], obs,
                by.x = "Sequence ID", by.y = "Cellname")
retained <- joined
pbmc <- retained[grep("B", retained$MainID), ]
stopifnot(nrow(pbmc) > 0)
pbmc$combo <- paste0(pbmc$Clean_IGHD, "_", pbmc$Clean_IGHJ)
counts <- as.data.frame(table(pbmc$MainID, pbmc$combo), stringsAsFactors = FALSE)
names(counts) <- c("MainID", "combo", "frequency")
counts$sample_total <- ave(counts$frequency, counts$MainID, FUN = sum)
counts$percentage <- counts$frequency / counts$sample_total * 100
counts$Post_Conception_Age_Weeks <- as.numeric(substring(counts$MainID, 2, 5))
stopifnot(!anyNA(counts$Post_Conception_Age_Weeks), all(counts$sample_total > 0))
statistics <- do.call(rbind, lapply(split(counts, counts$combo), function(x) {
  test <- suppressWarnings(cor.test(x$Post_Conception_Age_Weeks, x$percentage, method = "spearman"))
  data.frame(combo = x$combo[1], n_samples = nrow(x), rho = unname(test$estimate), p_value = test$p.value)
}))
statistics$p_adjusted_bh <- p.adjust(statistics$p_value, method = "BH")
statistics$age_basis <- "Post-conception age parsed from MainID"
statistics$direction <- ifelse(is.na(statistics$rho), "", ifelse(statistics$rho > 0, "+", "-"))
statistics$significant_raw_p_0_05 <- !is.na(statistics$p_value) & statistics$p_value < .05
statistics$significant_bh_q_0_05 <- !is.na(statistics$p_adjusted_bh) & statistics$p_adjusted_bh <= .05
statistics$marker <- ifelse(is.na(statistics$p_adjusted_bh) | statistics$p_adjusted_bh > 0.05,
                            "", ifelse(statistics$rho > 0, "+", "-"))
overall <- as.data.frame(table(pbmc$combo), stringsAsFactors = FALSE)
names(overall) <- c("combo", "frequency")
overall$percentage <- overall$frequency / sum(overall$frequency) * 100
d_order <- c("IGHD1-1", "IGHD2-2", "IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6",
             "IGHD1-7", "IGHD2-8", "IGHD3-9", "IGHD3-10", "IGHD4-11", "IGHD5-12",
             "IGHD6-13", "IGHD1-14", "IGHD2-15", "IGHD3-16", "IGHD4-17", "IGHD5-18",
             "IGHD6-19", "IGHD1-20", "IGHD2-21", "IGHD3-22", "IGHD4-23", "IGHD5-24",
             "IGHD6-25", "IGHD1-26", "IGHD7-27")
grid <- expand.grid(Dgene = d_order, Jgene = paste0("IGHJ", 1:6), stringsAsFactors = FALSE)
grid$combo <- paste0(grid$Dgene, "_", grid$Jgene)
plot_data <- merge(grid, overall, by = "combo", all.x = TRUE, sort = FALSE)
plot_data <- merge(plot_data, statistics, by = "combo", all.x = TRUE, sort = FALSE)
plot_data$frequency[is.na(plot_data$frequency)] <- 0
plot_data$percentage[is.na(plot_data$percentage)] <- 0
plot_data$marker[is.na(plot_data$marker)] <- ""
plot_data <- plot_data[order(match(plot_data$Jgene, paste0("IGHJ", 1:6)), match(plot_data$Dgene, d_order)), ]
stopifnot(abs(sum(overall$percentage) - 100) < 1e-8)
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
write.csv(plot_data, file.path(args[3], "M05i_IGH_DJ_plot_data.csv"), row.names = FALSE)
write.csv(statistics, file.path(args[3], "M05i_IGH_DJ_all_statistics.csv"), row.names = FALSE)
write.csv(counts, file.path(args[3], "M05i_IGH_DJ_sample_percentages.csv"), row.names = FALSE)
write.csv(pbmc[, c("Sequence ID", "MainID", "Clean_IGHD", "Clean_IGHJ")],
          file.path(args[3], "M05i_IGH_DJ_cell_input.csv"), row.names = FALSE)
