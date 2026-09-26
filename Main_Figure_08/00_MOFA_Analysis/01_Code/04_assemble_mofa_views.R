.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
root <- normalizePath(file.path(.scf_start_dir, "../../.."), winslash = "/")
source(file.path(root, "00_Settings/00_setup_R.R"))
scf_check_R_packages("reshape2")
output <- file.path(.scf_start_dir, "../02_Data")
files <- c("Single_Cell_Features.csv.gz", "Protein_Flow_Bulk_Features.csv.gz", "Repertoire_Features.csv.gz")
data <- do.call(rbind, lapply(files, function(file) read.csv(gzfile(file.path(output, file)))))
if (anyDuplicated(data[, c("sample_id", "type", "variable")])) stop("Duplicate sample-feature entries")
data$ident <- paste0(data$type, "_0_", data$variable)
wide <- reshape2::dcast(data, sample_id ~ ident, value.var = "value")
rownames(wide) <- wide$sample_id
wide$sample_id <- NULL
wide <- wide[rowSums(!is.na(wide)) > 0, , drop = FALSE]

inverse_normal_rank <- function(x) {
  observed <- !is.na(x)
  x[observed] <- qnorm(rank(x[observed], ties.method = "average") / (sum(observed) + 1))
  x
}

normalized <- data.frame(lapply(wide, inverse_normal_rank), row.names = rownames(wide))
normalized$sample_id <- rownames(normalized)
data <- reshape2::melt(normalized, id.vars = "sample_id")
data$type <- sub("_0_.*", "", data$variable)
data$variable <- sub(".*_0_", "", data$variable)
write.csv(data, gzfile(file.path(output, "MOFA_Normalized_Features.csv.gz")), row.names = FALSE)

data <- data[!grepl("__RPS|__RPL|^RP", data$variable), ]
is_single_cell <- data$type == "single_cell"
data$type[is_single_cell] <- sub("(.*__).*$", "\\1", data$variable[is_single_cell])
data$variable <- gsub("γδ", "gd", gsub("Naïve", "Naive", gsub("Vδ", "VD", data$variable, fixed = TRUE), fixed = TRUE), fixed = TRUE)
data$type <- gsub("Naïve", "Naive", gsub("\\.", "x", data$type), fixed = TRUE)
samples <- sort(unique(data$sample_id))
views <- lapply(unique(data$type), function(view) {
  d <- reshape2::dcast(data[data$type == view, ], variable ~ sample_id, value.var = "value")
  rownames(d) <- d$variable
  d$variable <- NULL
  for (sample in setdiff(samples, colnames(d))) d[[sample]] <- NA_real_
  as.matrix(d[, samples, drop = FALSE])
})
names(views) <- gsub("__|\\.", "", unique(data$type))
saveRDS(views, file.path(output, "MOFA_View_Inputs.rds"))
summary <- data.frame(view = names(views), features = vapply(views, nrow, integer(1)),
                      samples = vapply(views, ncol, integer(1)))
write.csv(summary, file.path(output, "MOFA_View_Summary.csv"), row.names = FALSE)
