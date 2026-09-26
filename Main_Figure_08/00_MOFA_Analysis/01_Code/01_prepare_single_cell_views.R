.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
root <- normalizePath(file.path(.scf_start_dir, "../../.."), winslash = "/")
source(file.path(root, "00_Settings/00_setup_R.R"))
input <- file.path(root, "00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA/01_Source_Data")
output <- file.path(.scf_start_dir, "../02_Data")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(file.path(root, "00_Global_Data/Fetal_Immune_Atlas_Sample_Metadata.csv"))
assays <- c("scRNAseq_Available", "Flow_Cytometry_Available", "Plasma_Olink_Available",
            "Stimulated_Olink_Available", "PHA_Bulk_RNA_Available")
metadata <- metadata[metadata$Main_Organ == "PBMC" & rowSums(metadata[, assays] == 1, na.rm = TRUE) == length(assays), ]
metadata <- metadata[order(metadata$Post_Conception_Age_Weeks, metadata$MainID), ]
samples <- metadata$MainID
write.csv(metadata, file.path(output, "MOFA_Sample_Metadata.csv"), row.names = FALSE)

expression <- read.csv(gzfile(file.path(input, "Single_Cell_Pseudobulk.csv.gz")), row.names = 1, check.names = FALSE)
expression <- t(as.matrix(expression))
selection <- read.csv(file.path(input, "Single_Cell_Selected_Genes.csv"), check.names = FALSE)
sample_ids <- sub("^([^_]+_[^_]+)_.*", "\\1", colnames(expression))
cell_types <- sub("^[^_]+_[^_]+_", "", colnames(expression))

quantile_normalize <- function(x) {
  ranks <- apply(x, 2, rank, ties.method = "min")
  means <- rowMeans(apply(x, 2, sort))
  result <- apply(ranks, 2, function(rank) means[rank])
  dimnames(result) <- dimnames(x)
  result
}

features <- lapply(unique(selection$cell_type), function(cell_type) {
  keep <- cell_types == cell_type & sample_ids %in% samples
  x <- expression[, keep, drop = FALSE]
  colnames(x) <- sample_ids[keep]
  library_scale <- colSums(x) / mean(colSums(x))
  library_scale[library_scale == 0] <- 1
  x <- sweep(x, 2, library_scale, "/")
  genes <- selection$gene[selection$cell_type == cell_type]
  x <- x[rownames(x) %in% genes, , drop = FALSE]
  if (!nrow(x) || ncol(x) < 2) stop("Insufficient expression data for ", cell_type)
  x <- quantile_normalize(log2(x + 1))
  aligned <- matrix(0, nrow(x), length(samples), dimnames = list(rownames(x), samples))
  aligned[, colnames(x)] <- x
  label <- gsub(" ", "_", gsub("CXCR5\\+", "CXCR5high", gsub("CXCR5-", "CXCR5low", cell_type)))
  data.frame(sample_id = rep(samples, each = nrow(aligned)),
             variable = rep(make.names(paste0(label, "__", rownames(aligned))), times = length(samples)),
             value = as.vector(aligned), type = "single_cell")
})
write.csv(do.call(rbind, features), gzfile(file.path(output, "Single_Cell_Features.csv.gz")), row.names = FALSE)
