args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
script <- if (length(files)) tail(files, 1)[[1]] else sub("^--file=", "", args[[1]])
here <- dirname(normalizePath(script))
root <- normalizePath(file.path(here, "../.."))
source(file.path(root, "00_Settings/00_setup_R.R"))

scf_check_R_packages(c("Matrix", "SeuratObject"))
data_dir <- file.path(here, "02_Processed")
metadata <- read.csv(file.path(data_dir, "M03_Cell_Metadata.csv"), row.names = 1, check.names = FALSE)
genes <- readLines(file.path(data_dir, "M03_Genes.tsv"))
read_matrix <- function(name) {
  connection <- gzfile(file.path(data_dir, paste0("M03_", name, ".mtx.gz")))
  on.exit(close(connection))
  x <- as(Matrix::readMM(connection), "CsparseMatrix")
  dimnames(x) <- list(genes, rownames(metadata))
  x
}
counts <- read_matrix("counts")
object <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = metadata)
expression <- read_matrix("log_normalized")
rownames(expression) <- rownames(object)
object <- SeuratObject::SetAssayData(object, slot = "data", new.data = expression)
saveRDS(object, file.path(root, "00_Global_Data/Fetal_Immune_Atlas_Seurat.rds"))
