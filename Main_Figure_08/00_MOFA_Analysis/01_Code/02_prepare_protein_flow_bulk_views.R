.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
root <- normalizePath(file.path(.scf_start_dir, "../../.."), winslash = "/")
source(file.path(root, "00_Settings/00_setup_R.R"))
scf_check_R_packages("readxl")
global <- file.path(root, "00_Global_Data")
output <- file.path(.scf_start_dir, "../02_Data")
metadata <- read.csv(file.path(output, "MOFA_Sample_Metadata.csv"))
samples <- metadata$MainID

as_long <- function(x, view) {
  x <- x[, samples, drop = FALSE]
  data.frame(sample_id = rep(colnames(x), each = nrow(x)),
             variable = rep(rownames(x), times = ncol(x)), value = as.vector(x), type = view)
}

read_olink <- function(sheet, view) {
  raw <- as.data.frame(readxl::read_excel(file.path(global, "Fetal_Immune_Atlas_Olink_NPX.xlsx"),
                                        sheet = sheet, col_names = FALSE, col_types = "text"))
  protein_columns <- which(grepl("^OID", as.character(raw[6, ])))
  sample_rows <- match(samples, raw[[1]])
  if (anyNA(sample_rows)) stop("Missing Olink samples in sheet: ", sheet)
  x <- matrix(as.numeric(as.matrix(raw[sample_rows, protein_columns])),
              nrow = length(samples), dimnames = list(samples, as.character(raw[4, protein_columns])))
  # Round NPX to two decimals, with half values rounded away from zero.
  x <- sign(x) * floor(abs(x) * 100 + 0.5 + 1e-10) / 100
  as_long(t(x), view)
}

plasma <- read_olink("plasma", "Plama_Olink")
stimulated <- read_olink("cell supernatant", "Cell_stimulate_Olink")
flow <- read.csv(file.path(global, "Shared_Inputs/Main_Figure_08/00_MOFA/01_Source_Data/Flow_Parent_Frequencies.csv"),
                 row.names = 1, check.names = FALSE)
flow <- as_long(as.matrix(flow), "Sflow_Freq_inParent")
bulk <- read.csv(file.path(global, "Fetal_Immune_Atlas_PHA_Bulk_TPM.csv"), row.names = 1, check.names = FALSE)
bulk <- as.matrix(bulk[, samples, drop = FALSE])
bulk <- bulk[rowSums(bulk) > 10, , drop = FALSE]
bulk <- log2(bulk + 1)
bulk <- bulk[rowSums(bulk) > 100, , drop = FALSE]
bulk <- as_long(bulk, "Cell_stimulate_BulkRNA")
write.csv(rbind(plasma, stimulated, flow, bulk),
          gzfile(file.path(output, "Protein_Flow_Bulk_Features.csv.gz")), row.names = FALSE)
