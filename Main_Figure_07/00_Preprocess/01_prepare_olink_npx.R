args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
script <- if (length(files)) tail(files, 1)[[1]] else sub("^--file=", "", args[[1]])
here <- dirname(normalizePath(script))
root <- normalizePath(file.path(here, "../.."))
source(file.path(root, "00_Settings/00_setup_R.R"))

scf_check_R_packages("readxl")
output <- file.path(here, "02_Processed")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
metadata <- read.csv(file.path(root, "00_Global_Data/Fetal_Immune_Atlas_Sample_Metadata.csv"), check.names = FALSE)
workbook <- file.path(root, "00_Global_Data/Fetal_Immune_Atlas_Olink_NPX.xlsx")
for (sheet in c("plasma", "cell supernatant")) {
  raw <- as.data.frame(readxl::read_excel(workbook, sheet = sheet, col_names = FALSE, col_types = "text"))
  columns <- which(grepl("^OID", as.character(raw[6, ])))
  available <- if (sheet == "plasma") "Plasma_Olink_Available" else "Stimulated_Olink_Available"
  samples <- metadata[metadata$Main_Organ == "PBMC" & metadata[[available]] == 1, ]
  rows <- match(samples$MainID, raw[[1]])
  if (anyNA(rows)) stop("Missing samples in Olink sheet: ", sheet)
  x <- matrix(as.numeric(as.matrix(raw[rows, columns])), nrow = length(rows))
  x <- sign(x) * floor(abs(x) * 100 + 0.5 + 1e-10) / 100
  result <- data.frame(MainID = rep(samples$MainID, times = length(columns)),
    Post_Conception_Age_Weeks = rep(samples$Post_Conception_Age_Weeks, times = length(columns)),
    Protein = rep(as.character(raw[4, columns]), each = length(rows)), NPX = as.vector(x))
  name <- if (sheet == "plasma") "M07_Plasma_NPX.csv" else "M07_Stimulated_NPX.csv"
  write.csv(result, file.path(output, name), row.names = FALSE)
}
