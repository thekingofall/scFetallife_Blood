# Shared R library configuration and package-version checks.
.scf_setup_sources <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_setup_source <- if (length(.scf_setup_sources)) tail(.scf_setup_sources, 1)[[1]] else NULL
.scf_setup_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_setup_dir <- if (!is.null(.scf_setup_source)) {
  dirname(normalizePath(.scf_setup_source, mustWork = TRUE))
} else if (length(.scf_setup_args)) {
  dirname(normalizePath(sub("^--file=", "", .scf_setup_args[[1]]), mustWork = TRUE))
} else {
  getwd()
}
.scf_R_reference <- read.delim(
  file.path(.scf_setup_dir, "R_package_versions.tsv"),
  colClasses = "character", check.names = FALSE
)

scf_configure_R <- function(library_path = Sys.getenv(
  "SCF_R_LIB", unset = "/data5/Lipeng/mambaforge-pypy3/envs/S5/lib/R/library"
)) {
  if (nzchar(library_path) && dir.exists(library_path)) {
    .libPaths(unique(c(library_path, .libPaths())))
  }
  invisible(.libPaths())
}

scf_check_R_packages <- function(packages = NULL, strict = TRUE,
                                 reference = .scf_R_reference) {
  if (is.null(packages)) packages <- setdiff(reference$package, "R")
  packages <- unique(c("R", as.character(packages)))
  installed <- vapply(packages, function(package) {
    if (package == "R") return(as.character(getRversion()))
    if (package %in% loadedNamespaces()) {
      return(as.character(getNamespaceVersion(package)))
    }
    if (!length(find.package(package, quiet = TRUE))) return(NA_character_)
    tryCatch(as.character(utils::packageVersion(package)),
             error = function(e) NA_character_)
  }, character(1))
  expected <- reference$version[match(packages, reference$package)]
  versions_match <- mapply(function(actual, recorded) {
    if (is.na(actual) || is.na(recorded)) return(FALSE)
    utils::compareVersion(actual, recorded) == 0L
  }, installed, expected)
  status <- ifelse(is.na(installed), "missing",
    ifelse(is.na(expected), "not_recorded",
      ifelse(versions_match, "match", "different_version")))
  result <- data.frame(package = packages, installed = unname(installed),
                       reference = expected, status = status,
                       stringsAsFactors = FALSE)
  different <- result$status == "different_version"
  if (any(different)) {
    message("Versions differ from the configured reference versions: ", paste(
      paste0(result$package[different], " ", result$installed[different],
             " (reference ", result$reference[different], ")"), collapse = "; "))
  }
  missing <- result$package[result$status == "missing"]
  if (length(missing) && strict) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Install these packages before running this script.", call. = FALSE)
  }
  invisible(result)
}

if (!interactive()) {
  options(device = function(file = NULL, ...) grDevices::pdf(file = file, ...))
}
scf_configure_R()
if (sys.nframe() == 0L) {
  print(scf_check_R_packages(strict = FALSE), row.names = FALSE)
}

.scf_project_root <- normalizePath(file.path(.scf_setup_dir, ".."), winslash = "/")
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")
