.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("MOFA2"))
suppressPackageStartupMessages(library(MOFA2))
root <- normalizePath(file.path(.scf_start_dir, "../../.."), winslash = "/")
data_list <- readRDS(file.path(.scf_start_dir, "../02_Data/MOFA_View_Inputs.rds"))
model <- create_mofa(data_list)
data_options <- get_default_data_options(model)
data_options$scale_views <- TRUE
model_options <- get_default_model_options(model)
model_options$num_factors <- 16
training_options <- get_default_training_options(model)
training_options$seed <- 2024
training_options$maxiter <- 50000
training_options$weight_views <- FALSE
model <- prepare_mofa(model, data_options = data_options, model_options = model_options,
                      training_options = training_options)
output_dir <- file.path(root, "Main_Figure_08/00_MOFA_Analysis/02_Data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
run_mofa(model, outfile = file.path(output_dir, "MOFA_Model.hdf5"))
