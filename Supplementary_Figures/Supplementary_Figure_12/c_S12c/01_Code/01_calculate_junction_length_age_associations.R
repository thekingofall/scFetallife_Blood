.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "gridExtra", "dplyr", "tidyr", "purrr", "broom"))

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
panel_dir <- normalizePath(file.path(.scf_start_dir,".."),winslash="/")
data_dir <- file.path(panel_dir,"03_PlotData")
dir.create(data_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(panel_dir,"02_Figures"),recursive=TRUE,showWarnings=FALSE)
chain <- "IGL"
regions <- c("3'V-REGION","N-REGION","5'J-REGION")
Dfcalculate_correlation <- function(df, group_var, var1, var2) {
  results_df <- df %>%
    group_by(.data[[group_var]]) %>%
    nest() %>%
    mutate(correlation = map(data, ~cor.test(.x[[var1]], .x[[var2]], method = "spearman")) %>% map(broom::tidy)) %>%
    unnest(correlation) %>%
    select(.data[[group_var]], estimate, p.value)
colnames(results_df)<-c('lengths','R_value', 'p_value')
  return(results_df)
}

input_files <- c("IGL_3V_REGION_length_distribution.csv","IGL_N_REGION_length_distribution.csv","IGL_5J_REGION_length_distribution.csv")
input_dir <- file.path(.scf_shared,"Supplementary_Figures/Supplementary_Figure_12/01_Junction_Length")

distributions <- list();statistics <- list()
for (i in seq_along(regions)) {
  data <- read.csv(file.path(input_dir,input_files[i]),check.names=FALSE)
  data$Post_Conception_Age_Weeks <- as.numeric(substring(data$MainID,2,5))
  data$Ntlen <- data$length_nt
  data$Freq_ratio <- data$frequency_percent
  result <- suppressWarnings(Dfcalculate_correlation(data,"Ntlen","Post_Conception_Age_Weeks","Freq_ratio"))
  result$region <- regions[i]
  distributions[[i]] <- data
  statistics[[i]] <- result
}
write.csv(do.call(rbind,distributions),file.path(data_dir,paste0(chain,"_junction_length_distributions.csv")),row.names=FALSE)
write.csv(do.call(rbind,statistics),file.path(data_dir,paste0(chain,"_junction_length_age_statistics.csv")),row.names=FALSE)
