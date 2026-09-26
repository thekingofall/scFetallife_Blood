.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "tidyr", "tibble", "stringr", "Biobase", "Mfuzz", "e1071", "purrr", "reshape2", "ComplexHeatmap"))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(tibble);library(stringr);library(Biobase);library(Mfuzz);library(purrr);library(reshape2);library(ComplexHeatmap)})
root <- normalizePath(file.path(.scf_start_dir, "../../../.."), winslash = "/")
panel <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
plot_data <- file.path(panel, "03_PlotData")
dir.create(plot_data, showWarnings = FALSE)
input <- readRDS(file.path(root, "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_17/00_Counts_and_Gene_Lengths/PBMC_Counts_and_Gene_Lengths.rds"))
tpm <- as.data.frame(apply(input$counts, 2, function(x) {
  rate <- log(x) - log(input$gene_length)
  exp(rate - log(sum(exp(rate))) + log(1e6))
}))
metadata <- data.frame(MainID = colnames(tpm))
metadata$Post_Conception_Age_Weeks <- as.numeric(str_extract(metadata$MainID, "[0-9]+[.][0-9]+"))
metadata$Grouping_Age_Weeks <- round(metadata$Post_Conception_Age_Weeks)
metadata$stage <- case_when(metadata$Grouping_Age_Weeks < 22 ~ "16-22PCW",
                            metadata$Grouping_Age_Weeks < 28 ~ "22-28PCW",
                            metadata$Grouping_Age_Weeks < 34 ~ "28-34PCW",
                            metadata$Grouping_Age_Weeks <= 40 ~ "34-40PCW")
long <- tpm %>% rownames_to_column("gene") %>% gather("MainID", "value", -gene) %>%
  left_join(metadata, by = "MainID")
means <- long %>% group_by(gene, stage) %>% summarise(mean_value = mean(value), .groups = "drop")
wide <- as.data.frame(spread(means, key = stage, value = mean_value))
rownames(wide) <- wide$gene
wide <- wide[, -1]
wide <- wide[rowSums(wide) > 0, ]
source(file.path(root, "Main_Figure_02/e_M02e/01_Code/01_define_expression_clustering.R"))
cm <- clusterData(exp = wide, cluster.method = "mfuzz", cluster.num = 16)
saveRDS(cm, file.path(plot_data, "PBMC_Bulk_16_Modules.rds"))
saveRDS(wide, file.path(plot_data, "PBMC_Bulk_Stage_Mean_TPM.rds"))
write.csv(wide, file.path(plot_data, "PBMC_Bulk_Stage_Mean_TPM.csv"))
write.csv(metadata, file.path(plot_data, "PBMC_Bulk_Stage_Assignments.csv"), row.names = FALSE)
write.csv(as.data.frame(table(cm$wide.res$cluster)), file.path(plot_data, "PBMC_Bulk_Module_Gene_Counts.csv"), row.names = FALSE)
