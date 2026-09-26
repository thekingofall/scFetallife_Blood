.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

options(stringsAsFactors=FALSE)
tpm <- as.matrix(read.csv(file.path(.scf_global,"Fetal_Immune_Atlas_PHA_Bulk_TPM.csv"),row.names=1,check.names=FALSE))
sample_meta <- read.csv(file.path(.scf_global,"Fetal_Immune_Atlas_Sample_Metadata.csv"))
sample_meta <- sample_meta[match(colnames(tpm),sample_meta$MainID),]
stopifnot(!anyNA(sample_meta$MainID),!anyNA(sample_meta$Post_Conception_Age_Weeks))
map <- read.delim(file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_16/a_S16a/03_PlotData/S16_assay_gene_mapping.tsv"))
map$assay[map$assay=="STAMPB"] <- "STAMBP"
log_tpm <- log2(tpm[map$gene_symbol,,drop=FALSE]+1)
colnames(log_tpm) <- sample_meta$MainID
result <- map
tests <- lapply(seq_len(nrow(log_tpm)),function(i)
  suppressWarnings(cor.test(seq_len(ncol(log_tpm)),log_tpm[i,],method="spearman")))
result$spearman_r <- vapply(tests,function(x)unname(x$estimate),numeric(1))
result$p_value <- vapply(tests,function(x)x$p.value,numeric(1))
result$q_value_BH <- p.adjust(result$p_value,"BH")
result <- result[order(-result$spearman_r),]
result$rank_by_r <- seq_len(nrow(result))
long <- do.call(rbind,lapply(seq_len(nrow(map)),function(i)data.frame(
  assay=map$assay[i],gene_symbol=map$gene_symbol[i],MainID=sample_meta$MainID,
  Gestational_Age_Weeks=sample_meta$Gestational_Age_Weeks,
  Post_Conception_Age_Weeks=sample_meta$Post_Conception_Age_Weeks,
  sample_order=seq_len(nrow(sample_meta)),log2_TPM_plus_1=as.numeric(log_tpm[i,]))))
write.table(result,file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_16/a_S16a/03_PlotData/S16_spearman_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(long,file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_16/a_S16a/03_PlotData/S16_log2TPM_long.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(sample_meta,file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_16/a_S16a/03_PlotData/S16_sample_metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
