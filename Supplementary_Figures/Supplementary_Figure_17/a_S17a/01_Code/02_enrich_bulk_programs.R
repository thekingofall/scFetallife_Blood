.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("clusterProfiler", "org.Hs.eg.db"))
suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db)})
root <- normalizePath(file.path(.scf_start_dir, "../../../.."), winslash = "/")
panel <- normalizePath(file.path(.scf_start_dir, ".."), winslash = "/")
plot_data <- file.path(panel, "03_PlotData")
dir.create(plot_data, showWarnings = FALSE)
cm <- readRDS(file.path(plot_data, "PBMC_Bulk_16_Modules.rds"))
programs <- read.delim(file.path(plot_data, "PBMC_Bulk_Temporal_Programs.tsv"))
results <- lapply(seq_len(nrow(programs)), function(i) {
  genes <- unique(cm$wide.res$gene[cm$wide.res$cluster == programs$cluster[i]])
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  enrichment <- enrichGO(gene = unique(ids$ENTREZID), OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                         pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  d <- as.data.frame(enrichment)
  d$program <- programs$program[i]
  d$cluster <- programs$cluster[i]
  d
})
write.csv(do.call(rbind, results), file.path(plot_data, "PBMC_Bulk_Program_GO_BP.csv"), row.names = FALSE)
