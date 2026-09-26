.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("MOFA2", "dplyr", "clusterProfiler", "org.Hs.eg.db"))
suppressPackageStartupMessages({library(MOFA2);library(dplyr);library(clusterProfiler);library(org.Hs.eg.db)})
root <- normalizePath(file.path(.scf_start_dir, "../../.."), winslash = "/")
input <- file.path(root, "00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA")
output <- file.path(.scf_start_dir, "../02_Data")
dir.create(output, showWarnings = FALSE)
model <- load_model(file.path(input, "Fetal_Immune_Atlas_MOFA.hdf5"), verbose = FALSE)
views <- names(readRDS(file.path(input, "MOFA_View_Inputs.rds")))[c(4, 8:25)]
weights <- get_weights(model, factors = 3, views = views, as.data.frame = TRUE)
weights$sign <- ifelse(weights$value > 0, "positive", "negative")
weights <- weights[weights$value != 0, ]
top <- weights %>% group_by(view, sign) %>% top_n(100, wt = abs(value)) %>% ungroup()
top$gene <- gsub("_Cell_stimulate_BulkRNA", "", sub(".*__", "", top$feature))
top$Group <- paste(top$sign, top$view, sep = "_")
results <- lapply(unique(top$Group), function(group) {
  genes <- unique(top$gene[top$Group == group])
  enriched <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = "BH", minGSSize = 1,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  d <- as.data.frame(enriched)
  if (!nrow(d)) return(NULL)
  d$Group <- group
  d
})
all_terms <- bind_rows(results)
selected <- all_terms %>% filter(!grepl("cytoplasmic translation", Description)) %>%
  group_by(Group) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 1) %>% ungroup()
write.csv(top, file.path(output, "MOFA_Factor3_Top100_Features.csv"), row.names = FALSE)
write.csv(all_terms, file.path(output, "MOFA_Factor3_GO_BP.csv"), row.names = FALSE)
write.csv(selected, file.path(output, "MOFA_Factor3_Selected_GO_BP.csv"), row.names = FALSE)
