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
counts <- read.table(file.path(.scf_shared,"Supplementary_Figures/Supplementary_Figure_17/01_Bulk_Expression/CellAll_count.txt"),header=TRUE)
counts$Geneid <- substr(counts$Geneid,1,15)
colnames(counts) <- gsub(".sorted.bam","",colnames(counts))
colnames(counts) <- gsub("Cell_A03.Alignment.","Z",colnames(counts))
for (prefix in c("D","B","C")) colnames(counts) <- gsub(paste0("Z",prefix),prefix,colnames(counts))
gene_map <- read.csv(file.path(.scf_shared,"00_Common/02_Gene_Annotation/GRCh38_108_gene_names.csv"))
merged <- merge(counts,gene_map,by.x="Geneid",by.y="gene_id")
merged <- merged[!is.na(merged$gene_name),]
merged <- merged[!duplicated(merged$gene_name),]
sample_names <- c("Z48","Z67","Z99","Z3","Z79","Z43","Z53","Z47","Z7","Z38","Z98","Z33","Z70","Z52","Z50","Z96")
sample_meta <- read.csv(file.path(.scf_global,"Fetal_Immune_Atlas_Sample_Metadata.csv"))
sample_key <- function(x) sub("^([A-Z])0+([0-9])", "\\1\\2", toupper(x))
sample_meta <- sample_meta[match(sample_key(sample_names),sample_key(sample_meta$Name)),]
stopifnot(!anyNA(sample_meta$Name),!anyNA(sample_meta$Post_Conception_Age_Weeks))
expdata <- as.matrix(merged[,sample_names]);rownames(expdata)<-merged$gene_name
countToTpm <- function(counts,effLen) {
  rate <- log(counts)-log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate-denom+log(1e6))
}
tpm <- apply(expdata,2,function(x)countToTpm(x,merged$Length))
colnames(tpm) <- sample_meta$MainID
write.csv(tpm, file.path(.scf_global,"Fetal_Immune_Atlas_PHA_Bulk_TPM.csv"))
