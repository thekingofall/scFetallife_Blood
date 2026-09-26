.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("dplyr", "ggplot2"))

.scf_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.scf_code_dir <- .scf_start_dir
dir.create(file.path(.scf_code_dir, "../02_Figures"), recursive = TRUE, showWarnings = FALSE)
.scf_project_root <- normalizePath(file.path(.scf_code_dir, "../../.."), winslash = "/", mustWork = FALSE)
.scf_global <- file.path(.scf_project_root, "00_Global_Data")
.scf_shared <- file.path(.scf_global, "Shared_Inputs")

library(dplyr)
library(ggplot2)
panel_dir <- normalizePath(file.path(.scf_start_dir, ".."), winslash="/")
data_dir <- file.path(panel_dir, "03_PlotData")
dir.create(data_dir, recursive=TRUE, showWarnings=FALSE)
lineages <- c("PRECURSOR","B_CELL","T/ILC","NK","MYELOID","DC","MK/ERY","OTHERS")
tissues <- c("PBMC","Liver","Thymus","Spleen")
adata1_obs <- read.csv(file.path(.scf_shared,"Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv"),check.names=FALSE)
cell_types <- unique(adata1_obs$Last_cell_type_num)
cell_types <- cell_types[order(as.numeric(sub("_.*", "", cell_types)))]
grouped_data_result <- adata1_obs %>%
  group_by(MainID,Last_cell_type_num,Cell_lineage) %>%
  summarise(count=n(),.groups="drop") %>%
  group_by(MainID) %>% mutate(week_total=sum(count),percentage=count/week_total) %>%
  ungroup() %>% mutate(Body=substr(MainID,1,1))
body_count <- grouped_data_result %>% distinct(MainID,Body) %>% count(Body,name="body_count")
grouped_data_result <- left_join(grouped_data_result,body_count,by="Body")
resultw <- grouped_data_result %>% group_by(Body,Last_cell_type_num,Cell_lineage) %>%
  summarise(WeightedAverage=weighted.mean(percentage,week_total*body_count),.groups="drop")
resultw$Last_cell_type_num <- factor(resultw$Last_cell_type_num,levels=cell_types)
resultw$Cell_lineage <- factor(resultw$Cell_lineage,levels=lineages)
resultw$Body <- factor(resultw$Body,levels=c("B","L","T","S"),labels=tissues)
write.table(resultw,file.path(data_dir,"M01e_cell_type_tissue_fractions.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
