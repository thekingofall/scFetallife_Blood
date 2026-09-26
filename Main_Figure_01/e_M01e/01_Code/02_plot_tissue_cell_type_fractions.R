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
resultw <- read.delim(file.path(data_dir,"M01e_cell_type_tissue_fractions.tsv"))
resultw$Body <- factor(resultw$Body,levels=tissues)
cell_types <- unique(resultw$Last_cell_type_num)
resultw$Last_cell_type_num <- factor(resultw$Last_cell_type_num,levels=cell_types[order(as.numeric(sub("_.*","",cell_types)))])
resultw$Cell_lineage <- factor(resultw$Cell_lineage,levels=lineages)
P21=ggplot(resultw, aes(x= Last_cell_type_num, y=WeightedAverage,fill=Body)) +
  geom_col(position="fill") +
  theme_minimal() +
  xlab("Cell Type") +
  ylab("Weighted Average") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +

facet_grid(. ~ Cell_lineage, space = "free", scales = "free", switch = "y") +
  scale_fill_manual(values= c('#C71000FF',"#f49128","#023f75","#5A9599FF"))+
  theme_bw(base_size = 16, base_family="Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+xlab("")+ylab("Relative abundance (%)")  +
scale_y_continuous(labels=seq(0,100, by=25)) +xlab("")

ggsave(file.path(panel_dir,"02_Figures/M01e_tissue_cell_type_fractions.pdf"),plot=P21,width=16,height=8,device=cairo_pdf)
ggsave(file.path(panel_dir,"02_Figures/M01e_tissue_cell_type_fractions.png"),plot=P21,width=16,height=8,dpi=300,bg="white")
