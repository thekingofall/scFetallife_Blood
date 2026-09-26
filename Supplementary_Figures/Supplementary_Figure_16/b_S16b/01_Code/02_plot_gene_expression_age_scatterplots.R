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

library(ggplot2)
panel_dir <- normalizePath(file.path(.scf_start_dir,".."),winslash="/")
data_dir <- file.path(panel_dir,"03_PlotData")
result <- read.delim(file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_16/a_S16a/03_PlotData/S16_spearman_results.tsv"))
long <- read.delim(file.path(.scf_shared, "Supplementary_Figures/Supplementary_Figure_16/a_S16a/03_PlotData/S16_log2TPM_long.tsv"))
original <- read.csv(file.path(panel_dir,"03_Source_Data/PBMC_tpm_olinkRP_original_submission.csv"))
result <- result[order(-result$spearman_r),]
colorname2 <- c(
  "#46A040", "#00AF99", "#FFC179", "#98D9E9", "#F6313E", "#FFA300",
  "#333366", "#FF5A00", "#663366", "#FF6666", "#8F1336", "#0081C9",
  "#001588", "#CC0033", "#CC9966", "#CC0033", "#999933", "#009966",
  "#CCCC33", "#CCFF99", "#333399", "#993333", "#490C65", "#BA7FD0",
  "#A6CEE3", "#1F78B4", "#DE77AE", "#B2DF8A", "#006D2C", "#868686",
  "#B5AD64", "#9DA8E2", "#91C392", "#FF9900", "#339966", "#46A040",
  "#00AF99", "#FFC179", "#98D9E9", "#F6313E", "#FFA300", "#333366",
  "#FF5A00", "#663366", "#FF6666", "#8F1336", "#0081C9", "#001588",
  "#CC0033", "#CC9966", "#CC0033", "#999933", "#009966", "#CCCC33",
  "#CCFF99", "#333399", "#993333"
)
colorname3 <- rep(colorname2, 5)


original <- original[!is.na(original$Assay) & nzchar(original$Assay), , drop = FALSE]
original_assay_levels <- sort(unique(original$Assay))
stopifnot(length(original_assay_levels) <= length(colorname3))
original_assay_colors <- setNames(
  colorname3[seq_along(original_assay_levels)],
  original_assay_levels
)

gene_to_original_assay <- tapply(
  original$Assay,
  original$OlinkGene,
  function(x) x[[1]]
)
result$original_assay <- unname(gene_to_original_assay[result$gene_symbol])
result$plot_color <- unname(original_assay_colors[result$original_assay])

if (anyNA(result$plot_color)) {
  missing_genes <- result$gene_symbol[is.na(result$plot_color)]
  stop("Missing gene colors for: ", paste(missing_genes, collapse = ", "))
}

long$assay <- factor(long$assay,levels=result$assay)
result$assay <- factor(result$assay,levels=levels(long$assay))
P1 <- ggplot(long,aes(x=Post_Conception_Age_Weeks,y=log2_TPM_plus_1,color=assay))+
  geom_point()+geom_smooth(method="lm",se=TRUE,formula=y~x)+
  facet_wrap(~assay,scales="free_y",ncol=7)+
  scale_color_manual(values=setNames(result$plot_color,as.character(result$assay)))+
  scale_x_continuous(breaks=c(16,28,40),limits=c(16,40))+
  geom_text(data=result,aes(x=Inf,y=Inf,label=paste("r = ",round(spearman_r,2),"\np = ",round(p_value,3),"\nq = ",round(q_value_BH,3))),hjust=1,vjust=1)+
  labs(x="pcw",y="Log2(TPM + 1)")+theme_bw(base_family="Arial")+
  theme(legend.position="none",panel.grid=element_blank(),strip.background=element_rect(fill="white",color="black"),aspect.ratio=0.8)
ggsave(file.path(panel_dir,"02_Figures/S16b_gene_expression_age_scatterplots.pdf"),P1,width=14,height=18,device=cairo_pdf)
ggsave(file.path(panel_dir,"02_Figures/S16b_gene_expression_age_scatterplots.png"),P1,width=14,height=18,dpi=250,bg="white")
