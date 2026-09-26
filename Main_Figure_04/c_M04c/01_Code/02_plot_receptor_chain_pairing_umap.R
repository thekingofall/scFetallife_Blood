.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("ggplot2", "ggrastr"))
library(ggplot2)
library(ggrastr)
panel_dir <- normalizePath(file.path(.scf_start_dir,".."))
receptor <- "TCR"
data <- read.csv(file.path(panel_dir,"03_PlotData",paste0(receptor,"_chain_pairing_UMAP.csv.gz")))
colors <- c("ambiguous"="#8A252B", "extra VDJ"="#8A4198", "extra VJ"="#145096",
            "orphan VDJ"="#F6C619", "orphan VJ"="#E9412F", "single pair"="#008EA0", "two full chains"="#266B69")
data$chain_pairing <- factor(data$chain_pairing,levels=names(colors))
data <- data[order(data$chain_pairing),]
p <- ggplot(data,aes(UMAP1,UMAP2,color=chain_pairing))+
  geom_point_rast(size=0.10,stroke=0,alpha=0.88,raster.dpi=300)+
  scale_color_manual(values=colors)+
  coord_fixed()+labs(x="UMAP1",y="UMAP2")+
  theme_classic(base_size=10,base_family="Arial")+
  theme(legend.position="none",axis.text=element_blank(),axis.ticks=element_blank(),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.7),axis.line=element_blank())
out <- file.path(panel_dir,"02_Figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(out,paste0(receptor,"_chain_pairing_UMAP.pdf")),p,width=2.3,height=2.3,device=cairo_pdf)
ggsave(file.path(out,paste0(receptor,"_chain_pairing_UMAP.png")),p,width=2.3,height=2.3,dpi=300,bg="white")
