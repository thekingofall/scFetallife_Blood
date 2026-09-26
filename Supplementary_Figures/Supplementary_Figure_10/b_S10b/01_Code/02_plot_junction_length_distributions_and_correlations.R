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
chain <- "TRB"
regions <- c("3'V-REGION","N1-REGION","D-REGION","N2-REGION","5'J-REGION")
create_pval_plot <- function(df,titlex) {
 
  df$pvalue_cat <- cut(df$p_value,
                       breaks = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                       
                       labels = c('< 0.0001', '0.0001-0.001', '0.001-0.01',
                                  '0.01-0.05', '> 0.05'),
                       include.lowest = TRUE)

  colors <- c('< 0.0001' = "#C71000B2", '0.0001-0.001' =  "#FF6348B2", 
              '0.001-0.01' = "#FF95A8B2", '0.01-0.05' = "#8A4198B2" ,
              '> 0.05' ="#008EA0B2" )

 
  plot <- ggplot(df, aes(x = factor(lengths), y = R_value, fill = pvalue_cat)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = colors) +
      labs(x = "Lengths", y = "R Value(spearman) ", fill = "P Value") +
      theme_bw(base_family="Arial") +
      ggtitle(titlex) + 
      theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           panel.border = element_rect(linetype = "solid", colour = "black", size = 1.5)) + 
      scale_y_continuous(breaks = seq(-1, 1, by = 0.2), labels = seq(-1, 1, by = 0.2))
  
  return(plot)
}



data <- read.csv(file.path(data_dir,paste0(chain,"_junction_length_distributions.csv")))
statistics <- read.csv(file.path(data_dir,paste0(chain,"_junction_length_age_statistics.csv")))
plots <- list()
for (region in regions) {
  z <- data[data$region==region,]
  distribution <- ggplot(z)+
    geom_point(aes(Ntlen,Freq_ratio,color=Post_Conception_Age_Weeks),shape=1)+
    geom_smooth(aes(Ntlen,Freq_ratio,group=Post_Conception_Age_Weeks,color=Post_Conception_Age_Weeks),se=FALSE,linewidth=0.3)+
    theme_linedraw(base_family="Arial")+
    labs(x="Length (nt)",y="ratio (%)",color="pcw",title=region)+
    theme(plot.title=element_text(hjust=0.5,size=13,face="bold"),panel.border=element_rect(color="black",linewidth=1.5))+
    scale_color_gradientn(colours=rev(colorRampPalette(c("#C71000B2","#FF6F00B2","#6a73cf","#00AF99"))(100)))
  correlation <- create_pval_plot(statistics[statistics$region==region,],region)+labs(x="Length (nt)",y="Spearman r")
  plots <- c(plots,list(distribution,correlation))
}
combined <- arrangeGrob(grobs=plots,ncol=2)
stem <- paste0(chain,"_junction_length_distributions_and_age_associations")
ggsave(file.path(panel_dir,"02_Figures",paste0(stem,".pdf")),combined,width=12,height=3.3*length(regions),device=cairo_pdf)
ggsave(file.path(panel_dir,"02_Figures",paste0(stem,".png")),combined,width=12,height=3.3*length(regions),dpi=300,bg="white")
