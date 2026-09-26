.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("RColorBrewer"))

col1=c("#C71000FF","#8A4198FF",  "#023f75", "#008EA0FF",
 "#266b69","#f6c619", "#FF6348FF","#c29f62")


col2=c("#8A4198FF",  "#145096", "#008EA0FF", "#e9412f","#266b69","#f6c619", "#3178a2","#6a73cf")
col21=c("#3D3B25FF","#8A4198FF",  "#145096", "#008EA0FF", "#e9412f","#266b69","#f6c619", "#3178a2","#6a73cf")

col3=c("#C71000FF", "#008EA0FF", "#8A4198FF", "#5A9599FF",
"#FF6348FF", "#84D7E1FF", "#FF95A8FF", "#3D3B25FF", "#ADE2D0FF",
"#1A5354FF", "#3F4041FF")
colorgan=c('#C71000FF',"#f49128","#023f75","#5A9599FF")


colname4=c( "#C71000FF", "#008EA0FF", "#8A4198FF", "#5A9599FF",
"#FF6348FF", "#84D7E1FF", "#FF95A8FF", "#3D3B25FF", "#ADE2D0FF",
"#1A5354FF", "#3F4041FF","#fa6e01","#972b1d", "#e6a84b", "#4c211b", "#ff717f", "#009966",
"#c62d17", "#023f75", "#ea894e", "#266b69", "#eb4601", "#f6c619",
"#f49128", "#194a55", "#c29f62", "#83ba9e","#187c65", "#A6CEE3",
"#223e9c", "#aebea6", "#edae11", "#c74732",
"#6a73cf", "#edd064", "#0eb0c8", "#f2ccac", "#a1d5b9", "#e1abbc",
"#46A040", "#00AF99", "#FFC179", "#98D9E9", "#F6313E", "#FFA300",
"#333366", "#663366", "#FF6666", "#8F1336", "#0081C9",  "#CC0033",
"#CC9966", "#CC0033", "#999933", "#CCCC33", "#CCFF99", "#333399","#001588",
 "#490C65", "#BA7FD0", "#1F78B4", "#DE77AE" , "#2f2f2f",
"#006D2C", "#868686", "#9DA8E2", "#91C392", "#FF9900", "#339966","#993333")


color5=c("#8dd3c7","#00ffff","#00b1bf","#ff1593","#80b1d3","#b3de69","#e78bc4","#bbbbbb","#e62722","#a3cce2",
        "#fde8d2","#fa9274","#e9412f","#f5fad0","#fec48b","#edf8b3","#c1e6bc","#c02750","#111a3d","#b76235",
        "#145096")
col1=c("#8A4198FF",  "#145096","#f6c619","#e9412f", "#008EA0FF",
 "#266b69", "#3178a2")

col11=c("#8a252b","#8A4198FF",  "#145096","#f6c619","#e9412f", "#008EA0FF",
 "#266b69", "#3178a2")
seqcolor2_byr=c("#1a2d69","#2685c9","#96c52c","#d2dd36","#f4c91e","#e73718")
seqcolor3_wr=c("#fef3ed","#fcb398","#fb7858","#e84433","#b11518")
seqcolor4_bwr=c("#2166ac","#518bbd","#fef3ed","#d25e55","#b2182b")
seqcolor5_rwb=c('#C6403D', "#f6c619","white", '#246BAE')

cols01 <- c("#f49128", "#194a55", "#187c65", "#f26115", "#c29f62", "#83ba9e")
cols02 <- c("#c62d17", "#023f75", "#ea894e", "#266b69", "#eb4601", "#f6c619")
cols03 <- c("#fa6e01", "#2f2f2f", "#972b1d", "#e6a84b", "#4c211b", "#ff717f")
cols04 <- c("#223e9c", "#b12b23", "#aebea6", "#edae11", "#0f6657", "#c74732")
cols05 <- c("#6a73cf", "#edd064", "#0eb0c8", "#f2ccac", "#a1d5b9", "#e1abbc")

all_colors <- c(cols01, cols02, cols03, cols04, cols05)
all_colors
colorname2=c("#46A040","#00AF99" ,"#98D9E9" ,"#F6313E", "#FFA300","#2f2f2f",  "#FFC179",
"#FF5A00", "#663366","#FF6666","#8F1336", "#0081C9", "#001588", "#CC0033",
"#CC9966","#CC0033","#999933","#009966","#CCCC33","#CCFF99","#0eb0c8","#993333","#333366",
"#490C65", "#BA7FD0","#A6CEE3", "#1F78B4", "#DE77AE", "#B2DF8A", "#006D2C", "#B5AD64",
"#9DA8E2","#91C392", "#194a55", "#187c65", "#c29f62","#f49128","#c62d17", "#333329",
"#023f75", "#ea894e", "#266b69",  "#e1abbc","#f6c619","#fa6e01", "#972b1d", "#e6a84b", "#4c211b", "#ff717f","#223e9c",
"#aebea6", "#edae11", "#c74732", "#6a73cf", "#edd064" ,
"#0eb0c8", "#f2ccac" ,"#868686","#339966", "#83ba9e", "#b12b23", "#0f6657","#f26115","#eb4601")
colorname3=colorname2
length(colorname3)


colorsseq1 <-c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594")
colorsseq2<-c("#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594")


colbgyor1=c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142")

cellt <- c('#67ba54','#f8ed68','#d596c1','#fbe925','#ecb8d4','#5abb88',
      '#d7d7d6','#66c3af','#f5a5a1','#d9e27b','#bf9ac6','#93c73b',
      '#9075b2','#a0cfed','#f9ab5d','#27b5ad','#f38861','#22c1f2',
      '#fdd271','#72cbdd','#dfc521','#97cd8b')


library(RColorBrewer)
colors_main <- c("#A5D7E5", "#FEC567")
colors_multi <- c("#a5d7e5",
                  "#fec567",
                  "#db7f9b",
                  "#b1babe",
                  "#afdd91",
                  "#c48ec0",
                  "#be8f87",
                  "#b2ecd4",
                  "#ffefba",
                  "#e3c6ea",
                  "#acb791",
                  "#f4b9c1",
                  "#a9c7cd",
                  "#77a1bf",
                  "#e0c79b",
                  "#fb7272",
                  "grey50")
colors_insulin <- c("#badce6", "#54a2b8")
colors_mval <- c("5.71" = "#a6dbbb",
                 "7" = "#3cb67b",
                 "9.29" = "#359566",
                 "12.4" = "#2e7651")
colors_adi <- c("Adipocyte 1" = "#BFABC3",
                "Adipocyte 2" = "#907198",
                "Adipocyte 3" = "#62376E")

color_low <- "#441153"
color_low2 <- "#62376e"
color_low3 <- "#772A82"
color_high <- "#3CB67B"
color_high2 <- "#1D793D"
colors_ramp <- colorRampPalette(colors = c(color_low, color_high))


color_func_blue <- colorRampPalette(colors = c("#C6DBEF", "#075a84"))
color_func_orange <- colorRampPalette(colors = c("#F3E55C", "#E8602D"))
color_func_pink <- colorRampPalette(colors = c("#bfabc3", "#62376e"))
color_func_green <- colorRampPalette(colors = c("#a6dbbb", "#359566"))
color_func_purp <- colorRampPalette(colors = c("#edbbcb", "#c480a3"))
