.scf_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.scf_source_file <- if (length(.scf_source_files)) tail(.scf_source_files, 1)[[1]] else NULL
.scf_cli_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.scf_start_dir <- if (!is.null(.scf_source_file)) dirname(normalizePath(.scf_source_file)) else if (length(.scf_cli_file)) dirname(normalizePath(sub("^--file=", "", .scf_cli_file[[1]]))) else getwd()
source(file.path(.scf_start_dir, "../../..", "00_Settings/00_setup_R.R"))
scf_check_R_packages(c("Biobase", "ComplexHeatmap", "dplyr", "e1071", "Mfuzz", "purrr", "reshape2", "WGCNA"))

globalVariables(c('.', 'cluster', 'cluster2', 'cluster_name','modulecol'))
clusterData <- function(exp = NULL,
                        scaleData = TRUE,
                        cluster.method = c("mfuzz","kmeans","wgcna"),
                        object = NULL,
                        min.std = 0,
                        cluster.num = NULL,
                        subcluster = NULL,
                        seed = 5201314){
  ComplexHeatmap::ht_opt(message = FALSE)


  cluster.method <- match.arg(cluster.method)

  if(cluster.method == "mfuzz"){


    myset <- Biobase::ExpressionSet(assayData = as.matrix(exp))
    myset <- Mfuzz::filter.std(myset,min.std = min.std,visu = FALSE)


    if(scaleData == TRUE){
      myset <- Mfuzz::standardise(myset)
    }else{
      myset <- myset
    }

    cluster_number <- cluster.num
    m <- Mfuzz::mestimate(myset)


    set.seed(seed)


    mfuzz1 <- function(eset,centers,m,...){
      cl <- e1071::cmeans(Biobase::exprs(eset),centers = centers,method = "cmeans",m = m,...)
    }

    mfuzz_res <- mfuzz1(myset, c = cluster_number, m = m)


    mtx <- Biobase::assayData(myset)
    mtx <- mtx$exprs
    raw_cluster_anno <- cbind(mtx,cluster = mfuzz_res$cluster)


    mem <- cbind(mfuzz_res$membership,cluster2 = mfuzz_res$cluster) %>%
      data.frame(check.names = FALSE) %>%
      dplyr::mutate(gene = rownames(.))


    lapply(1:cluster.num, function(x){
      ms <- mem %>% dplyr::filter(cluster2 == x)
      res <- data.frame(membership = ms[[x]],gene = ms$gene,cluster2 = ms$cluster2,
                        check.names = FALSE)
    }) %>% do.call('rbind',.) -> membership_info


    dnorm <- cbind(myset@assayData$exprs,cluster = mfuzz_res$cluster) %>%
      data.frame(check.names = FALSE) %>%
      dplyr::mutate(gene = rownames(.))


    final_res <- merge(dnorm,membership_info,by = 'gene') %>%
      dplyr::select(-cluster2) %>%
      dplyr::arrange(cluster)


    if(!is.null(subcluster)){
      final_res <- final_res %>% dplyr::filter(cluster %in% subcluster)
    }


    df <- reshape2::melt(final_res,
                         id.vars = c('cluster','gene','membership'),
                         variable.name = 'cell_type',
                         value.name = 'norm_value')


    df$cluster_name <- paste('cluster ',df$cluster,sep = '')


    cltn <- table(final_res$cluster)
    cl.info <- data.frame(table(final_res$cluster))
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn],")",sep = ''))
    }) -> df


    df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",1:nrow(cl.info),
                                                             " (",cl.info$Freq,")",sep = ''))


    return(list(wide.res = final_res,
                long.res = df,
                type = cluster.method,
                geneMode = "none",
                geneType = "none"))
  }else if(cluster.method == "kmeans"){


    exp <- filter.std(exp,min.std = min.std,visu = FALSE)


    if(scaleData == TRUE){
      hclust_matrix <- exp %>% t() %>% scale() %>% t()
    }else{
      hclust_matrix <- exp
    }


    set.seed(seed)
    ht <- ComplexHeatmap::Heatmap(hclust_matrix,
                                  show_row_names = F,
                                  show_row_dend = F,
                                  show_column_names = F,
                                  row_km = cluster.num)


    ht = ComplexHeatmap::draw(ht)
    row.order = ComplexHeatmap::row_order(ht)


    purrr::map_df(1:length(names(row.order)),function(x){
      data.frame(od = row.order[[x]],
                 id = as.numeric(names(row.order)[x]),
                 check.names = FALSE)
    }) -> od.res

    cl.info <- data.frame(table(od.res$id),check.names = FALSE)


    m <- hclust_matrix[od.res$od,]


    wide.r <- m %>%
      data.frame(check.names = FALSE) %>%
      dplyr::mutate(gene = rownames(.),
                    cluster = od.res$id) %>%
      dplyr::arrange(cluster)


    if(!is.null(subcluster)){
      wide.r <- wide.r %>% dplyr::filter(cluster %in% subcluster)
    }


    df <- reshape2::melt(wide.r,
                         id.vars = c('cluster','gene'),
                         variable.name = 'cell_type',
                         value.name = 'norm_value')


    df$cluster_name <- paste('cluster ',df$cluster,sep = '')


    cltn <- table(wide.r$cluster)
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn],")",sep = ''))
    }) -> df


    df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",1:nrow(cl.info),
                                                             " (",cl.info$Freq,")",sep = ''))

    return(list(wide.res = wide.r,
                long.res = df,
                type = cluster.method,
                geneMode = "none",
                geneType = "none"))
  }else if(cluster.method == "wgcna"){

    net <- object
    cinfo <- data.frame(cluster = net$colors + 1,
                        modulecol = WGCNA::labels2colors(net$colors),
                        check.names = FALSE)

    expm <- data.frame(t(scale(exp)))
    expm$gene <- rownames(expm)
    final.res <- data.frame(cbind(expm,cinfo)) %>%
      dplyr::arrange(cluster)


    if(!is.null(subcluster)){
      final.res <- final.res %>% dplyr::filter(cluster %in% subcluster)
    }


    cl.info <- data.frame(table(final.res$cluster))


    df <- reshape2::melt(final.res,
                         id.vars = c('cluster','gene','modulecol'),
                         variable.name = 'cell_type',
                         value.name = 'norm_value')


    df$cluster_name <- paste('cluster ',df$cluster,sep = '')


    cltn <- table(final.res$cluster)
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn]," ",unique(tmp$modulecol),")",sep = ''))
    }) -> df


    df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",1:nrow(cl.info),
                                                             " (",cl.info$Freq,")",sep = ''))


    return(list(wide.res = final.res,
                long.res = df,
                type = cluster.method,
                geneMode = "none",
                geneType = "none"))
  }else{
    message("supply with mfuzz, kmeans or wgcna !")
  }
}


"exps"
