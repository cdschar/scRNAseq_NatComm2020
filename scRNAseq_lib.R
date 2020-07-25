#functions for custom plotting of Monocle objects

library(monocle)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(tibble)
library(viridis)



monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

# Plots clusters of cells .
plot_cell_clusters <- function(cds, 
                               x=1, 
                               y=2, 
                               color_by="Cluster", 
                               markers=NULL, 
                               show_cell_names=FALSE, 
                               cell_size=1.5,
                               cell_name_size=2,
                               limits=c(0,1),
                               breaks = NULL,
                               colors = NULL,
                               color_by_scale = FALSE,
                               ...){
  if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 0){
    stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
  }
  
  #set up variables
  gene_short_name <- NULL
  sample_name <- NULL
  data_dim_1 <- NULL
  data_dim_2 <- NULL
  
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info <- pData(cds)
  
  #extract tSNE coordinates for each cell into a data.frame
  tSNE_dim_coords <- reducedDimA(cds)
  data_df <- data.frame(t(tSNE_dim_coords[c(x,y),]))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- colnames(cds)
  data_df <- merge(data_df, lib_info, by.x="sample_name", by.y="row.names")
  
  #if markers variable used...
  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      cds_subset <- cds[row.names(markers_fData),]
      if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        integer_expression <- TRUE
      }
      else {
        integer_expression <- FALSE
        
      }
      if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        
        if (is.null(sizeFactors(cds_subset))) {
          stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
        }
        cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
      }
      else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
      }
      markers_exprs <- cds_exprs
      #markers_exprs <- reshape2::melt(as.matrix(cds_exprs))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  
  #other options for plotting by markers...
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + facet_wrap(~feature_label) 
  }else{
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
  }

#set up color scheme
if(is.null(colors) == TRUE) {
  #cols <- c("#4393C3","#4393C3","#F7FCB9","#D6604D", "#D6604D")
  cols <- c("#4575B4","#4575B4","#F7FCB9","#D73027", "#D73027")   #blue, gray, red
  #cols <- c("#FEE5D9", "#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15")   #all red
  #cols <- c("#969696", "#006d2c") #grey, green
  }else{
    cols = colors #use user supplied color scheme
  }
  
  #scale color_by variable if needed
  if(color_by_scale == "log10") {
  	color_by = log10(color_by + 0.1)
   }else{
    color_by = color_by
  }
 
  # FIXME: setting size here overrides the marker expression funtionality. 
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    g <- g + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE) + 
      scale_colour_gradientn(colors=cols,na.value = "grey80", name = paste0("log10(value + 0.1)"), ...)
  }else if ( length(color_by)>1 ){
    g <- g + geom_point(aes(color = color_by), size=I(cell_size), na.rm = TRUE) +
    scale_colour_gradientn(colors=cols, limits = breaks, oob = scales::squish)
  }else if ( !color_by %in% c("UMI", "num_genes_expressed", "Pseudotime") ){
    g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
  }
  
  g <- g + 
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() + 
    xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(text = element_text(size = 15))
  g
}

#funtion to plot expression of multiple genes along pseudo time
plot_genes_time <- function(cds, genes, min_expr=NULL, label_by_short_name=TRUE,
                            trend_formula="~ sm.ns(Pseudotime, df=3)", panel_order=NULL,
                            relative_expr=TRUE, max_norm=TRUE, lwd=1, export = NULL){
  
  gene_table = fData(cds)
  gids = gene_table$id[gene_table$gene_short_name %in% genes]
  cds_subset = cds[gids, ]
  f_id <- NA
  Cell <- NA
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    integer_expression <- TRUE
  } else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  } else {
    cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  #cds_exprs$f_id <- as.character(cds_exprs$f_id)
  #cds_exprs$Cell <- as.character(cds_exprs$Cell)
  
  if (integer_expression) {
    cds_exprs$adjusted_expression <- cds_exprs$expression
  } else {
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    } else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  } else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
  model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                                       relative_expr = T, new_data = new_data)
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame("expectation"=model_expectation[x$f_id, x$Cell]))
  cds_exprs <- merge(cds_exprs, expectation)
  #cds_exprs$expectation <- expectation#apply(cds_exprs,1, function(x) model_expectation[x$f_id, x$Cell])
  
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  
  #exprs = as.matrix(exprs(cds)[gids, ])
  #pt = pData(cds)$Pseudotime
  #exprs = as.data.frame(t(rbind(exprs, pt)))
  #colnames(exprs)<-c(gene_table[gids,]$gene_short_name, "pt")
  #df1 = melt(exprs, id="pt")
  #colnames(df1)<-c("Pseudo_time", "Genes", "Expression")
  colnames(cds_exprs)[which(colnames(cds_exprs)=="gene_short_name")]="Genes"
  
  if(!is.null(export)){
    gene = unique(cds_exprs$Genes)[1]
    x1 = which(cds_exprs$Genes==gene)
    df = cds_exprs[x1,c("Cell","Size_Factor","num_genes_expressed","UMI","Cluster","Pseudotime","State")]
    for(gene in unique(cds_exprs$Genes)){
      x1 = which(cds_exprs$Genes==gene)
      gene_df = cds_exprs[x1,c("expression","expectation")]
      colnames(gene_df) = c( paste0(gene, "_expression"), paste0(gene, "_expectation") )
      df = cbind(df, gene_df)
    }
    write.table(df, export, row.names=F, sep="\t", quote=F)
  }
  
  if(max_norm){
    for(gene in unique(cds_exprs$Genes)){
      x1 = which(cds_exprs$Genes==gene)
      cds_exprs$expectation[x1] = cds_exprs$expectation[x1]/max(cds_exprs$expectation[x1])
    }
  }
  
  p = ggplot(cds_exprs, aes(Pseudotime, expectation, colour=Genes)) + geom_line(size=lwd) +
      ylab("Relative Expression")+monocle_theme_opts()
  return(p)

}


#funtion to plot average expression of multiple genes along pseudo time
plot_genes_avg_time <- function(cds, genes, min_expr=NULL, label_by_short_name=TRUE,
                            trend_formula="~ sm.ns(Pseudotime, df=3)", panel_order=NULL,
                            relative_expr=TRUE ){
  
  gene_table = fData(cds)
  gids = gene_table$id[gene_table$gene_short_name %in% genes]
  cds_subset = cds[gids, ]
  f_id <- NA
  Cell <- NA
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    integer_expression <- TRUE
  } else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  } else {
    cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  #cds_exprs$f_id <- as.character(cds_exprs$f_id)
  #cds_exprs$Cell <- as.character(cds_exprs$Cell)
  
  if (integer_expression) {
    cds_exprs$adjusted_expression <- cds_exprs$expression
  } else {
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    } else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  } else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
  model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                                       relative_expr = T, new_data = new_data)
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame("expectation"=model_expectation[x$f_id, x$Cell]))
  cds_exprs <- merge(cds_exprs, expectation)
  #cds_exprs$expectation <- expectation#apply(cds_exprs,1, function(x) model_expectation[x$f_id, x$Cell])
  
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  
  colnames(cds_exprs)[which(colnames(cds_exprs)=="gene_short_name")]="Genes"
  avg_exprs=aggregate(expectation~Cell, data=cds_exprs, FUN=mean)
  pt=aggregate(Pseudotime~Cell, data=cds_exprs, FUN=mean)
  df1=merge(avg_exprs, pt)
  
  p = ggplot(df1, aes(Pseudotime, expectation)) + geom_line(size=2, col="blue")+
    ylab("Relative Expression")+monocle_theme_opts()
  return(p)
  
}

data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
}

plot_genes_violin <- function( cds_subset, grouping = "State", min_expr = NULL, cell_size = 0.75, 
                              nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
                              plot_trend = FALSE, label_by_short_name = TRUE, relative_expr = TRUE, 
                              log_scale = TRUE, show.legend=False) 
{
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression = TRUE
  }
  else {
    integer_expression = FALSE
    relative_expr = TRUE
  }
  if (integer_expression) {
    cds_exprs = exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs = Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    #cds_exprs = reshape2::melt(round(as.matrix(cds_exprs)))
    cds_exprs = reshape2::melt(as.matrix(cds_exprs))
  }
  else {
    cds_exprs = exprs(cds_subset)
    cds_exprs = reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)) {
    min_expr = cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) = c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] = min_expr
  cds_pData = pData(cds_subset)
  
  # 
  # # Custom bit for adding in a group for 
  # if(! is.null(show_combined)) {
  #   for(combine_gene in show_combined) {
  #     cds_pData_all = subset(cds_pData, gene == combine_gene)
  #     cds_pData_all[, grouping] = paste("All", combine_gene)
  #     cds_pData = rbind(cds_pData, cds_pData_all)
  #   }
  # }
  
  cds_fData = fData(cds_subset)
  cds_exprs = merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs = merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression = log10(cds_exprs$expression)
  
  
  
  
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label = cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] = cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label = cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label = cds_exprs$f_id
  }
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label = factor(cds_exprs$feature_label, 
                                     levels = panel_order)
  }
  q = ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q = q + geom_violin(aes_string(fill = color_by), scale = "width", show.legend=show.legend)
  }
  else {
    q = q + geom_violin(scale = "width", show.legend=show.legend)
  }
  if (plot_trend == TRUE) {
    q = q + stat_summary(fun.data = "mean_cl_boot", 
                         size = 0.2)
    q = q + stat_summary(aes_string(x = grouping, y = "expression", 
                                    group = color_by), fun.data = "mean_cl_boot", 
                         size = 0.2, geom = "line")
  }
  q = q + facet_wrap(~feature_label, nrow = nrow, 
                     ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
     q = q + expand_limits(y = c(min_expr, 1))
  }
  
  
  q = q + ylab("Expression") + xlab(grouping) + stat_summary(fun.data=data_summary) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  if (log_scale == TRUE){
    
    q = q + scale_y_log10()
  }
  q
}

plot_gray <- function (cds, x = 1, y = 2, color_by = "Cluster", markers = NULL, 
    show_cell_names = FALSE, cell_size = 1.5, cell_name_size = 2, 
    ...) 
{
    if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 
        0) {
        stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
    }
    gene_short_name <- NULL
    sample_name <- NULL
    data_dim_1 <- NULL
    data_dim_2 <- NULL
    lib_info <- pData(cds)
    tSNE_dim_coords <- reducedDimA(cds)
    data_df <- data.frame(t(tSNE_dim_coords[c(x, y), ]))
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- colnames(cds)
    data_df <- merge(data_df, lib_info, by.x = "sample_name", 
        by.y = "row.names")
    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in% 
            markers)
        if (nrow(markers_fData) >= 1) {
            cds_subset <- cds[row.names(markers_fData), ]
            if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                "negbinomial.size")) {
                integer_expression <- TRUE
            }
            else {
                integer_expression <- FALSE
            }
            if (integer_expression) {
                cds_exprs <- exprs(cds_subset)
                if (is.null(sizeFactors(cds_subset))) {
                  stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
                }
                cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
                cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
            }
            else {
                cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
            }
            markers_exprs <- cds_exprs
            colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
            markers_exprs <- merge(markers_exprs, markers_fData, 
                by.x = "feature_id", by.y = "row.names")
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
            by.y = "cell_id")
    	data_df$value <- with(data_df, ifelse(value >= 0.01, value, NA))
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
            facet_wrap(~feature_label)
    }
    else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        g <- g + geom_point(aes(color = log10(value + 0.1)), 
            size = I(cell_size), na.rm = TRUE) + scale_color_viridis(option = "viridis", name = "log10(values + 0.1)", na.value = "grey80", end = 0.8)
    }
    else {
        g <- g + geom_point(aes_string(color = color_by), size = I(cell_size), 
            na.rm = TRUE)
    }
    g <- g + monocle_theme_opts() + xlab(paste("Component", x)) + 
        ylab(paste("Component", y)) + theme(legend.position = "top", 
        legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
        theme(panel.background = element_rect(fill = "white")) + 
        theme(text = element_text(size = 15))
    g
}

plot_cell_trajectory <- function(cds, 
                                 x=1, 
                                 y=2, 
                                 color_by="State", 
                                 show_tree=TRUE, 
                                 show_backbone=TRUE, 
                                 backbone_color="black", 
                                 markers=NULL, 
                                 use_color_gradient = FALSE,
                                 markers_linear = FALSE,
                                 show_cell_names=FALSE,
                                 show_state_number = FALSE,
                                 cell_size=1.5,
                                 cell_link_size=0.75,
                                 cell_name_size=2,
                                 state_number_size = 2.9,
                                 show_branch_points=TRUE,
                                 theta = 0,
                                 magic_array = c(NA),
                                 ...) {
  requireNamespace("igraph")
  gene_short_name <- NA
  sample_name <- NA
  sample_state <- pData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info_with_pseudo <- pData(cds)
  
  if (is.null(cds@dim_reduce_type)){
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  
  if (cds@dim_reduce_type == "ICA"){
    reduced_dim_coords <- reducedDimS(cds)
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree") ){
    reduced_dim_coords <- reducedDimK(cds)
  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  
  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  
  dp_mst <- minSpanningTree(cds)
  
  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }
  
  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    select_(source = "from", target = "to") %>%
    left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")
  
  data_df <- t(monocle::reducedDimS(cds)) %>%
    as.data.frame() %>%
    select_(data_dim_1 = x, data_dim_2 = y) %>%
    rownames_to_column("sample_name") %>%
    mutate(sample_state) %>%
    left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")
  
  return_rotation_mat <- function(theta) {
    theta <- theta / 180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }
  rot_mat <- return_rotation_mat(theta)
  
  cn1 <- c("data_dim_1", "data_dim_2")
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
  
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData),])))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    if(length(magic_array)>1){
      data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
      data_df$value = magic_array
    } else{
      data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    }
    if(use_color_gradient) {
      #cols <- c("#4393C3","#4393C3","#F7FCB9","#D6604D", "#D6604D")
      cols <- c("#4575B4","#4575B4","#F7FCB9","#D73027", "#D73027")   #blue, gray, red
      
      if(markers_linear){
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color= value), size=I(cell_size), na.rm = TRUE) + 
          scale_colour_gradientn(colors=cols) + facet_wrap(~feature_label)
        } else {
          g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE) + 
              scale_colour_gradientn(colors=cols) + facet_wrap(~feature_label)
        }
    } else {
      if(markers_linear){
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size= (value * 0.1))) + facet_wrap(~feature_label)
      } else {
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2, size=log10(value + 0.1))) + facet_wrap(~feature_label)
      }
    }
  } else {
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
  }
  if (show_tree){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=cell_link_size, linetype="solid", na.rm=TRUE, data=edge_df)
  }
  
  # FIXME: setting size here overrides the marker expression funtionality. 
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    if(use_color_gradient) {
      # g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    } else {
      g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    }
  }else {
    if(use_color_gradient) {
      # g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
    } else {
      g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
    }
  }
  
  
  if (show_branch_points && cds@dim_reduce_type == 'DDRTree'){
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- ica_space_df %>%
      slice(match(mst_branch_nodes, sample_name)) %>%
      mutate(branch_point_idx = seq_len(n()))
    
    g <- g +
      geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                 size=5, na.rm=TRUE, branch_point_df) +
      geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="branch_point_idx"),
                size=4, color="white", na.rm=TRUE, branch_point_df)
  }
  if (show_cell_names){
    g <- g + geom_text(aes(label=sample_name), size=cell_name_size)
  }
  if (show_state_number){
    g <- g + geom_text(aes(label = sample_state), size = state_number_size)
  }
  
  g <- g + 
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() + 
    xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) +
    theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}
