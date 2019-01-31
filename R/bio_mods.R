#' Function to investigate expression for custom modules
#'
#' This function exports a complex heatmap looking at the expression of different modules which can be customised
#' @param exp data frame containing the expression data
#' @param mod.list A list of modules. Each element contains the list of genes for a modules. The gene names must match the rownames in the exp dataframe.
#' @param meta Dataframe where each column contains an annotation/tract for samples in the heatmap. The order of samples in meta must match that of exp. 
#' @param mean.var Character defining the meta column to average (mean) over. 
#' @param cluster.rows The method to use for clustering of rows
#' @param cols Chacter vector, or named vector to fix the order, defining the colours of each mean.var group. 
#' @param show.names Show row names. Logical.
#' @param main Title of heatmap
#' @keywords heatmap, module
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @export

bio_mods <- function(exp, mod.list, meta, mean.var,cluster.rows=TRUE, cols=c("red", "blue", "green3"),  main="", show.names=TRUE){
  #req(length(mod.list) > 0 )
  
  #try(nrow(meta) == ncol(exp))
  mod.names <- names(mod.list)

  all.genes <- as.character(unlist(mod.list))
  exp.sub <- exp[match(all.genes, rownames(exp)), ]
  
  
  split <- c()
  for (i in mod.names){split <- c(split, rep(i, (length(mod.list[[i]][mod.list[[i]] %in% all.genes])) )) 		}
  split <- split[rowSums(is.na(exp.sub)) != ncol(exp.sub)]
  
  exp.sub <- exp.sub[rowSums(is.na(exp.sub)) != ncol(exp.sub),   ]

  # ungraded has been removed, to reintroduce see version 6.3
  exp.sub <- t(scale(t(exp.sub)))
  split <- split[rowSums(is.na(exp.sub)) != ncol(exp.sub)]
  exp.sub <- exp.sub[rowSums(is.na(exp.sub)) != ncol(exp.sub),   ]
  
  p <- c()
  r.list <- list()
  for(i in unique(meta[, mean.var])){
    p <- c(p, which(meta[, mean.var] == i))
    r.list[[i]] <- exp.sub[, which(meta[, mean.var] == i)]
  }
  
  exp.sub <- exp.sub[,p]
  cutoff <- (quantile(exp.sub, na.rm=TRUE)[5] - quantile(exp.sub, na.rm=TRUE)[4])/2 + quantile(exp.sub, na.rm=TRUE)[4]
  clim<-c(-1*cutoff, 0, cutoff)
  col_fun <- circlize::colorRamp2(clim, c("blue", "white", "red"))

  if(is.null(cols)){
    n <- length(unique(meta[, mean.var]))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    cols=setNames(sample(col_vector, n), unique(meta[, mean.var]))
  } else if(is.null(names(cols)) ) cols = setNames(cols, unique(meta[, mean.var]))

  cols = list(cols)
  names(cols) = mean.var
  ha2 <- ComplexHeatmap::HeatmapAnnotation(Pathotype=meta[p,mean.var], show_legend = FALSE, col=cols)
  
  hm.a.list <- list()
  for(i in unique(meta[, mean.var])){
    ccol = setNames(as.character(sapply(cols, "[[", i)), i)
    ccol = list(ccol)
    names(ccol) = mean.var
    hm.a.list[[i]] <- ComplexHeatmap::HeatmapAnnotation(Pathotype=meta[which(meta[, mean.var] == i), mean.var], 
                                                        show_legend = FALSE, col=ccol)
  }
  

  mean.list = lapply(1:length(unique(meta[, mean.var])),  function(x) x=data.frame())
  names(mean.list) = unique(meta[, mean.var])
  
  for(m in mod.names){
    for(i in unique(meta[, mean.var])){
      mean.list[[i]] <- rbind(  mean.list[[i]], t(data.frame(colMeans(r.list[[i]][split == m, ], na.rm=TRUE))))
    
      rownames(mean.list[[i]])[nrow(mean.list[[i]])] <- paste(m, "mean Z-score")
    }
  }
  
  ht_list = NULL #ComplexHeatmap::Heatmap(rbind(r.list[[1]], mean.list[[1]]))  ## Heatmap(...) + NULL gives you a HeatmapList object
  
  for(s in 1:length(unique(meta[, mean.var]))) {
    show.rows=FALSE
    if (s==length(unique(meta[, mean.var]))) show.rows = show.names
    ht_list = ht_list + ComplexHeatmap::Heatmap(rbind(r.list[[s]], mean.list[[s]]), cluster_rows=cluster.rows, show_column_names = FALSE, top_annotation_height = unit(10, 'mm'), top_annotation = hm.a.list[[s]],
                                                col=col_fun, split = c(split, rep(" ", nrow(mean.list[[s]]))), column_title=paste("(n=", ncol(r.list[[s]]), ")", sep=""),
                                                clustering_distance_rows = "pearson", #column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
                                                clustering_method_rows = "ward.D2", show_heatmap_legend = FALSE,  cluster_columns = TRUE, row_names_gp = grid::gpar(fontsize=6),
                                                row_dend_width = unit(10, "mm"), show_row_dend = FALSE, show_row_names=show.rows)
  }
  
  temp <- ComplexHeatmap::Heatmap(rbind(r.list[[1]], mean.list[[1]]), col=col_fun)

  anno_legend_list = lapply(ha2@anno_list, function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
  heatmap_legend = color_mapping_legend(temp@matrix_color_mapping, title="Z-score", plot = FALSE, color_bar="continuous")
  
  padding = unit.c(unit(3, "mm"), unit(0, "mm"), unit(c(0, 12), "mm"))

  draw(ht_list, heatmap_legend_list = c(list(heatmap_legend), anno_legend_list ),
         heatmap_legend_side = "left", column_title=main,
         column_title_gp = grid::gpar(fontsize = 18), gap=unit(c(1, 1, 1), "mm"), padding=padding)
  
  
  
  
  
}


