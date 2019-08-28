#' Function to investigate expression for custom modules
#'
#' This function exports a complex heatmap looking at the expression of different modules which can be customised
#' @param exp data frame containing the expression data
#' @param mod.list A list of modules. Each element contains the list of genes for a modules. The gene names must match the rownames in the exp dataframe.
#' @param meta Dataframe where each column contains an annotation/tract for samples in the heatmap. The order of samples in meta must match that of exp. 
#' @param split.var Character defining the meta column to average (mean) over. 
#' @param mean.subjects Logical to determine whether to add a row for the mean value for all subjects in a group 
#' @param cluster.rows The method to use for clustering of rows
#' @param cols Chacter vector, or named vector to fix the order, defining the colours of each mean.var group. 
#' @param show.names Show row names. Logical.
#' @param main Title of heatmap
#' @keywords heatmap, module
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import grid
#' @export

bio_mods <- function(exp, mod.list, meta, mean.var, cluster.rows=TRUE, cols=NULL,  main="", show.names=TRUE, mean.subjects=TRUE){
	
	mod.names <- names(mod.list)
	
	all.genes <- as.character(unlist(mod.list))
	exp.sub <- exp[match(all.genes, rownames(exp)), ]
	
	
	split <- c()
	for (i in mod.names){split <- c(split, rep(i, (length(mod.list[[i]][mod.list[[i]] %in% all.genes])) )) 		}
	split <- split[rowSums(is.na(exp.sub)) != ncol(exp.sub)]
	
	exp.sub <- exp.sub[rowSums(is.na(exp.sub)) != ncol(exp.sub),   ]
	
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
	names(mean.list) = c(as.character(unique(meta[, mean.var])))#, paste(unique(meta[, mean.var]), "mean"))
	
	for(m in mod.names){
		for(i in unique(meta[, mean.var])){
			temp <- rbind(mean.list[[i]], t(data.frame(colMeans(r.list[[i]][split == m, ], na.rm=TRUE))))
			rownames(temp)[nrow(temp)] <- paste(m, "mean Z-score")
			
			if(mean.subjects == TRUE) {
				temp = rbind(temp, rep(mean(as.numeric(temp[paste(m, "mean Z-score"),]), na.rm=T), ncol(temp)))
				rownames(temp)[nrow(temp)] <- paste(m, "mean over subjects")
			}
			mean.list[[i]] <- temp 
			mean.list[[i]] <- mean.list[[i]][c(which(grepl("over subjects", rownames(mean.list[[i]]))),  
																				 which(! grepl("over subjects", rownames(mean.list[[i]])))), ]
		}
	}
	
	
	ht_list = NULL 
	
	for(s in 1:length(unique(meta[, mean.var]))) {
		show.rows=FALSE
		if (s==length(unique(meta[, mean.var]))) show.rows = show.names
		df = rbind(r.list[[s]], mean.list[[s]])
		df = as.matrix(df, dimnames=dimnames(df))
		rownames(df) = gsub("\\.", " ", make.names(gsub("\\..*|\n", "", rownames(df)), unique=F))
		ht_list = ht_list + ComplexHeatmap::Heatmap(df, cluster_rows=cluster.rows, show_column_names = FALSE, top_annotation_height = unit(10, 'mm'),
																								top_annotation = hm.a.list[[s]],
																								col=col_fun, split = c(split, rep(" ", nrow(mean.list[[s]]))), 
																								column_title=paste("(n=", ncol(r.list[[s]]), ")", sep=""),
																								clustering_distance_rows = "pearson", 
																								clustering_method_rows = "ward.D2", show_heatmap_legend = FALSE,  cluster_columns = TRUE, 
																								row_names_gp = grid::gpar(fontsize=6),
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

