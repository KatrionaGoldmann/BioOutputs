library("gtools")
library(ComplexHeatmap)

#' Mean Expression across Groups
#'
#' This function calculates the mean expression (z-score) of sample types
#' @param exp Expression Data
#' @param var Vector classing samples by variables 
#' @return data frame with mean expression levels 
#' @keywords heatmap, fold change, expression, pvalues
bio_mean.by.var = function(exp, var){
	apply(exp, 2, function(x) tapply(x, var, mean))
}

#' Fold Change between Groups
#'
#' This function calculates fold change in expression between sample groups
#' @param exp Expression Data
#' @param var Vector classing samples by variables 
#' @return data frame with mean expression levels 
#' @keywords heatmap, fold change, expression, pvalues
bio_fold.change = function(exp, var){
	mbv = bio_mean.by.var(exp, var)
	comps = combinations(n=length(unique(var)), r=2, v=as.character(unique(var)))
	comps = lapply(1:nrow(comps), function(x) unlist(comps[x, ]))
	names(comps) = unlist(lapply(comps, function(x) paste(x[1], "vs", x[2])))
	
	fc.calc = function(x) {
		vals1 = mbv[x[1], ]
		vals2 = mbv[x[2], ]
		fc = vals1/vals2
	}
	fc = lapply(comps, fc.calc)
	fc.df <- data.frame(matrix(unlist(fc), nrow=length(fc[[1]]), byrow=T))
	colnames(fc.df) = names(fc)
	rownames(fc.df) = make.names(names(fc[[1]]), unique=TRUE)
	fc.df 
}

#' Pvalues between Groups
#'
#' This function calculates pvalues in expression between sample groups
#' @param exp Expression Data
#' @param var Vector classing samples by variables 
#' @return data frame with mean expression levels 
#' @keywords heatmap, fold change, expression, pvalues
bio_p.vals = function(exp, var){
	comps = combinations(n=length(unique(var)), r=2, v=as.character(unique(var)))
	comps = lapply(1:nrow(comps), function(x) unlist(comps[x, ]))
	names(comps) = unlist(lapply(comps, function(x) paste(x[1], "vs", x[2])))
	
	p.calc = function(x) {
		pv = list()
		for(mod in rownames(exp)){
			vals1 = exp[mod, var %in% c(x[1], x[2])]
			var.sub = droplevels(var[var %in% c(x[1], x[2])])
			
			pv[[mod]] =t.test(vals1 ~ var.sub)$p.value
		}
		unlist(pv)
	}
	pv = lapply(comps, p.calc)
	pv.df <- data.frame(matrix(unlist(pv), nrow=length(pv[[1]]), byrow=T))
	dimnames(pv.df) = list(names(pv[[1]]), names(pv))
	pv.df
}

#' Fold change and pvalues between groups
#'
#' This function exports a complex heatmap looking at the expression of different modules which can be customised
#' @param exp Expression Data
#' @param var Vector classing samples by variables 
#' @param stars whether pvales should be written as numeric or start (default=FALSE)
#' @param overlay pvalues on fold change Heatmap or besidde (default = TRUE)
#' @param logp Whether or not to log the pvalues
#' @param prefix Prefix to Heatmap titles
#' @param ... Other parameters to pass to Complex Heatmap
#' @return A list containing the mean expression for each group, the fold change between groups, the pvalues comparing the expression from different groups, a heatmap with mean expression, a heatmap containing the fold change and p-values. 
#' @keywords heatmap, fold change, expression, pvalues
#' @import ComplexHeatmap
#' @import gtools
#' @export
#' @examples
#' m.var = metadata.final[metadata.final$Pathotype %in% c("Fibroid", "Lymphoid"), ]
#' exp = t(rldMatrix[, match(m.var$SampleID..QMUL.or.Genentech., colnames(rldMatrix))])
#' var = droplevels(m.var$Pathotype)
#' 
#' fc = bio_fold.change(exp, var)
#' p = bio_p.vals(exp, var)

bio_fc_heatmap = function(exp, var, prefix="", stars=FALSE, overlay = TRUE, logp = FALSE) 
{
	
  fc = bio_fold.change(t(exp), var)
	p = bio_p.vals(exp, var)
	rownames(p) = rownames(fc) = rownames(exp)
	
	if (stars == TRUE) {
	  pstar = p
		pstar[pstar <= 0.01] = "***"
		pstar[pstar < 0.05 & p > 0.01] = "**"
		pstar[pstar <= 0.05] = "*"
		pstar[pstar > 0.05] = ""
	} else{
		if(logp == TRUE) p = sapply(p, function(x) -log(x), digits=2)
		p = sapply(p, function(x) format(x, digits=2))
		
	}
	
	if(overlay == TRUE){
	hm = Heatmap(fc, column_title=prefix, 
					cell_fun = function(j, i, x, y, w, h, col) {
						grid.text(format(p[i, j], digits=2), x, y)
					}) 
	} else{
		hm = 	Heatmap(fc, column_title=paste(prefix, "Fold Change")) + 
		  Heatmap(p, column_title=paste(prefix, "Pvalues"), col = colorRamp2(c(0, max(p, na.rm=T)), c("white", "seagreen4")), name="-log(p)") 
	}
	hm
}


