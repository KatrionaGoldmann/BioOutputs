#' Function to generate a volcano plot
#'
#' This function generates a volcano plot from a top table using ggplot. 
#' @param toptable A data frame containing p value and fold change columns for parameters compared across multiple groups. The p value column should be named "pvalue". 
#' @param fc.col The column name which stores the fold change. Should be in the log2 format (default="log2FC")
#' @param padj.col The column  which contains adjusted p-values. If NULL adjusted pvalues will be calculated
#' @param padj.method correction method. Options include: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). Default="fdr
#' @param padj.cutoff The cutoff for adjusted pvalues. This adds a horizontal line of significance (default=NULL)
#' @param fc.cutoff The log2(fold change) significance cut-off (default=1)
#' @param marker.colour Character vector of four colours to map to the volcano plot. In the order non-significanct, fold-change significant, pvalue significant, significant in fold-change and pvalues (default=c("grey60", "olivedrab", "salmon", "darkturquoise"))
#' @param label.p.cutoff The cutoff for adjusted pvalues for labelling (default=NULL). Not recommended if many significant rows. 
#' @param label.row.indices Indices of rows to be labelled (default=NULL)
#' @param label.colour Colour of labels (default="black")
#' @param legend.labs A character vector for theThe legend label names (default=c("Not Significant", "FC>fc.cutoff", "Padj<padj.cutoff", "FC>fc.cutoff& Padj<padj.cutoff"))
#' @param add.lines Whether to add dashed lines at fc.cutoff and padj.cutoff (default=TRUE)
#' @param line.colour The color of dashed significance lines (default="grey14")
#' @param main Plot title
#' @param xlims,ylims The plot limits
#' @keywords volcano
#' @import ggplot2
#' @import ggrepel
#' @export
#' @examples
#' bio_volcano(toptable)

bio_volcano<- function(toptable, fc.col="log2FC", padj.col=NULL, padj.method="fdr",
											padj.cutoff=0.05, fc.cutoff=1, 
											marker.colour=c("grey60", "olivedrab", "salmon", "darkturquoise"), 
											label.p.cutoff=NULL, label.row.indices=NULL, label.colour="black", legend.labs=NULL,
											add.lines=TRUE, line.colour="grey14", 
											main="Volcano Plot", xlims=NULL, ylims=NULL){
	
	if(! fc.col %in% colnames(toptable)) stop("Could not find fold change column in toptable")
	if(is.null(padj.col)) toptable$padj <- p.adjust(toptable$pvalue, method=padj.method) else{
	if(! padj.col %in% colnames(toptable) ) stop("Could not find the Padj column")}
	
	# Introduce different forms of significance
	toptable$Significance <- "Not Significant"
	toptable$Significance[(abs(toptable[, fc.col]) > fc.cutoff) & !is.na(toptable$pvalue)] <- paste("FC>", fc.cutoff, sep="")
	toptable$Significance[(toptable$padj < padj.cutoff)] <- paste("Padj<", padj.cutoff, sep="")
	toptable$Significance[(toptable$padj<padj.cutoff) & (abs(toptable[, fc.col])>fc.cutoff)] <- paste("FC>", fc.cutoff, "& Padj<", padj.cutoff, sep="")
	
	avail <- unique(toptable$Significance)

	toptable$Significance <- factor(toptable$Significance, levels=c("Not Significant",  paste("FC>", fc.cutoff, sep=""), paste("Padj<", padj.cutoff, sep=""), paste("FC>", fc.cutoff, "& Padj<", padj.cutoff, sep="")))

	marker.colour = setNames(marker.colour, levels(toptable$Significance))
	if(is.null(legend.labs)) legend.labs=levels(droplevels(toptable$Significance))
	
	toptable$logFC <- toptable[, fc.col]
	
	plot <- ggplot(toptable, aes(x=logFC, y=-log10(pvalue), color=Significance) ) +
		
		geom_point(alpha=1/2, size=1.) +	scale_color_manual(labels=legend.labs, values=marker.colour) +
		
		# Set the x and y limits. 
		scale_x_continuous(limits=xlims) + scale_y_continuous(limits=ylims) +
		
		# Modify various aspects of the plot text and legend
		theme_classic() +
		theme(plot.title=element_text(hjust=0.5)) +
		labs(title = main, x=~Log[2]~"(Fold Change)", y="-"~Log[10]~"(P)")
	
	# Add lines 
	if(add.lines == TRUE){
		plot <- plot + geom_hline(yintercept=min(-log10(toptable$pvalue[toptable$padj < padj.cutoff]), na.rm=T), size= 0.3, linetype="longdash", colour="grey14") 
		plot <- plot + geom_vline(xintercept=c(-fc.cutoff, fc.cutoff), size= 0.3, linetype="longdash", colour="grey14") 
	}
	
	# Add labels
	if(! is.null(label.row.indices)	){
		plot <- plot + geom_text_repel(data=toptable[label.row.indices, ], aes(label= rownames(toptable)[label.row.indices]), size=2, vjust=0.5, hjust = -0.2, color=label.colour) 
	}
	
	if(! is.null(label.p.cutoff)){
		plot <- plot + geom_text_repel(data=toptable[toptable$padj < label.p.cutoff, ], aes(label= rownames(toptable)[toptable$padj < label.p.cutoff] ), size=2, vjust=0.5, hjust = -0.2, color=label.colour) 
	}
	
	plot
}
