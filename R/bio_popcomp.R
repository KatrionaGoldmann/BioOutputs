#' Function to compare multiple populations
#'
#' This function generates a volcano plot from a top table using ggplot. 
#' @param data A data frame containing x, y, and population data
#' @param xlab Column name or index of the x variable
#' @param ylab Column name or index of the y variable
#' @param poplab Column name or index of the population variable
#' @keywords volcano
#' @import ggplot2
#' @import ggpubr
#' @import ggrepel
#' @export
#' @examples
#' bio_volcano(toptable)

pop="pop"

# Function to create brackets
bracketsGrob <- function(...){
	l <- list(...)
	e <- new.env()
	e$l <- l
	grid:::recordGrob(  {	do.call(grid.brackets, l)}, e)
}

# Significance values
star <- function(pval) {
	if (pval <= 0.0001) {	return("****")}
	if (pval <= 0.001) {	return("*** ")}
	if (pval <= 0.01) {	return("**  ")}		
	if (pval <= 0.05) {	return("*   ")}
	else {return("    ")
	}
}

#df <- data.frame("x"=rep(c(1:10), 3), "y"=c(1:30)/10, "pop"=rep(c("a", "b", "c"), each=10))

bio_popcomp <- function(x){
	
	lmList(y ~ x|pop, data=df)
	
	
	ggplot(df, aes(x, y, group=pop, color=pop)) +
		geom_point() +
		geom_line()
	
	
	
	
	
	
	temp2 <- t(data.frame(temp[["Pr(>F)"]]))
	if(is.null(title)) {title <- col}
	dimnames(temp2) <- list(c(paste(title, "pvalues")), rownames(temp))
	aov.test <<- rbind(aov.test, temp2)
	
	
	
	
	p.out <- p.out[p.out$group1 == 0 | p.out$group2 == 0, ]
	p.out[p.out$group1 == 0, c("group1", "group2")] <- p.out[p.out$group1 == 0, c("group2", "group1")]
	p.vals <<- rbind(p.vals, cbind(p.out, "Comparison"=rep(col, nrow(p.out)), "Condition"=rep(title, nrow(p.out))))
	p.out.t <- data.frame(p.out[! grepl("ns|NS", p.out$p.signif), c("group1", "group2")])
	p.out.t <- split(as.character(unlist(p.out.t)), seq(nrow(p.out.t)))
	
	treated.text <- c(p.out$p.signif[p.out$group1 == 16 & p.out$Treatment == "Treated"], p.out$p.signif[p.out$group1 == 24 & p.out$Treatment == "Treated"], 
										p.out$p.signif[p.out$group1 == 36 & p.out$Treatment == "Treated"], p.out$p.signif[p.out$group1 == 48 & p.out$Treatment == "Treated"])
	non.treated.text <- c(p.out$p.signif[p.out$group1 == 16 & p.out$Treatment == "Untreated"], p.out$p.signif[p.out$group1 == 24 & p.out$Treatment == "Untreated"], 
												p.out$p.signif[p.out$group1 == 36 & p.out$Treatment == "Untreated"], p.out$p.signif[p.out$group1 == 48 & p.out$Treatment == "Untreated"])
	
	
	max <- 1.01*1.05*max(df.sum$y +df.sum$se, na.rm=T )
	min <- 0.99*min(df.sum$y -df.sum$se, na.rm=T )
	
	max.p <- min(c( ((df.sum$y[df.sum$Time == 48 & df.sum$Treatment=="Treated"]-min)/(max-min)), ((df.sum$y[df.sum$Time == 48 & df.sum$Treatment=="Untreated"]-min)/(max-min)) ))
	min.p <- max(c((df.sum$y[df.sum$Time == 48 & df.sum$Treatment=="Treated"]-min)/(max-min), (df.sum$y[df.sum$Time == 48 & df.sum$Treatment=="Untreated"]-min)/(max-min)))
	
	b1 <- bracketsGrob(0.85, min.p, 0.85, max.p, h=0.05, lwd=1, col="maroon", type=2, curvature=0.5)
	
	hm <- NA
	if(col %in% names(healthy.mean)) hm <- healthy.mean[col]
	
	timepoints <- as.numeric(as.character(unique(df$Time)))
	plot <- ggplot(aes(x=as.numeric(Time), y=y, group=Treatment, colour=Treatment), data=df.sum) +
		geom_errorbar(aes(ymin=y-se, ymax=y+se), width=.7) + 
		ylab(gsub("_", " ", col)) + 
		xlab("Time") + 
		theme(plot.title=element_text(hjust=0.5)) +
		scale_x_continuous(breaks=c(0,16,24,36,48), labels=c(0,16,24,36,48), limits = c(0, 56)) +
		scale_y_continuous(limits = c(c(min(c(hm, df.sum$y-df.sum$se), na.rm=T), 1.05*max(c(hm, df.sum$y+df.sum$se), na.rm=T)))) +
		scale_color_manual(values=c("salmon", "steelblue3")) +
		ggtitle(gsub("_", " ", title)) +
		geom_line(show.legend=NA) +
		geom_point() + 
		annotate("text", x = timepoints[timepoints != 0], y=Inf, label=gsub("ns", "", treated.text), colour="salmon", size=6,  hjust = 0.5, vjust = 2) +
		annotate("text", x = timepoints[timepoints != 0], y=Inf, label=gsub("ns", "", non.treated.text), colour="steelblue3", size=6,  hjust = 0.5, vjust = 3) +
		theme(legend.position="bottom", legend.box = "horizontal") 
	
	
	if(temp2[, "Time:Treatment"] < 0.05){	
		stars <- paste()
		plot <- plot + annotation_custom(b1) +  annotate("text", x = 53, y = ((min)+((min.p + (max.p-min.p)/2)*(max-min))), label=paste("'", star(temp2[, "Time:Treatment"]), "'", sep=""), parse =TRUE, color="maroon", size=6, hjust = 0)  
	}
	mean.df <<- rbind(mean.df, data.frame("Param"=rep(col, nrow(df.sum)), df.sum))
	
	
	if(mean.line == TRUE & col %in% names(healthy.mean)){
		plot <- plot + geom_segment(mapping=aes(x=0, y=hm, xend=48, yend=hm, linetype="Healthy Subjects"), size=0.5, color="navy") + 
			scale_linetype_manual("",values=c("Healthy Subjects"="dashed"))
		
	}
	
	plot
	plots.all[[paste("both", col)]] <<- plot
}
