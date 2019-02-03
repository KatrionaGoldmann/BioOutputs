#' Function to compare two treatment groups
#'
#' This function plots a line graph comparing two treatment groups over time 
#' @param df Data frame containing the data for both groups
#' @param y.col Column name corresponding to the y-axis values in df
#' @param x.col Column name corresponding to the x-axis values in df
#' @param group.col Column name which corresponding to the group of values in df
#' @param id.col Column name which corresponds to subject ID
#' @param main Title of pot
#' @param cols Character vector for colours of lines
#' @param p.col Colour of p-value text
#' @keywords line 
#' @import ggplot2
#' @import Rmisc
#' @import ggpubr
#' @import pBrackets
#' @export

bio_treatmentGroups <- function(df, group.col, y.col="y", x.col="x", id.col="id", main="", cols=NULL, p.col = "maroon" ){
 
  if(is.null(cols)) {cols = c("steelblue", "salmon")}
      
  df <- df[, c(x.col, y.col, group.col, id.col)]
  colnames(df) <- c("x.var", "y.var", "g.var", "id")
  df$g.var <- as.character(df$g.var)
  df$y.var <- as.numeric(as.character(df$y.var))
  df$x.var <- as.numeric(as.character(df$x.var))
  
  df <- df[! is.na(df$y.var), ]
  df <- df[! is.na(df$x.var), ]
  df <- df[! is.na(df$g.var), ]

  t1 <- anova(lmerTest::lmer(y.var ~ x.var + g.var + g.var*x.var + (1 | id), data=df))
  t2 <- t(data.frame(t1[["Pr(>F)"]]))
  if(is.null(main)) {main <- paste(y.col, "among groups")}
  dimnames(t2) <- list(c(paste(main, "pvalues")), rownames(t1))

  df.sum <- summarySE(df, measurevar="y.var", groupvars=c("x.var", "g.var"), na.rm=T)

  max <- 1.01*1.05*max(df.sum$y +df.sum$se, na.rm=T )
  min <- 0.99*min(df.sum$y -df.sum$se, na.rm=T )
  
  max.p <- min((df.sum$y.var[df.sum$x.var == max(df.sum$x.var[df.sum$g.var==unique(df$g.var)[1]], na.rm=T) & df.sum$g.var==unique(df$g.var)[1]]-min)/(max-min))
  min.p <- min((df.sum$y.var[df.sum$x.var == max(df.sum$x.var[df.sum$g.var==unique(df$g.var)[2]], na.rm=T) & df.sum$g.var==unique(df$g.var)[2]]-min)/(max-min))
  
  b1 <- bracketsGrob(0.8, 1.2*max(c(min.p, max.p)), 0.8, 1.2*min(c(min.p, max.p)), h=0.05, lwd=1, col=p.col, type=2, curvature=0.5)
  
  timepoints <- as.numeric(as.character(unique(df$x.var)))
  timepoint <- sort(timepoints)
  plot <- ggplot(aes(x=as.numeric(x.var), y=y.var, group=g.var, colour=g.var), data=df.sum) +
    geom_errorbar(aes(ymin=y.var-se, ymax=y.var+se), width=.7) + 
    ylab(y.col) + 
    xlab(x.col) + 
    theme(plot.title=element_text(hjust=0.5)) +
    scale_x_continuous(limits = c(min(timepoints, na.rm=T), 1.22*max(timepoints, na.rm=T))) +
    scale_color_manual(name = group.col, values=cols) +
    ggtitle(main) +
    geom_line(show.legend=NA) +
    geom_point() + 
    theme(legend.position="bottom", legend.box = "horizontal") 

  plot <- plot + annotation_custom(b1) + annotate("text", x = 1.1*max(timepoints, na.rm=T), y = ((min)+((min.p + (max.p-min.p)/2)*(max-min))), 
                                           label=paste("'p=", format(t2[, 3], digits=3), "'", sep=""), parse =TRUE, color=p.col, size=4, hjust = 0)  
 
  plot
}


bracketsGrob <- function(...){
  l <- list(...)
  e <- new.env()
  e$l <- l
  grid:::recordGrob(  {	do.call(grid.brackets, l)}, e)
}
