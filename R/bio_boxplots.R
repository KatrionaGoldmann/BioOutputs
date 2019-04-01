#' Comparison Boxplots
#' @description Creates boxplots showing the significance between groups. Adapted from kassambara/ggpubr. 
#' @param data 
#' @param x,y x and y variable names for drawing.
#' @param p.cutoff plot p-value if above p.cutoff threshold. To include all comparisons set as NULL. 
#' @param stars Logical. Whether significance shown as numeric or stars
#' @param method a character string indicating which method to be used for comparing means.c("t.test", "wilcox.test")
#' @param star.vals a list of arguments to pass to the function symnum for symbolic number coding of p-values. 
#' For example, the dafault is symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c('****', '***', '**', '*', 'ns')).
#' In other words, we use the following convention for symbols indicating statistical significance:
#' ns: p > 0.05; *: p <= 0.05; **: p <= 0.01; ***: p <= 0.001; ****: p <= 0.0001
#' @import ggpubr
#' @import gtools
#' @examples
#' bio_boxplots(iris, x="Species", y= "Sepal.Width", p.cutoff = 0.0001)
#' bio_boxplots(iris, x="Species", y= "Sepal.Width", NULL, stars=TRUE)

bio_boxplots = function(data, x, y, p.cutoff=0.5, stars=FALSE, method="t.test", star.vals=NULL){
  if(is.null(p.cutoff)) p.cutoff = 1
  
  df = data.frame("x"=data[, x], "y"=data[, y])
  mc = combinations(n=length(unique(df$x)), r=2,v=unique(as.character(df$x)),repeats.allowed=F)
  mc = lapply(1:nrow(mc), function(x) as.character(mc[x,]))
  
  # Plot only the significant comparisons
  keep = c()
  for(i in 1:length(mc)){
    
    # Only include if enough variables
    if( all( table(droplevels(df$x[df$x %in% mc[[i]]])) > 1) ){
      keep = c(keep, compare_means(y ~ x, data = df[df$x %in% mc[[i]], ], method=method)$p <= p.cutoff)
    } else{keep = c(keep, FALSE)}
  }
  mc = mc[keep]
  
  if(stars==TRUE) {
    if(is.null(star.vals)){
      symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
    } else{symnum.args=star.vals}
  } else{symnum.args=list()}
  
  ggboxplot(df, x = "x", y = "y", color = "x", add = "jitter") +
    labs(x=x, y=y) +
    stat_compare_means(comparisons = mc, method="t.test", symnum.args=symnum.args) + 
    theme(legend.title = element_blank(), legend.position="none")
}
