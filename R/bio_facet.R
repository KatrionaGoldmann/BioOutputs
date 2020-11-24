#' Function to create gg-style facet objects 
#'
#' @param text The facet title
#' @param colour Facet background colour
#' @param fontcolour Font colour
#' @param angle Text angle/direction (0=horizontal, 90=vertical)
#' @param align Text alignment
#' @param size Font size
#' @keywords plane, scatter
#' @import ggplot2
#' @export
#' @examples
#' blank = bio_facet("", "white") 
#' col1 = bio_facet("Sepal.Length", "grey80", angle=0) + theme(plot.margin=margin(0, 2, 0, 2))
#' col2 = bio_facet("Sepal.Width", "grey80", angle=0) + theme(plot.margin=margin(0, 2, 0, 2))
#' row = bio_facet("Species", "grey80") 
#' length=ggplot(iris, aes(x=Species, y=Sepal.Length)) + geom_boxplot()
#' width=ggplot(iris, aes(x=Species, y=Sepal.Width)) + geom_boxplot()
#' ggarrange(blank, col1, col2, row, length, width, ncol=3, heights=c(0.05, 1), widths=c(0.1, 1, 1))



bio_facet = function(text, colour="grey80", fontcolour="black", angle=90, align="center", size=5){
  p = ggplot(df=NULL)  + 
    theme_void() + 
    annotate("text", x=0, y=1, label=text, angle=angle, colour=fontcolour, size=size) + 
    theme(panel.background=element_rect(fill=colour, color=colour)) 
  if(align == "left") p = p + lims(x=c(0, 1))
  if(align == "right") p = p + lims(x=c(-1, 0))
  return(p)
}


