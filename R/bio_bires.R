#' Function to plot if results above or below cut-off line
#'
#' This function catagorises data points into above or below a cut-off plane
#' @param y x column name in data
#' @param x y column name in data
#' @param df.data Data frame containing x and y columns
#' @param df.plane Date frame modelling the plane
#' @param x.plane column name for the x axis in df.plane
#' @param y.plane column name for the y axis in df.plane
#' @param stepwise logical whether to plot the cutoff plane as stepwise or smoothed
#' @param colours colour vector for higher, lower and plane values (default=c("green", "red", "grey) respectively)
#' @param inc.equal logical whethere points on the line should be counted as above (dafault=TRUE) 
#' @param labels label for the markers (default=c("above", "below"))
#' @param type type of plot for data (options include point (dafault), line, stepwise)
#' @keywords plane, scatter
#' @import ggplot2
#' @export
#' @examples
#' data(beavers)
#' df.plane = beaver1
#' df.data = beaver2
#' df.plane$temp <- df.plane$temp +0.5
#' bio_bires(x="time", y="temp", df.data, x.plane="time", y.plane="temp", df.plane)


bio_bires <- function(x, y, df.data, x.plane, y.plane, df.plane, stepwise=TRUE, colours=c("red", "green", "grey"), inc.equal=TRUE, 
                       labels=c("above", "below"), type = "point"){
	colnames(df.data)[match(c(x, y), colnames(df.data))] = c("x", "y")
	colnames(df.plane)[match(c(x.plane, y.plane), colnames(df.plane))] = c("x", "y")
	
	df.data$col = NA
	# for each x find the nearest x in the plane
	for (i in 1:nrow(df.data)){
	  xval = df.data$x[i]
	  yval = df.data$y[i]
	  yval.plane = df.plane$y[which(abs(df.plane$x[df.plane$x <= xval]-xval)==min(abs(df.plane$x[df.plane$x <= xval]-xval), na.rm=T)  )]
	  if (yval > yval.plane) {df.data$col[i] = labels[1]} else {df.data$col[i] = labels[2]}
	  if (inc.equal == TRUE & yval == yval.plane) df.data$col[i] = labels[1]
	}
  
	df.data = df.data[order(df.data$x), ]
  p = ggplot(df.data, aes(x=x, y=y, col=col, group=1)) 
  if(type == "point"){
    p = p+geom_point() 
  } else if (type == "line") {p = p  + geom_line() } else {p = p  + geom_step()}
  
  p = p + scale_colour_manual(name = "Cutoff" ,values = colours[1:2]) +
  xlab(x) +
  ylab(y) +

  if(stepwise == TRUE){
    geom_step(aes(x, y, col="plane"), data = df.plane, col=colours[3])
  } else{ geom_line(aes(x=x, y=y), data = df.plane, inherit.aes=FALSE, color=colours[3]) }
  p
}
