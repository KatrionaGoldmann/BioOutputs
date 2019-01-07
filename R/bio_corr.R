
NULL
#' Correlation plot
#' @description Create a correlation plot. Taken from kassambara/ggpubr just changed the default arguments
#' @param x,y x and y variables for drawing.
#' @param color,fill point colors.
#' @param shape point shape. See \code{\link{show_point_shapes}}.
#' @param point logical value. If TRUE, show points.
#' @param rug logical value. If TRUE, add marginal rug.
#' @param add allowed values are one of "none", "reg.line" (for adding linear
#'   regression line) or "loess" (for adding local regression fitting).
#' @param add.params parameters (color, size, linetype) for the argument 'add';
#'   e.g.: add.params = list(color = "red").
#' @param conf.int logical value. If TRUE, adds confidence interval.
#' @param conf.int.level Level controlling confidence region. Default is 95\%.
#'   Used only when add != "none" and conf.int = TRUE.
#' @param fullrange should the fit span the full range of the plot, or just the
#'   data. Used only when add != "none".
#' @param ellipse logical value. If TRUE, draws ellipses around points.
#' @param ellipse.level the size of the concentration ellipse in normal
#'   probability.
#' @param ellipse.type Character specifying frame type. Possible values are
#'  \code{"convex"}, \code{"confidence"} or types supported by
#'  \code{\link[ggplot2]{stat_ellipse}()} including one of \code{c("t", "norm",
#'  "euclid")} for plotting concentration ellipses.
#'
#'  \itemize{ \item \code{"convex"}: plot convex hull of a set o points. \item
#'  \code{"confidence"}: plot confidence ellipses arround group mean points as
#'  \code{\link[FactoMineR]{coord.ellipse}()}[in FactoMineR]. \item \code{"t"}:
#'  assumes a multivariate t-distribution. \item \code{"norm"}: assumes a
#'  multivariate normal distribution. \item \code{"euclid"}: draws a circle with
#'  the radius equal to level, representing the euclidean distance from the
#'  center. This ellipse probably won't appear circular unless
#'  \code{\link[ggplot2]{coord_fixed}()} is applied.}
#' @param ellipse.alpha Alpha for ellipse specifying the transparency level of
#'   fill color. Use alpha = 0 for no fill color.
#' @param ellipse.border.remove logical value. If TRUE, remove ellipse border lines.
#' @param mean.point logical value. If TRUE, group mean points are added to the
#'   plot.
#' @param mean.point.size numeric value specifying the size of mean points.
#' @param star.plot logical value. If TRUE, a star plot is generated.
#' @param star.plot.lty,star.plot.lwd line type and line width (size) for star
#'   plot, respectively.
#' @param label the name of the column containing point labels. Can be also a
#'   character vector with length = nrow(data).
#' @param font.label a vector of length 3 indicating respectively the size
#'   (e.g.: 14), the style (e.g.: "plain", "bold", "italic", "bold.italic") and
#'   the color (e.g.: "red") of point labels. For example \emph{font.label =
#'   c(14, "bold", "red")}. To specify only the size and the style, use
#'   font.label = c(14, "plain").
#' @param font.family character vector specifying font family.
#' @param label.select character vector specifying some labels to show.
#' @param repel a logical value, whether to use ggrepel to avoid overplotting
#'   text labels or not.
#' @param label.rectangle logical value. If TRUE, add rectangle underneath the
#'   text, making it easier to read.
#' @param cor.coef logical value. If TRUE, correlation coefficient with the
#'   p-value will be added to the plot.
#' @param cor.coeff.args a list of arguments to pass to the function
#'   \code{\link{stat_cor}} for customizing the displayed correlation
#'   coefficients. For example: \code{cor.coeff.args = list(method = "pearson",
#'   label.x.npc = "right", label.y.npc = "top")}.
#' @param cor.method method for computing correlation coefficient. Allowed
#'   values are one of "pearson", "kendall", or "spearman".
#' @param cor.coef.coord numeric vector, of length 2, specifying the x and y
#'   coordinates of the correlation coefficient. Default values are NULL.
#' @param cor.coef.size correlation coefficient text font size.
#' @param ggp a ggplot. If not NULL, points are added to an existing plot.
#' @param show.legend.text logical. Should text be included in the legends? NA,
#'   the default, includes if any aesthetics are mapped. FALSE never includes,
#'   and TRUE always includes.
#' @param ... other arguments to be passed to \code{\link[ggplot2]{geom_point}}
#'   and \code{\link{ggpar}}.
#' @import dplyr
#' @details The plot can be easily customized using the function ggpar(). Read
#'   ?ggpar for changing: \itemize{ \item main title and axis labels: main,
#'   xlab, ylab \item axis limits: xlim, ylim (e.g.: ylim = c(0, 30)) \item axis
#'   scales: xscale, yscale (e.g.: yscale = "log2") \item color palettes:
#'   palette = "Dark2" or palette = c("gray", "blue", "red") \item legend title,
#'   labels and position: legend = "right" \item plot orientation : orientation
#'   = c("vertical", "horizontal", "reverse") }
#' @seealso \code{\link{stat_cor}}, \code{\link{stat_stars}}, \code{\link{stat_conf_ellipse}} and \code{\link{ggpar}}.
#' @examples
#' # Load data
#' data("mtcars")
#' df <- mtcars
#' df$cyl <- as.factor(df$cyl)
#' head(df[, c("wt", "mpg", "cyl")], 3)
#'
#' # Basic plot
#' # +++++++++++++++++++++++++++
#' ggscatter(df, x = "wt", y = "mpg",
#'    color = "black", shape = 21, size = 3, # Points color, shape and size
#'    add = "reg.line",  # Add regressin line
#'    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#'    conf.int = TRUE, # Add confidence interval
#'    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
#'    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
#'    )
#'
#' # loess method: local regression fitting
#' ggscatter(df, x = "wt", y = "mpg",
#'    add = "loess", conf.int = TRUE)
#'
#'
#' # Control point size by continuous variable values ("qsec")
#' ggscatter(df, x = "wt", y = "mpg",
#'    color = "#00AFBB", size = "qsec")
#'
#'
#' # Change colors
#' # +++++++++++++++++++++++++++
#' # Use custom color palette
#' # Add marginal rug
#' ggscatter(df, x = "wt", y = "mpg", color = "cyl",
#'    palette = c("#00AFBB", "#E7B800", "#FC4E07") )
#'
#'
#'
#'
#' # Add group ellipses and mean points
#' # Add stars
#' # +++++++++++++++++++
#' ggscatter(df, x = "wt", y = "mpg",
#'    color = "cyl", shape = "cyl",
#'    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#'    ellipse = TRUE, mean.point = TRUE,
#'    star.plot = TRUE)
#'
#'
#' # Textual annotation
#' # +++++++++++++++++
#' df$name <- rownames(df)
#' ggscatter(df, x = "wt", y = "mpg",
#'    color = "cyl", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#'    label = "name", repel = TRUE)
#'
#'
#' @export
# bio_corr <- function(data, x, y, combine = FALSE, merge = FALSE,
# 											color = "black", fill = "lightgray", palette = NULL,
# 											shape = 19, size = 2, point = TRUE,  rug = FALSE,
# 											title = NULL, xlab = NULL, ylab = NULL,
# 											facet.by = NULL, panel.labs = NULL, short.panel.labs = TRUE,
# 											add = c("none", "reg.line", "loess"), add.params = list(),
# 											conf.int = FALSE, conf.int.level = 0.95, fullrange = TRUE,
# 											ellipse = FALSE, ellipse.level = 0.95,
# 											ellipse.type = "norm", ellipse.alpha = 0.1,
# 											ellipse.border.remove = FALSE,
# 											mean.point = FALSE, mean.point.size = ifelse(is.numeric(size), 2*size, size),
# 											star.plot = FALSE, star.plot.lty = 1, star.plot.lwd = NULL,
# 											label = NULL,  font.label = c(12, "plain"), font.family = "",
# 											label.select = NULL, repel = FALSE, label.rectangle = FALSE,
# 											cor.coef = FALSE, cor.coeff.args = list(), cor.method = "pearson", cor.coef.coord = c(NULL, NULL), cor.coef.size = 4,
# 											ggp = NULL, show.legend.text = NA,
# 											ggtheme = theme_pubr(),
# 											...){
# 	
# 	
# 	# Default options
# 	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 	.opts <- list(
# 		combine = combine, merge = merge,
# 		color = color, fill = fill, palette = palette,
# 		title = title, xlab = xlab, ylab = ylab,
# 		facet.by = facet.by, panel.labs = panel.labs, short.panel.labs = short.panel.labs,
# 		shape = shape, size = size, point = point,  rug = rug,
# 		add = add, add.params = add.params,
# 		conf.int = conf.int, conf.int.level = conf.int.level, fullrange = fullrange,
# 		ellipse = ellipse, ellipse.level = ellipse.level,
# 		ellipse.type = ellipse.type, ellipse.alpha = ellipse.alpha,
# 		ellipse.border.remove = ellipse.border.remove,
# 		mean.point = mean.point, mean.point.size = mean.point.size,
# 		star.plot = star.plot, star.plot.lty = star.plot.lty, star.plot.lwd = star.plot.lwd,
# 		label = label, font.label = font.label, font.family = font.family,
# 		label.select = label.select, repel = repel, label.rectangle = label.rectangle,
# 		cor.coef = cor.coef, cor.coeff.args = cor.coeff.args, cor.method = cor.method,
# 		cor.coef.coord = cor.coef.coord, cor.coef.size = cor.coef.size,
# 		ggp = ggp, show.legend.text = show.legend.text, ggtheme = ggtheme, ...)
# 	
# 	if(!missing(data)) .opts$data <- data
# 	if(!missing(x)) .opts$x <- x
# 	if(!missing(y)) .opts$y <- y
# 	
# 	# User options
# 	#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 	.user.opts <- as.list(match.call(expand.dots = TRUE))
# 	.user.opts[[1]] <- NULL # Remove the function name
# 	# keep only user arguments
# 	for(opt.name in names(.opts)){
# 		if(is.null(.user.opts[[opt.name]]))
# 			.opts[[opt.name]] <- NULL
# 	}
# 	
# 	font.label <- .parse_font(font.label) %>% .compact()
# 	font.label$color <- ifelse(is.null(font.label$color), color, font.label$color)
# 	.opts$font.label <- font.label
# 	
# 	.opts$fun <- ggscatter_core
# 	if(missing(ggtheme) & (!is.null(facet.by) | combine))
# 		.opts$ggtheme <- theme_pubr(border = TRUE)
# 	p <- do.call(.plotter, .opts)
# 	if(.is_list(p) & length(p) == 1) p <- p[[1]]
# 	return(p)
# 	
# }


bio_corr <- function(data, x, y,
													 color = "black", fill = "lightgray", palette = NULL,
													 shape = 19, size = 2, point = TRUE,  rug = FALSE,
													 title = NULL, xlab = NULL, ylab = NULL,
													 add = "reg.line", add.params = list(color="red", fill = "lightgray"),
													 conf.int = TRUE, conf.int.level = 0.95, fullrange = TRUE,
													 ellipse = FALSE, ellipse.level = 0.95,
													 ellipse.type = "norm", ellipse.alpha = 0.1,
													 ellipse.border.remove = FALSE,
													 mean.point = FALSE, mean.point.size = ifelse(is.numeric(size), 2*size, size),
													 star.plot = FALSE, star.plot.lty = 1, star.plot.lwd = NULL,
													 label = NULL,  font.label = c(12, "plain"), font.family = "",
													 label.select = NULL, repel = FALSE, label.rectangle = FALSE,
													 cor.coef = TRUE, cor.coeff.args = list(), cor.method = "pearson", cor.coef.coord = c(NULL, NULL), cor.coef.size = 4,
													 ggp = NULL, show.legend.text = NA,
													 ggtheme = theme_classic(),
													 ...)
{
	
	library(RCurl)
	
	#Load in the scripts from ggpubr
	for(i in c("stat_cor", "utilities", "utilities_label", "utilities_color", "utilities_base", "set_palette")){
		script <- getURL(paste("https://raw.githubusercontent.com/kassambara/ggpubr/master/R/", i, ".R", sep=""), ssl.verifypeer = FALSE)
		eval(parse(text = script))
	}
	
	add <- match.arg(add)
	add.params <- .check_add.params(add, add.params, error.plot = "", data, color, fill, ...)
	
	if(length(label) >1){
		if(length(label) != nrow(data))
			stop("The argument label should be a column name or a vector of length = nrow(data). ",
					 "It seems that length(label) != nrow(data)")
		else data$label.xx <- label
		label <- "label.xx"
	}
	# label font
	font.label <- .parse_font(font.label)
	font.label$size <- ifelse(is.null(font.label$size), 12, font.label$size)
	font.label$color <- ifelse(is.null(font.label$color), color, font.label$color)
	font.label$face <- ifelse(is.null(font.label$face), "plain", font.label$face)
	
	if(is.null(ggp)) p <- ggplot(data, aes_string(x, y))
	else p <- ggp
	
	if(point) p <- p +
		.geom_exec(geom_point, data = data, x = x, y = y,
							 color = color, fill = fill, size = size,
							 shape = shape, ...)
	
	# Adjust shape when ngroups > 6, to avoid ggplot warnings
	if(shape %in% colnames(data)){
		ngroups <- length(levels(data[, shape]))
		if(ngroups > 6) p <- p + scale_shape_manual(values=1:ngroups, labels = levels(data[, shape]))
	}
	
	# Add marginal rug
	# +++++++++++
	if(rug) p <- p + .geom_exec(geom_rug, data = data,
															color = color, size = size/2)
	
	# Add reg line or loess
	# ++++++++++++
	if(add %in% c("reg.line", "loess")){
		add <- ifelse(add == "reg.line", stats::lm, stats::loess)
		if(is.null(add.params$linetype)) add.params$linetype <- "solid"
		
		.args <- .geom_exec(NULL, data = data,
												se = conf.int, level = conf.int.level,
												color = add.params$color, fill = add.params$fill,
												linetype = add.params$linetype, size = add.params$size,
												fullrange = fullrange)
		
		mapping <- .args$mapping
		option <- .args$option
		option[["method"]] <- add
		option[["mapping"]] <- do.call(ggplot2::aes_string, mapping)
		p <- p + do.call(geom_smooth, option)
	}
	
	
	# Add ellipses
	# +++++++++++
	if(ellipse){
		grp <- intersect(unique(c(color, fill, shape)), colnames(data))[1]
		# NO grouping variable
		if(is.na(grp)) {
			grp <- factor(rep(1, nrow(data)))
			grp_name <- "group"
			data$group <- grp
		}
		# Case of grouping variable
		else {
			grp_name <- grp
			data[, grp_name] <- as.factor(data[, grp_name])
		}
		
		if (ellipse.type == 'convex')
			p <- p + .convex_ellipse(data, x, y, grp_name, color, fill, ellipse.alpha,
															 ellipse.border.remove = ellipse.border.remove)
		else if(ellipse.type == "confidence")
			p <- p + .confidence_ellipse(data, x, y, grp_name, color, fill,
																	 alpha = ellipse.alpha, level = ellipse.level,
																	 ellipse.border.remove = ellipse.border.remove)
		else if (ellipse.type %in% c('t', 'norm', 'euclid'))
			p <- p + .stat_ellipse(data, x, y, grp_name, color = color, fill = fill,
														 alpha = ellipse.alpha, type = ellipse.type, level = ellipse.level,
														 ellipse.border.remove = ellipse.border.remove)
	}
	# /ellipse
	
	# Add mean points
	# +++++++++
	if(mean.point) {
		p <- p + .geom_exec(stat_mean, data = data,
												color = color, shape = shape, fill = fill,
												size = mean.point.size)
	}
	
	# Star plots
	# ++++++++++++
	if(star.plot){
		p <- p + .geom_exec(stat_stars, data = data,
												color = color, linetype = star.plot.lty, size = star.plot.lwd)
	}
	
	#/ star plots
	
	# Add textual annotation
	# ++++++
	alpha <- 1
	if(!is.null(list(...)$alpha)) alpha <- list(...)$alpha
	if(!is.null(label)) {
		lab_data <- data
		# Select some labels to show
		if(!is.null(label.select))
			lab_data  <- subset(lab_data, lab_data[, label, drop = TRUE] %in% label.select,
													drop = FALSE)
		
		if(repel){
			ggfunc <- ggrepel::geom_text_repel
			if(label.rectangle) ggfunc <- ggrepel::geom_label_repel
			p <- p + .geom_exec(ggfunc, data = lab_data, x = x, y = y,
													label = label, fontface = font.label$face,
													size = font.label$size/3, color = font.label$color,
													alpha = alpha, family = font.family,
													box.padding = unit(0.35, "lines"),
													point.padding = unit(0.3, "lines"),
													force = 1, show.legend = show.legend.text, seed=123)
		}
		else{
			ggfunc <- geom_text
			vjust  <- -0.7
			if(label.rectangle) {
				ggfunc <- geom_label
				vjust <- -0.4
			}
			p <- p + .geom_exec(ggfunc, data = lab_data, x = x, y = y, color = color,
													label = label, fontface = font.label$face, family = font.family,
													size = font.label$size/3, color = font.label$color,
													vjust = vjust, alpha = alpha, show.legend = show.legend.text)
		}
	}
	
	# Add correlation coefficient
	if(cor.coef){
		
		if(!missing(cor.method))
			cor.coeff.args$method <- cor.method
		if(!missing(cor.coef.size))
			cor.coeff.args$size <- cor.coef.size
		if(!missing(cor.coef.coord)){
			cor.coeff.args$label.x <- cor.coef.coord[1]
			cor.coeff.args$label.y <- cor.coef.coord[2]
		}
		p <- p + do.call(stat_cor, cor.coeff.args)
	}
	
	p <- ggpar(p, palette = palette, ggtheme = ggtheme,
						 title = title, xlab = xlab, ylab = ylab,...)
	if(font.family != "")
		p <- p + theme(text = element_text(family = font.family))
	p
}

.brewerpal <- function(){
	c(
		# sequential
		'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges',
		'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
		'YlGn', 'YlGnBu YlOrBr', 'YlOrRd',
		#Divergent
		'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
		# Qualitative
		'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'
	)
}
.ggscipal <- function(){
	# Scientific Journal and Sci-Fi Themed Color Palettes for ggplot2
	# ggsci package: https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
	c("npg", "aaas", "nejm", "lancet", "jama", "jco", "ucscgb", "d3", "locuszoom",
		"igv", "uchicago", "startrek", "tron", "futurama", "rickandmorty", "simpsons")
}
.is_color <- function(x) {
	sapply(x, function(X) {
		tryCatch(is.matrix(grDevices::col2rgb(X)),
						 error = function(e) FALSE)
	})
}
.ggfill <- function(palette = NULL, ...) {
	fill_palette(palette = palette, ...)
}

.is_empty <- function(x){
	length(x) == 0
}

.set_axis_limits <- function(xlim = NULL, ylim = NULL){
	if(!is.null(xlim) | !is.null(ylim)) coord_cartesian(xlim, ylim)
}


.add_item <- function(.list, ...){
	pms <- list(...)
	for(pms.names in names(pms)){
		.list[[pms.names]] <- pms[[pms.names]]
	}
	.list
}

.set_legend <- function(p, legend = NULL,
												legend.title = NULL, font.legend = NULL)
{
	if(is.null(legend.title)) legend.title = waiver()
	font <- .parse_font(font.legend)
	
	if(!is.null(legend)) p <- p + theme(legend.position = legend)
	
	if(!.is_empty(legend.title)){
		
		if(.is_list(legend.title)) p <- p + do.call(ggplot2::labs, legend.title)
		else p <- p +
				labs(color = legend.title, fill = legend.title, linetype = legend.title, shape = legend.title)
	}
	
	if(!is.null(font)){
		p <- p + theme(
			legend.text = element_text(size = font$size,
																 face = font$face, colour = font$color),
			legend.title = element_text(size = font$size,
																	face = font$face, colour = font$color)
		)
	}
	
	p
}

.get_gg_xy_variables <- function(p){
	. <- NULL
	x <- p$mapping['x'] %>% as.character() %>% gsub("~", "", .)
	y <- p$mapping['y'] %>% as.character() %>% gsub("~", "", .)
	xy <- c(x, y)
	names(xy) <- c("x", "y")
	return(xy)
}

.set_ticksby <- function(p, xticks.by = NULL, yticks.by = NULL)
{
	.data <- p$data
	# .mapping <- as.character(p$mapping)
	.mapping <- .get_gg_xy_variables(p)
	
	if(!is.null(yticks.by)) {
		y <- .data[, .mapping["y"]]
		ybreaks <- seq(0, max(y, na.rm = TRUE), by = yticks.by)
		p <- p + scale_y_continuous(breaks = ybreaks)
	}
	else if(!is.null(xticks.by)) {
		x <- .data[, .mapping["x"]]
		xbreaks <- seq(0, max(x, na.rm = TRUE), by = xticks.by)
		p <- p + scale_x_continuous(breaks = xbreaks)
	}
	p
}
.set_ticks <-
	function(ticks = TRUE, tickslab = TRUE, font.tickslab = NULL,
					 xtickslab.rt = NULL, ytickslab.rt = NULL,
					 font.xtickslab = font.tickslab, font.ytickslab = font.tickslab)
	{
		
		. <- xhjust <- NULL
		if(!is.null(xtickslab.rt)) {
			if(xtickslab.rt > 5) xhjust <- 1
		}
		else xhjust <- NULL
		
		if (ticks)
			ticks <-
				element_line(colour = "black")
		else
			ticks <- element_blank()
		
		if (is.null(font.xtickslab)) font.x <- list()
		else font.x <- .parse_font(font.xtickslab)
		if (is.null(font.ytickslab)) font.y <- list()
		else font.y <- .parse_font(font.ytickslab)
		
		if (tickslab) {
			xtickslab <- font.x %>% .add_item(hjust = xhjust, angle = xtickslab.rt) %>%
				do.call(element_text, .)
			ytickslab <- font.y %>% .add_item(angle = ytickslab.rt) %>% do.call(element_text, .)
		}
		else {
			xtickslab <- element_blank()
			ytickslab <- element_blank()
		}
		theme(
			axis.ticks = ticks, axis.text.x = xtickslab, axis.text.y = ytickslab
		)
	}

# Add convex ellipse
# data a data frame
# x,y: x and y variables
# grp_name: grp variable
.convex_ellipse <- function(data, x, y, grp_name, color = "black", fill = "lightgray", alpha = 0.1,
														ellipse.border.remove = FALSE ){
	
	grp_levels <- levels(data[, grp_name])
	if(length(grp_levels) == 1) .geom_exec(geomfunc = stat_chull, data = data,
																				 color = color, fill = fill, alpha = alpha,
																				 geom = "polygon")
	else {
		if( ellipse.border.remove) color <- NULL
		else color = grp_name
		.geom_exec(geomfunc = stat_chull, data = data,
							 color = color, fill = grp_name, alpha = alpha,
							 geom = "polygon")
	}
}

# Confidence ellipse
.confidence_ellipse <- function(data, x, y, grp_name, color = "black", fill = "lightgray",
																alpha = 0.1, level = 0.95, ellipse.border.remove = FALSE){
	grp_levels <- levels(data[, grp_name])
	if(length(grp_levels) == 1) {
		mapping <- aes_string(x = x, y = y)
		stat_conf_ellipse(mapping = mapping, data = data,
											color = color, fill = fill, alpha = alpha,
											level = level, geom = "polygon")
	}
	else {
		mapping = aes_string(x = x, y = y, colour = grp_name, fill = grp_name)
		if(ellipse.border.remove ) mapping = aes_string(x = x, y = y,  fill = grp_name)
		stat_conf_ellipse(mapping = mapping, data = data,
											level = level, alpha = alpha,
											geom = 'polygon')
	}
}

# Add ggplot2 stat ellipse
.stat_ellipse <- function(data, x, y, grp_name, color = "black", fill = "lightgray",
													alpha = 0.1, type = "norm", level = 0.95, ellipse.border.remove = FALSE)
{
	grp_levels <- levels(data[, grp_name])
	if(length(grp_levels) == 1){
		mapping <- aes_string(x = x, y = y)
		ggplot2::stat_ellipse(mapping = mapping, data = data,
													level = level, type = type,
													colour = color, fill = fill, alpha = alpha,
													geom = 'polygon')
	}
	else{
		mapping = aes_string(x = x, y = y, colour = grp_name, group = grp_name, fill = grp_name)
		if(ellipse.border.remove) mapping = aes_string(x = x, y = y, colour = NULL, group = grp_name, fill = grp_name)
		ggplot2::stat_ellipse(mapping = mapping, data = data,
													level = level, type = type, alpha = alpha,
													geom = 'polygon')
	}
}


geom_exec <- function (geomfunc = NULL, data = NULL, position = NULL, ...) {
	params <- list(...)
	
	mapping <-
		list() # option to pass to mapping aes() or aes_string()
	option <- list() # option to the geom_*()
	
	allowed_options <- c(
		# general
		"x", "y", "color", "colour", "linetype", "fill", "size", "shape", "width",
		"alpha", "na.rm", "lwd", "pch", "cex", "position", "stat", "geom",
		"show.legend", "inherit.aes", "fun.args", "fontface",
		# boxplot
		"outlier.colour", "outlier.shape", "outlier.size",
		"outlier.stroke", "notch", "notchwidth", "varwidth",
		# dot plot
		"binwidth", "binaxis", "method", "binpositions",
		"stackdir", "stackratio", "dotsize",
		# Violin
		"trim", "draw_quantiles", "scale",
		# error
		"ymin", "ymax", "xmin", "xmax",
		# text
		"label", "hjust", "vjust", "fontface", "angle", "family",
		# text.repel
		"segment.size", "force",
		# smooth
		"se", "level", "fullrange",
		"conf.int.level",
		# straightline
		"xintercept", "yintercept",
		# histograms
		"bins",
		# rug
		"sides",
		# segment
		"arrow", "xend", "yend",
		# stat_summary,
		"fun.data", "fun.y", "fun.ymin", "fun.ymax"
		
	)
	
	columns <- colnames(data)
	for (key in names(params)) {
		value <- params[[key]]
		if (is.null(value)) {
			
		}
		else if (unlist(value)[1] %in% columns & key %in% allowed_options) {
			mapping[[key]] <- value
			
		}
		else if (key %in% allowed_options) {
			option[[key]] <- value
		}
		else if (key =="group")  mapping[[key]] <- value # for line plot
		# else warnings("Don't know '", key, "'")
	}
	if (!is.null(position))
		option[["position"]] <- position
	option[["data"]] <- data
	if(is.null(geomfunc)){
		res <- list(option = option, mapping = mapping)
	}
	else{
		option[["mapping"]] <- do.call(ggplot2::aes_string, mapping)
		res <- do.call(geomfunc, option)
	}
	res
}

.check_add.params <- function(add, add.params, error.plot, data, color, fill,  ...){
	if(color %in% names(data) & is.null(add.params$color))  add.params$color <- color
	if(fill %in% names(data) & is.null(add.params$fill))  add.params$fill <- fill
	if(is.null(add.params$color)) add.params$color <- color
	if(is.null(add.params$fill) & ("crossbar" %in% error.plot | "boxplot" %in% add | "violin" %in% add)) add.params$fill <- fill
	if(is.null(add.params$fill)) add.params$fill <- add.params$color
	#else add.params$fill <- add.params$color
	if(!is.null(list(...)$shape) & is.null(add.params$shape)) add.params$shape <- list(...)$shape
	add.params
}

.parse_font <- function(font){
	if(is.null(font)) res <- NULL
	else if(inherits(font, "list")) res <- font
	else{
		# matching size and face
		size <- grep("^[0-9]+$", font, perl = TRUE)
		face <- grep("plain|bold|italic|bold.italic", font, perl = TRUE)
		if(length(size) == 0) size <- NULL else size <- as.numeric(font[size])
		if(length(face) == 0) face <- NULL else face <- font[face]
		color <- setdiff(font, c(size, face))
		if(length(color) == 0) color <- NULL
		res <- list(size=size, face = face, color = color)
	}
	res
}

ggpar <- function(p, palette = NULL, gradient.cols = NULL,
									main = NULL, submain = NULL, caption = NULL, xlab = NULL, ylab = NULL,
									title = NULL, subtitle = NULL,
									font.main = NULL, font.submain = NULL, font.x = NULL, font.y = NULL, font.caption = NULL,
									font.title = NULL, font.subtitle = NULL, font.family = "",
									xlim = NULL, ylim = NULL,
									xscale = c("none", "log2", "log10", "sqrt"),
									yscale = c("none", "log2", "log10", "sqrt"),
									format.scale = FALSE,
									legend = NULL,
									legend.title = NULL, font.legend = NULL,
									ticks = TRUE, tickslab = TRUE, font.tickslab = NULL,
									font.xtickslab = font.tickslab, font.ytickslab = font.tickslab,
									x.text.angle = NULL, y.text.angle = NULL,
									xtickslab.rt = x.text.angle, ytickslab.rt = y.text.angle,
									xticks.by = NULL, yticks.by = NULL,
									rotate = FALSE,
									orientation = c("vertical", "horizontal", "reverse"),
									ggtheme = NULL,
									...)
{
	
	original.p <- p
	if(rotate) orientation <- "horizontal"
	
	if(is.ggplot(original.p)) list.plots <- list(original.p)
	else if(is.list(original.p)) list.plots <- original.p
	else stop("Can't handle an object of class ", class (original.p))
	if(!is.null(title)) main <- title
	if(!is.null(subtitle)) submain <- subtitle
	if(!is.null(font.title)) font.main <- font.title
	if(!is.null(font.subtitle)) font.submain <- font.subtitle
	if(is.numeric(palette)) palette <- grDevices::palette()[palette]
	
	
	for(i in 1:length(list.plots)){
		p <- list.plots[[i]]
		if(is.ggplot(p)){
			p <- p + .ggcolor(palette)+
				.ggfill(palette)
			if(!is.null(ggtheme)) p <- p + ggtheme # labs_pubr() +
			if(!is.null(gradient.cols)) p <- p + .gradient_col(gradient.cols)
			
			p <- p +.set_ticks(ticks, tickslab, font.tickslab,
												 xtickslab.rt, ytickslab.rt,
												 font.xtickslab = font.xtickslab, font.ytickslab = font.ytickslab)
			p <- .set_ticksby(p, xticks.by, yticks.by)
			p <- p + .set_axis_limits(xlim, ylim)
			p <-.set_legend(p, legend, legend.title, font.legend)
			p <- .set_scale(p, xscale = xscale, yscale = yscale, format.scale = format.scale)
			p <- .labs(p, main, xlab, ylab,
								 font.main, font.x, font.y,
								 submain = submain, caption = caption, font.submain = font.submain, font.caption = font.caption)
			p <- .set_orientation(p, orientation)
			if(font.family != "")
				p <- p + theme(text = element_text(family = font.family))
			list.plots[[i]] <- p
		}
		
	}
	
	if(is.ggplot(original.p)) list.plots[[1]]
	else list.plots
}

.set_orientation <-
	function(p, orientation = c("vertical", "horizontal", "reverse")) {
		ori <- match.arg(orientation)
		if (ori == "horizontal") p + coord_flip()
		else if (ori == "reverse")
			p + scale_y_reverse()
		else p
	}

.labs <- function(p, main = NULL, xlab = NULL, ylab = NULL,
									font.main = NULL, font.x = NULL, font.y = NULL,
									submain = NULL, caption = NULL,
									font.submain = NULL, font.caption = NULL)
{
	
	font.main <- .parse_font(font.main)
	font.x <- .parse_font(font.x)
	font.y <- .parse_font(font.y)
	font.submain <- .parse_font(font.submain)
	font.caption <- .parse_font(font.caption)
	
	if(is.logical(main)){
		if(!main) main <- NULL
	}
	
	if(is.logical(submain)){
		if(!submain) submain <- NULL
	}
	
	if(is.logical(caption)){
		if(!caption) caption <- NULL
	}
	
	
	if (!is.null(main)) {
		p <- p + labs(title = main)
	}
	
	if (!is.null(submain)) {
		p <- p + labs(subtitle = submain)
	}
	
	if (!is.null(caption)) {
		p <- p + labs(caption = caption)
	}
	
	if (!is.null(xlab)) {
		if (xlab == FALSE)
			p <- p + theme(axis.title.x = element_blank())
		else
			p <- p + labs(x = xlab)
	}
	
	if (!is.null(ylab)) {
		if (ylab == FALSE)
			p <- p + theme(axis.title.y = element_blank())
		else
			p <- p + labs(y = ylab)
	}
	
	if (!is.null(font.main))
		p <-
		p + theme(
			plot.title = element_text(
				size = font.main$size,
				lineheight = 1.0, face = font.main$face, colour = font.main$color
			)
		)
	if (!is.null(font.submain))
		p <-
		p + theme(
			plot.subtitle = element_text(
				size = font.submain$size,
				lineheight = 1.0, face = font.submain$face, colour = font.submain$color
			)
		)
	if (!is.null(font.caption))
		p <-
		p + theme(
			plot.caption = element_text(
				size = font.caption$size,
				lineheight = 1.0, face = font.caption$face, colour = font.caption$color
			)
		)
	if (!is.null(font.x))
		p <-
		p + theme(axis.title.x = element_text(
			size = font.x$size,
			face = font.x$face, colour = font.x$color
		))
	if (!is.null(font.y))
		p <-
		p + theme(axis.title.y = element_text(
			size = font.y$size,
			face = font.y$face, colour = font.y$color
		))
	
	
	
	p
}

.set_scale <- function (p, xscale = c("none", "log2", "log10", "sqrt"),
												yscale = c("none", "log2", "log10", "sqrt"),
												format.scale = FALSE)
{
	
	xscale <- match.arg(xscale)
	yscale <- match.arg(yscale)
	.x <- ".x"
	
	if(format.scale){
		if(!requireNamespace("scales")) stop("The R package 'scales' is required.")
		
		if(yscale == "log2"){
			p <- p + scale_y_continuous(trans = scales::log2_trans(),
																	breaks = scales::trans_breaks("log2", function(x) 2^x),
																	labels = scales::trans_format("log2", scales::math_format(2^.x)))
		}
		else if(yscale == "log10"){
			p <- p + scale_y_continuous(trans = scales::log10_trans(),
																	breaks = scales::trans_breaks("log10", function(x) 10^x),
																	labels = scales::trans_format("log10", scales::math_format(10^.x)))
		}
		
		if(xscale == "log2"){
			p <- p + scale_x_continuous(trans = scales::log2_trans(),
																	breaks = scales::trans_breaks("log2", function(x) 2^x),
																	labels = scales::trans_format("log2", scales::math_format(2^.x)))
		}
		else if(xscale == "log10"){
			p <- p + scale_x_continuous(trans = scales::log10_trans(),
																	breaks = scales::trans_breaks("log10", function(x) 10^x),
																	labels = scales::trans_format("log10", scales::math_format(10^.x)))
		}
		
	}
	
	else{
		if(xscale != "none")  p <- p + scale_x_continuous(trans = xscale)
		if(yscale != "none") p <- p + scale_y_continuous(trans = yscale)
	}
	p
}

color_palette <- function(palette = NULL, ...) {
	brewerpal <- .brewerpal()
	ggscipal <- .ggscipal()
	
	res <- NULL
	if (is.null(palette))
		palette <- ""
	if (length(palette) == 1) {
		if (palette %in% brewerpal)
			ggplot2::scale_color_brewer(..., palette = palette)
		else if (palette %in% ggscipal)
			.scale_color_ggsci(palette = palette)
		else if (palette == "grey")
			ggplot2::scale_color_grey(..., start = 0.8, end = 0.2)
		else if (palette == "hue")
			ggplot2::scale_color_hue(...)
		else if(.is_color(palette))
			ggplot2::scale_color_manual(..., values = palette)
	}
	else if (palette[1] != "")
		ggplot2::scale_color_manual(..., values = palette)
}


fill_palette <- function(palette = NULL, ...){
	
	brewerpal <- .brewerpal()
	ggscipal <- .ggscipal()
	
	res <- NULL
	if (is.null(palette))
		palette <- ""
	if (length(palette) == 1) {
		if (palette %in% brewerpal)
			ggplot2::scale_fill_brewer(..., palette = palette)
		else if (palette %in% ggscipal)
			.scale_fill_ggsci(palette = palette)
		else if (palette == "grey")
			ggplot2::scale_fill_grey(..., start = 0.8, end = 0.2)
		else if (palette == "hue")
			ggplot2::scale_fill_hue(...)
		else if(.is_color(palette))
			ggplot2::scale_fill_manual(..., values = palette)
	}
	else if (palette[1] != "")
		ggplot2::scale_fill_manual(..., values = palette)
}


.ggcolor <- function(palette = NULL, ...) {
	color_palette(palette = palette, ...)
}