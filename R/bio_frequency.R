#' Function to generate a frequency table
#'
#' This function generates a frequency table from factor or character vector columns in a data frame
#' @param data A data frame containing columns to be counted
#' @param columns Column names or indices to be counted in data
#' @param freq.percent Whether the table should include frequency counts, percentages or both (options = c("freq", "percent", "both")). Default="both"
#' @param include.na Include NA values (options are TRUE/FALSE, default=TRUE)
#' @param remove.vars Character vector of variables not to be included in the counts (e.g. remove.vars = c("") remove blanks from the count)
#' @keywords frequency, counts, table
#' @export
#' @examples
#' bio_frequency(mtcars, columns=c("vs", "am"))



bio_frequency <- function(data, columns, freq.percent="both", include.na=TRUE, remove.vars=NULL){
	
	# Check the freq.perc option is set up correctly
	if(! freq.percent %in% c("both", "freq", "percent")) stop('freq.percent must be in c("both", "freq", "percent")')
	
	#Check if column levels are identical
	if(length(columns) > 1) {
		check <- sapply(data[, columns], function(x) {
			t <- sort(as.character(unique(x)))
			t <- t[! t %in% remove.vars]
			if(include.na == FALSE) t <- t[! is.na(t)]
			as.character(t)
		})
		if(class(check) == "matrix") {check <- split(check, rep(1:ncol(check), each = nrow(check)))}
		if(all(sapply(2:length(check),function(x) isTRUE(all.equal(check[[x-1]],check[[x]])))) == FALSE) stop("levels of 'columns' must be identical (once remove.vars omitted)")
	}
	
	# Convert to a character vector for handling
	data[, columns]=sapply(data[, columns], FUN = function(x) as.character(x))
	tot.row <- nrow(data)
	
	# Convert empty or NA values to missing
	data[, columns][is.na(data[, columns])]="Missing.na"

	all.vars <- as.character(unique(unlist(data[, columns])))
	all.vars <- sort(all.vars[! all.vars %in% remove.vars])
	if(include.na == FALSE) all.vars <- all.vars[all.vars != "Missing.na"]
	df <- data.frame(matrix(NA, nrow=1, ncol=length(all.vars)+1) )
	dimnames(df)=list("", c(all.vars, "Total"))
	for(i in columns){
		data.col <- data[, i]
		
		# Remove unwanted variables
		for(rem in remove.vars) data.col <- data.col[data.col != "rem"]
		if(include.na == FALSE) data.col <- data.col[data.col != "Missing.na"]
		
		breakdown <- table(data.col)
		if(any(! all.vars %in% names(breakdown))){
			t <- all.vars[which(! all.vars %in% names(breakdown))]
			temp <- data.frame(matrix(0, nrow=1, ncol=length(t), dimnames=list(c(i), c(t))))
			breakdown[t] <- temp
		}
		breakdown <- breakdown[all.vars]
		tot <- sum(unlist(breakdown), na.rm=T)
		if(freq.percent=="both") {out <- t(data.frame(paste(as.character(breakdown), " (", round(100*unlist(breakdown)/tot), "%)", sep="")))
		} else if(freq.percent == "freq") {out <- t(data.frame(paste(as.character(breakdown), sep="")))
		} else if(freq.percent == "percent") {out <- t(data.frame(round(100*as.numeric(as.character(unlist(breakdown)))/tot)))}
		
		dimnames(out) <- list(c(i), c(names(breakdown)))
		out <- cbind(out, "Total" = paste("n =", tot))
		df <- rbind(df, out[, match(colnames(df), colnames(out))])
	}
	df <-df[2:nrow(df), ]
	rownames(df) <- columns
	
	# Rearrange the columns
	if("Missing.na" %in% colnames(df)) df <- df[, c( colnames(df)[! colnames(df) %in% c("Missing.na", "Total")], "Missing.na", "Total")]
	colnames(df)[colnames(df) == "Missing.na"] <- "NA"
	df
}







