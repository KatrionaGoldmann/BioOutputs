#' Function to capitalize phrases
#'
#' This function puts character vectors into title appropriate capitalization. 
#' @param titles A character vector of phrases to be converted to titles
#' @param exeption.words A character vector of words with case to be forced (for example abbreviations and roman numerals)
#' @param replace.chars A named list of characters to replace in title. This works in order of appearance. E.g. c("\\."=" ") replaces fullstops with spaces. 
#' @keywords titles
#' @export
#' @examples
#' shakespeare_plays =c("all's well that ends well", "as you like it","comedy of errors","love's labour's lost", "measure for measure", "merchant of venice","merry wives of windsor","midsummer night's dream","much ado about nothing","taming of the shrew", "tempest", "twelfth night","two gentlemen of verona", "winter's tale", "henry iv, part i","henry iv, part ii", "henry v", "henry vi, part i","henry vi, part ii", "henry vi, part iii", "henry viii","king john", "pericles","richard ii", "richard iii", "antony and cleopatra","coriolanus","cymbeline", "hamlet","julius caesar", "king lear", "macbeth (the scottish play)", "othello", "romeo and juliet","timon of athens", "titus andronicus", "troilus and cressida")
#' exception.words= c("I", "II", "III", "IV", "V", "VI") # Pass in exceptions
#' 
#' bio_capitalize(shakespear_plays, exception.words)


bio_capitalize = function(titles, exception.words=c(), replace.chars=c("\\."=" ")){
  
  titles = as.character(titles)
  
  # words not to be capitalized according to Chicago Style
  not.caps = c("a", "aboard", "about", "above", "across", "after", "against", "along", "amid", "among", "an", 
               "and", "anti", "around", "as", "at", "before", "behind", "below", "beneath", "beside", "besides", 
               "between", "beyond", "but", "by", "concerning", "considering", "despite", "down", "during", "except", 
               "excepting", "excluding", "following", "for", "from", "in", "inside", "into", "like", "minus", "near", 
               "of", "off", "on", "onto", "opposite", "or", "outside", "over", "past", "per", "plus", "regarding", 
               "round", "save", "since", "so", "than", "the", "through", "to", "toward", "towards", "under", 
               "underneath", "unlike", "until", "up", "upon", "versus", "via", "with", "within", "without", "yet")
  
  cap = function(x){
    # Remove unwanted characters
    for(i in names(replace.chars)){ x = gsub(i, replace.chars[[i]], x)}
    
    # Split into words
    words = unlist(strsplit(x, "[[:space:]]|(?=[,.!?(])", perl=TRUE))
    words = words[words != ""]
    
    # Capitalize words
    title.words = sapply(words, function(y) {
      if(! tolower(y) %in% not.caps) {  
        output = paste0(toupper(substr(y, 1, 1)), tolower(substr(y, 2, nchar(y))))
      } else { output = tolower(y)}
      
      # Now account for the exception words
      if(tolower(output) %in% tolower(exception.words)){
        output = exception.words[which(tolower(exception.words) == tolower(output))]
      }
      output
    })
    
    # Capitalise first letter even if in list of non-caps
    if(! words[1] %in% exception.words){
      title = paste0(toupper(substr(paste(title.words, collapse=" "), 1, 1)), substr(paste(title.words, collapse=" "), 2, nchar(paste(title.words, collapse=" "))))
      } else{ title = paste(title.words, collapse=" ")}
    
    # Correct white spaces surrounding punctuation
    title = gsub(" (?=[,.!?])", "", title, perl=TRUE)
    title = gsub("\\( ", "(", title, perl=TRUE)
  }
  
  
  sapply(titles, cap)
  
}






