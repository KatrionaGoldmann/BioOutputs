% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bio_capitalize.R
\name{bio_capitalize}
\alias{bio_capitalize}
\title{Function to capitalize phrases}
\usage{
bio_capitalize(
  titles,
  exception.words = c(),
  replace.chars = c(`\\\\.` = " "),
  assume.abbrev = TRUE
)
}
\arguments{
\item{titles}{A character vector of phrases to be converted to titles}

\item{replace.chars}{A named list of characters to replace in title. This works in order of appearance. E.g. c("\\."=" ") replaces fullstops with spaces.}

\item{assume.abbrev}{A logical whether to assume words without vowels are abbrevations and should be forced capitalised}

\item{exeption.words}{A character vector of words with case to be forced (for example abbreviations and roman numerals)}
}
\description{
This function puts character vectors into title appropriate capitalization.
}
\examples{
shakespeare_plays =c("all's well that ends well", "as you like it","comedy of errors","love's labour's lost", "measure for measure", "merchant of venice","merry wives of windsor","midsummer night's dream","much ado about nothing","taming of the shrew", "tempest", "twelfth night","two gentlemen of verona", "winter's tale", "henry iv, part i","henry iv, part ii", "henry v", "henry vi, part i","henry vi, part ii", "henry vi, part iii", "henry viii","king john", "pericles","richard ii", "richard iii", "antony and cleopatra","coriolanus","cymbeline", "hamlet","julius caesar", "king lear", "macbeth (the scottish play)", "othello", "romeo and juliet","timon of athens", "titus andronicus", "troilus and cressida")
exception.words= c("I", "II", "III", "IV", "V", "VI") # Pass in exceptions

bio_capitalize(shakespear_plays, exception.words)
}
\keyword{titles}
