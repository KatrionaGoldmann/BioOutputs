#' Function to switch Gene Ids
#'
#' This function switches gene (or snp) ids using biomart 
#' @param IDs list of the Ids you want to convert
#' @param IDFrom What format these IDs are in (default Ensembl)
#' @param IDTo What format you want the IDs converted to (default gene names)
#' @param mart The biomart to use. Typically, for humans you will want ensembl (default). Alternatives can be found at listEnsembl()
#' @param dataset you want to use. To see the different datasets available within a biomaRt you can e.g. do: mart = useEnsembl('ENSEMBL_MART_ENSEMBL'), followed by listDatasets(mart).
#' @param attributes list of variables you want output
#' @keywords ids
#' @import biomaRt
#' @export
#' @examples
#' IDs = c("TNF", "A1BG", "FOX3")
#' bio_geneid(IDs, IDFrom='hgnc_symbol', IDTo = 'ensembl_transcript_id')


bio_geneid <- function(IDs, 
											 IDFrom='ensembl_transcript_id', 
											 IDTo='hgnc_symbol', 
											 mart="ensembl", 
											 dataset="hsapiens_gene_ensembl", 
											 attributes=c('chromosome_name','start_position','end_position')){

	ensembl = useEnsembl(biomart=mart, dataset=dataset)
	out.ids = getBM(attributes=c(IDFrom, IDTo, attributes), 
											filters =	IDFrom, 
											values =IDs, 
											mart = ensembl)
	
	out.ids
}



