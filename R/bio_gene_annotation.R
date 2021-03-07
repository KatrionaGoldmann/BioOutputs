library(biomaRt)
library(easyPubMed)
library(rentrez)
library(httr)
library(pbapply)
library(enrichR)
library(splitstackshape)

enrichr_library <- read.table("https://github.com/KatrionaGoldmann/BioOutputs/blob/master/data/enrichr_libraries.csv?raw=true", sep=",", 
                              header=TRUE, stringsAsFactors = FALSE)
enrichrdbs <-  listEnrichrDbs()

enrichr_update <- enrichrdbs$libraryName[! enrichrdbs$libraryName %in% 
                                           enrichr_library$libraryName]

if(length(enrichr_update) > 0){
  warning(paste("enrichr library csv needs to be updated. Missing:", 
                paste(enrichr_update, collapse=", ")))
}
remove(enrichr_update, enrichrdbs)

# wrapper function to find all info about a list of genes
#' @param genes
#' @param verbose Whether to split out the progress
#' @param gene_types Highlight interesting genes? 
#' @param gene_summaries Whether to find a summary of the gene info
#' @param associated_diseases Whether to find the associated diseases 
#' @param publications Do you want pulication info? (This is time consuming
#' if you have a lot of genes)
#' @param disease_cutoff The DisGeNET Score cutoff
#' @param diseases The MeSH disease classes to use 
#' ("CD20"=Immune System Diseases)
#' @param publication_keywords Words to search for publications with the genes 
#' @param publication_split Whether to use AND or OR when searching the 
#' publication keywords (default is OR)
gene_summary <- function(genes, 
                         verbose=TRUE, 
                         gene_types=TRUE, 
                         gene_summaries=TRUE, 
                         associated_diseases=TRUE, 
                         publications=FALSE,
                         gene_types_list = immune_types,
                         disease_cutoff=0, 
                         diseases=c("C20", "C05", "C10", "C17"), 
                         publication_keywords=c("Rheuamtoid", "Autoimmune"), 
                         publication_split="OR"){
  df = data.frame("Gene"=genes)
  
  if(gene_types){
    print("Annotating EMR's favourites...")
    df$Type <- what_are_you(genes, gene_types_list)$Type
  }
  
  if(gene_summaries){
    print("Getting gene summaries...")
    temp <- tell_me_about_yourself(genes)
    df <- cbind(df, temp[match(df$Gene, temp$Gene), 
                         c("description", "summary")])
  }
  
  if(associated_diseases){
    print("Finding associated diseases...")
    temp = what_diseases(genes, cutoff=disease_cutoff, 
                         disease_classes = diseases, verbose=FALSE)
    df <- cbind(df, temp[match(df$Gene, temp$Gene), 2])
    colnames(df)[ncol(df)] = "Associated_diseases"
  }
  
  if(publications){
    print("Getting publications from PubMed...")
    temp <- got_any_publications(genes, 
                                 keywords = publication_keywords, 
                                 split=publication_split)
    df <- cbind(df, temp[match(df$Gene, temp$Gene), 2])
    colnames(df)[ncol(df)] = "Publications"
  }
  
  df[is.na(df)] = ""
  rownames(df) = df$Gene
  
  return(df)
}

#' returns pathways from enrichR
#' @param genes List of gene names
#' @param drop_terms terms in databases not to search. If NULL only keep_terms 
#' used
#' @param keep_terms terms in databases to search. If NULL only all except 
#' drop_terms used
#' @param dbs database of enrichR pathways. If NULL all are selected minus 
#' those removed by drop_terms. See listEnrichrDbs().  
#' @param cutoff pvalue cutoff
#' @param min_N minimum number of genes to consider
#' @param remove_old Whether to remove legacy libraries ()
#' @param libraries Which type of libraries to check for enrichment. Options 
#' include c('Transcription', 'Pathways', 'Ontologies', 'Diseases_Drugs', 
#' 'Cell_Types', 'Misc'))
#' @param plot_cap Maximum number of terms to plot
#' @param plot_pcutoff Pvalue cutoff for plotting

enriched_pathways <- function(genes, 
                              drop_terms=c("User", "Enrichr"), 
                              keep_terms=c(),
                              dbs=NULL, 
                              cutoff=0.05, 
                              min_N=2, 
                              remove_old = TRUE, 
                              libraries = c('Transcription', 'Pathways'), 
                              plot_cap=20, 
                              plot_pcutoff=0.01,
                              plot_colour = "Combined.Score", 
                              plot_x="P.value", 
                              plot_vline=NA){
  
  if(any(! libraries %in% c('Transcription', 'Pathways', 'Ontologies', 
                            'Diseases_Drugs', 'Cell_Types', 'Misc'))){
    stop(paste("libraries must be in c('Transcription', 'Pathways',", 
               "'Ontologies', 'Diseases_Drugs', 'Cell_Types', 'Misc')"))
  }
  
  
  
  if(is.null(dbs)) {
    dbs <- enrichr_library 
    pathways <- apply(data.frame(dbs[, libraries], stringsAsFactors = FALSE), 1, 
                      function(r) any(r=="x"))
    dbs <- dbs[pathways, ]
  }
  if(remove_old){
    dbs <- dbs[dbs$Legacy != "x", ]
  }
  if(length(drop_terms) > 0) {
    dbs <- dbs[! grepl(paste0(drop_terms, collapse="|"), dbs$libraryName, 
                       ignore.case = TRUE), ]
  }
  if(length(keep_terms) > 0) {
    dbs <- dbs[grepl(paste0(keep_terms, collapse="|"), dbs$libraryName, 
                     ignore.case = TRUE), ]
  }
  
  enriched <- enrichr(genes, dbs$libraryName)
  
  temp = do.call(rbind, enriched)
  temp$Library = gsub("\\..*", "", rownames(temp))
  if(length(drop_terms) > 0) {
    temp <- temp[! grepl(paste0(drop_terms, collapse="|"), temp$Term, 
                         ignore.case = TRUE), ]
  }
  if(length(keep_terms) > 0) {
    temp <- temp[! grepl(paste0(keep_terms, collapse="|"), temp$Term, 
                         ignore.case = TRUE), ]
  }
  
  temp = temp[, c('Term', 'Library',  'Overlap', 'P.value', 'Adjusted.P.value', 
                  'Old.P.value', 'Old.Adjusted.P.value', 'Odds.Ratio', 
                  'Combined.Score', 'Genes')]
  
  temp$N = as.numeric(gsub("/.*", "", temp$Overlap))
  temp = temp[temp$N >= min_N, ]
  temp = temp[temp$P.value <= cutoff, ]
  
  temp = temp[order(temp$P.value), ]
  
  if(nrow(temp)<=1) stop("Need more than one pathway")
  plot_df <- temp[temp$P.value < plot_pcutoff, ]
  plot_df[, plot_x] = as.numeric(as.character(plot_df[, plot_x]))
  plot_df = plot_df[order(plot_df[, plot_x]), ]
  plot_df = plot_df[1:(min(nrow(plot_df), plot_cap, na.rm=T)), ]
  plot_df$count = splitstackshape::getanID(plot_df$Term)$.id 
  plot_df$Term[plot_df$count != 1] = paste0(plot_df$Term, " (", plot_df$count, ")")[plot_df$count != 1]
  plot_df$Term <- factor(plot_df$Term, levels=rev(unique(plot_df$Term)))
  plot_df$x=plot_df[, plot_x]
  plot_df$col=plot_df[, plot_colour]
  
  if(grepl("P.value", plot_x)) {
    plot_df$x = -log10(plot_df$x)
    plot_x=paste0("-log(", plot_x, ")")
  }
  
  plot = ggplot(plot_df, aes(x=Term, y=x, fill=col)) +
    geom_bar(stat="identity") + 
    coord_flip() +
    theme_minimal() + 
    labs(x="", y=plot_x, fill=plot_colour) + 
    scale_fill_continuous(low="red", high="midnightblue", 
                          guide=guide_colorbar(reverse=grepl("P.value", plot_colour)))
  
  if(! is.na(plot_vline)) {
    plot <- plot + geom_hline(yintercept = plot_vline, 
                              colour="grey60", linetype="dashed")
  }
  
  temp_genes = sort(table(unlist(strsplit(temp$Genes,  ";"))), decreasing=TRUE)
  
  
  return(list(enrichment=temp, enrichment_genes_table=temp_genes, plot=plot))
}





# Immune type genes EMR are interested in 
immune_types = list("ADAM"=c("ADAM", ""),
                    "ALDOA" =c("ALDO", ""),
                    "BMP"=c("BMP", ""),
                    "CAMK"=c("CAMK", ""),
                    "CASP"=c("CASP", ""),
                    "CIDE"=c("CIDE", ""),
                    "Chemokines"=c("CCL|CXC", "signalling & communication"),
                    "Cluster of Differentiation"=c("^CD[0-9]", ""),
                    "CC"=c("^CC", ""),
                    "Collagens"=c("^COL", ""),
                    "Complement"=c('^C[0-9]', ""),
                    "Complement Receptors"=c('^CR[0-9]', ""),
                    "CREB"=c("^CRE", ""),
                    "CSF"=c("^CSF", ""),
                    "CILP"=c("^CILP", ""),
                    "C4"=c("C4", ""),
                    "Dual-specificity phosphatase"=c("^DUSP", ""),
                    "Eukaryotic initiation factors"=c("^EIF", ""),
                    "Fibroblast Growth Factors"=c("FGF", ""),
                    "FLG"=c("FLG", ""),
                    "FOS"=c("FOS", "regulators of cell proliferation, differentiation, and transformation"),
                    "GAS"=c("GAS", ""),
                    "HLA"=c("HLA", ""),
                    "HOX"=c("HOX", ""),
                    "ICAM"=c("^ICAM", ""),                    
                    "Immunoglobulins (heavy)"=c("^IGH", "Antibodies"),
                    "Immunoglobulins (lamda)"=c("^IGLV", "Antibodies"),
                    "Immunoglobulins (kappa)"=c("^IGK", "Antibodies"),
                    "Interferon induced transmembrane proteins" = c("IFITM", ""),
                    "Interferon (OAS)" = c("OAS", ""),
                    "Interferon (ISG)" = c("ISG", ""),
                    "Interferon regulatory factors" = c("IRF", ""),
                    "Interferon (IFN)" = c("IFN", ""),
                    "INF" = c("INF", ""),
                    "Interleukins and their receptors"=c("^IL[0-9]", ""),
                    "JAK" = c("^JAK", ""),
                    "KRT"=c("KRT", ""),
                    "Killer Cells"=c("KLR", "Natural killer (NK) cells are lymphocytes that can mediate lysis of certain tumor cells and virus-infected cells without previous activation"),
                    "LAM" = c("LAM", ""),
                    "LGAL" = c("LGAL", ""),
                    "MAPK"=c("^MAP[0-9]|MAPK", "kinase kinases: signalling, communication"),
                    "MMP"=c("MMP", ""),
                    "MS4"=c("MS4", ""),
                    "NAD(P)H:quinone acceptor oxidoreductase"=c("^NQO", ""),
                    "NFK"=c("NFK", ""),
                    "NLRP"=c("NLRP", ""),
                    "NOTCH"=c("NOTCH", ""),
                    "Phosphorlation Genes"=c("PHOSPHO", "Protein phosphorylation. Hydrolase superfamily"),
                    "Phosphatase"=c("^PTP", "critical regulators of signal transduction"),
                    "PDC"=c("^PDC", ""),
                    "PLEKH"=c("^PLEKH", ""),
                    "PPP"=c("^PPP", ""),
                    "PSORS"=c("PSO", ""),
                    "RAS"=c("RAS|RALA", "GTPases which act as molecular switches in the cell"),
                    "REL"=c("^REL", ""),
                    "RGS"=c("RGS", "Pain signalling"),
                    "S100"=c("^S100", "regulation of protein phosphorylation, transcription factors, Ca2+ homeostasis, the dynamics of cytoskeleton constituents, enzyme activities, cell growth and differentiation, and the inflammatory response"),
                    "SAA"=c("^SAA", ""),
                    "SLAM"=c("^SLAM", ""),
                    "SOCS"=c("^SOCS", ""),
                    "STAT"=c("^STAT", ""),
                    "TGF"=c("^TGF", ""),
                    "TLR"=c("TLR", ""),
                    "TNF"=c("TNF", ""),
                    "TP5"=c("^TP5", ""),
                    "TRAF"=c("TRAF", ""),
                    "TREX"=c("TREX", ""),
                    "Ubiquitin specific peptidases"=c("^USP", ""),
                    "Ubiquitin Family"=c("^UBF", ""),
                    "WNT"=c('WNT', ""))



#' genes EMR are interested in
#' @param genes A list of gene names
#' @param types A list of annotation type where each contains a 
#' character vector for grep command and general info about the gene
what_are_you <- function(genes, types=immune_types){ # What are you? 
  df <- data.frame("Gene"=genes)
  df$type <- ""
  for(i in 1:length(types)){
    df$Type[grepl(types[[i]][1], df$Gene)] <- names(types)[i]
    df$Function[grepl(types[[i]][1], df$Gene)] <- types[[i]][2]
  }
  return(df)
}

#'  full gene names
#'  @param genes A list of gene names
#'  @param name_type What format are the names
whats_your_name <- function(genes, name_type = "hgnc_symbol"){ 
  df <- data.frame("Gene"=genes)
  ensembl <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  df2 <- getBM(attributes=c("hgnc_symbol", "description"),
               filters = name_type, values=genes, mart=ensembl)
  
  df$Description <- df2$description[match(df$Gene, df2$hgnc_symbol)]
}


#' Gene summaries
#' @param genes A list of gene names
tell_me_about_yourself <- function(genes){
  all.df <- pblapply(genes, function(g) {
    r_search <- entrez_search(db="gene", term=g)
    if(length(r_search$ids) > 0){
      info <- entrez_summary(db="gene", id=r_search$ids, version="2.0")
      df = data.frame("Gene"=g, 
                      "name"=extract_from_esummary(info, "name"),
                      "description"=extract_from_esummary(info, "description"), 
                      "summary"=extract_from_esummary(info, "summary"),
                      stringsAsFactors = F)
      df = df[df$name == g, ]
      if(nrow(df) > 0){
        df = t(data.frame(apply(df, 2, function(x){
          x = x[x != ""]
          paste(unique(x), collapse="; ")
        })))
        colnames(df) = c("Gene", "name",  "description", "summary")
      }
    } else{
      df = data.frame("Gene"=g, "name"=g, "description"="", "summary"="")
    }
    df[is.na(df)] = ""
    return(df)
  })
  final.df = do.call(rbind, all.df)
}



#' What diseases are known to be associated with genes
#' @param genes A list of gene names
#' @param cutoff The DisGeNET Score cutoff
#' @param disease_classes The MeSH disease classes to use 
#' (C"D20"=Immune System Diseases)
#' @param verbose Whether to spit out plotgress
#' @references https://www.disgenet.org/api/#/GDA/gdaByGene
what_diseases <- function(genes, 
                          cutoff=0, 
                          disease_classes=c("C20", "C05", "C10", "C17"), 
                          verbose=T){
  df = data.frame("Gene"=genes, "Associated_diseases"="", stringsAsFactors = F)
  lapply(genes, function(g){
    if(verbose) print(g)
    cutoff_query = paste0("?min_score=", cutoff) 
    if(!is.null(disease_classes)) {
      disease_query = paste0(paste0("&disease_class=", disease_classes), 
                             collapse="") 
    }else {disease_query=""}
    call = paste0('https://www.disgenet.org/api/gda/gene/', g, 
                  cutoff_query, disease_query)
    
    get_diseases <- GET(call)
    get_diseases_text <- httr::content(get_diseases, "text", encoding = "UTF-8")
    if(get_diseases_text != ""){
      get_diseases_json <- fromJSON(get_diseases_text, flatten = FALSE)
      diseases <- unique(get_diseases_json$disease_name)
    } else{diseases=""}
    df$"Associated_diseases"[df$Gene == g] <<- paste(diseases, collapse="; ")
  })
  colnames(df)[2] = paste0("Diseases_minscore_", cutoff, "_class_", 
                           paste(disease_classes, collapse="_"))
  return(df)
}


#' What publications are available?
#' @param genes A list of gene names
#' @param keywords Words to search with the genes
#' @param split search keywords using AND or OR (default is OR)
got_any_publications <- function(genes, 
                                 keywords=c("rheumatoid", 
                                            "inflammatory", 
                                            "autoimmune"), 
                                 split="OR", 
                                 verbose=T){
  
  df <- pblapply(genes, function(g){
    my_query <- paste(g, ' AND (', 
                      paste(keywords, collapse=paste0(" ", split, " ")), ")")
    my_entrez_id <- get_pubmed_ids(my_query)
    hits <- fetch_pubmed_data(my_entrez_id)
    titles <- custom_grep(hits, "ArticleTitle", "char")
    data.frame("Gene"=g, "Publications"=paste(titles, collapse="; "))
  })
  do.call(rbind, df)
}

