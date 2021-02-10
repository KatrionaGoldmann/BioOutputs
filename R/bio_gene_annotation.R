library(biomaRt)
library(easyPubMed)
library(rentrez)
library(httr)
library(pbapply)


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
what_do_you_do <- function(genes, 
                 verbose=TRUE, 
                 gene_types=TRUE, 
                 gene_summaries=TRUE, 
                 associated_diseases=TRUE, 
                 publications=FALSE, 
                 disease_cutoff=0, 
                 diseases=c("C20", "C05", "C10", "C17"), 
                 publication_keywords=c("Rheuamtoid", "Autoimmune"), 
                 publication_split="OR"){
  df = data.frame("Gene"=genes)
  
  if(gene_types){
    print("Annotating EMR's favourites...")
    df$Type <- what_are_you(genes, c(immune_types, Others))$Type
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

# Immune type genes EMR are interested in 
immune_types = list("Collagens"=c("COL", ""),
                    "HOX"=c("HOX", ""),
                    "Cluster of Differentiation"=c("^CD[0-9]", ""),
                    "Interleukins and their receptors"=c("^IL[0-9]", ""),
                    "Chemokines"=c("CCL|CXC", "signalling & communication"),
                    "Complement"=c('^C[0-9]', ""),
                    "Complement Receptors"=c('^CR[0-9]', ""),
                    "Interferon (IFI)" = c("IFI", ""),
                    "Interferon (OAS)" = c("OAS", ""),
                    "Interferon (ISG)" = c("ISG", ""),
                    "Interferon (IRF)" = c("IRF", ""),
                    "Interferon (IFN)" = c("IFN", ""),
                    "TLR"=c("TLR", ""),
                    "ADAM"=c("ADAM", ""),
                    "STAT"=c("STAT", ""),
                    "GAS"=c("GAS", ""),
                    "Fibroblast Growth Factors"=c("FGF", ""),
                    "CAMK"=c("CAMK", ""),
                    "TRAF"=c("TRAF", ""),
                    "TNF"=c("TNF", ""),
                    "BMP"=c("BMP", ""),
                    "WNT"=c('WNT', ""),
                    "MS4A"=c("MS4A", ""),
                    "KRT"=c("KRT", ""),
                    "NLRP"=c("NLRP", ""),
                    "PSORS"=c("PSO", ""),
                    "NFK"=c("NFK", ""),
                    "ALDOA" =c("ALDO", ""),
                    "CIDE"=c("CIDE", ""),
                    "FLG"=c("FLG", ""),
                    "TP5"=c("^TP5", ""),
                    "MMP"=c("MMP", ""),
                    "NOTCH"=c("NOTCH", ""),
                    "USP"=c("USP", ""),
                    "CREB"=c("CRE", ""),
                    "MAPK"=c("MAPK", "kinase kinases: signalling, communication"),
                    "HLA"=c("HLA", ""),
                    "SAA"=c("SAA", ""),
                    "Phosphorlation Genes"=c("PHOSPHO", "Protein phosphorylation. Hydrolase superfamily"),
                    "FOS"=c("FOS", "regulators of cell proliferation, differentiation, and transformation"),
                    "RAS"=c("RAS|RALA", "GTPases which act as molecular switches in the cell"),
                    "RGS"=c("RGS", "Pain signalling"),
                    "Immunoglobulins (IGH)"=c("^IGH", "Antibodies"),
                    "Immunoglobulins (IGL)"=c("^IGLV", "Antibodies"),
                    "Immunoglobulins (IGK)"=c("^IGK", "Antibodies"))

Others = list("Other noteworthy genes"=c("GSK3|HYAL|UBF|KLF|GATA|PFKFB|ACAN|S100B", ""))

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
      info <- entrez_summary(db="gene", id=r_search$ids)
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

