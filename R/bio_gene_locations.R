#' Gene Locations
#'
#' This function produces ggplot and plotly plots showing the genes in a 
#' genomic range
#' @param chromosome Which chromosome are you looking at 
#' @param xrange chromosome position
#' @param subset_gene List of genes to plot/subset to within xrange
#' @param ggtheme ggplot theme to give plots. 
#' @param output_plots Logical whether to output plots
#' @param ens_version The EnsemblDB version. You will need this package 
#' installed. Options include: "EnsDb.Hsapiens.v75", "EnsDb.Hsapiens.v86"...
#' @param stat The plotting stat for autoplot (either 'gene_level' for gene level
#' annotation or 'tx_level' for transcript level annotation)
#' @param label_col The column to label genes/transcripts by from the EnsemblDB. 
#'  Options include those in `EnsDb.Hsapiens.v75@tables$gene` or 
#'  `EnsDb.Hsapiens.v75@tables$tx` (depending on stat), for example: 'gene_id', 
#'  'tx_id', 'symbol', 'entrezid', ...
#' @importFrom ggplot2 aes lims
#' @importFrom ggbio autoplot theme_null
#' @importFrom plotly ggplotly layout config
#' @importFrom gginnards extract_layers
#' @importFrom ensembldb genes transcripts
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples
#' bio_gene_locations()

bio_gene_locations <- function(chromosome = 6, 
                               xrange = c(28e6, 34e6), 
                               subset_genes=c(), 
                               font_size = 3,
                               ggtheme = theme_null(), 
                               output_plots = TRUE, 
                               ens_version = "EnsDb.Hsapiens.v75", 
                               stat = "gene_level", 
                               text_col = "black",
                               annotation_col = "black", 
                               label_col = ifelse(stat == "gene_level", 
                                                  "symbol", "tx_id")){
  
  require(ens_version, character.only = TRUE)
  edb <- get(ens_version)
  
  if(! stat %in% c("gene_level", "tx_level")) 
    stop("stat must be either 'gene_level' or 'tx_level'.")
  id_col <- ifelse(stat == "gene_level", "gene_id", "tx_id")
  
  if(length(subset_genes) > 0) {
    filt <- ~ (symbol %in% subset_genes) 
  } else { filt <- AnnotationFilterList() }
  
  # subset to genes of interest if input
  if(stat == "gene_level") { 
    TX =  data.frame(genes(edb, filter = filt))
  } else {TX =  data.frame(transcripts(edb, filter = filt)) }
  TX = TX[! is.na(TX$start), ]
  TX = TX[TX$seqnames == chromosome, ]
  TX = TX[TX$start > xrange[1], ]
  TX = TX[TX$end < xrange[2], ]
  TX$direction = TX$start - TX$end
  
  if(nrow(TX) == 0) stop(ifelse(length(subset_genes) == 0, 
                                "No genes in this range.", 
                                "No subset genes in this range."))
  
  if(! output_plots){
    return(TX)
  } else {
    
    if(nrow(TX) > 50) 
      message(paste0("----------------\n\n", 
                     "Beware: This region contains many ", 
                     ifelse(stat=="gene_level", "genes", "transcripts"), 
                     " (n=", nrow(TX), ") ", 
                     "therefore plotting may take some time. Try subsetting ", 
                     "to genes of interest if this is too time consuming (set ",
                     "output_plots=FALSE to view the genes/transcripts within ", 
                     "this region in a data frame). \n\n----------------\n"))
    
    if(stat == "gene_level"){
      p_txdb = autoplot(edb,
                        ~ (gene_id %in% TX$gene_id &
                             symbol %in% TX$symbol),
                        names.expr="", heights=1, 
                        stat="reduce", color=annotation_col)
    } else { 
      p_txdb = autoplot(edb, ~ (gene_id %in% TX$gene_id), color=annotation_col,
                        names.expr="", heights=1, stat="identity")
      
    }
    
    # get the gene information from the plot
    gt = extract_layers(attr(p_txdb, 'ggplot'), "GeomText")
    df.p = gt[[1]]$data
    df.p$.labels <- TX[match(df.p[, id_col], TX[, id_col]), label_col]
    
    # get the segments
    sl = extract_layers(attr(p_txdb, 'ggplot'), "GeomSegment")
    
    # ggplot output
    g.plot =  attr(p_txdb, 'ggplot') + ggtheme
    
    # create the names and positions
    df.p = df.p[order(df.p$stepping, df.p$midpoint), ]
    if(stat == "gene_level") df.p$y <- 1 else df.p$y <- df.p$stepping
    
    # ggplot gene location
    gglocation = g.plot + 
      geom_text_repel(data=df.p, max.overlaps=35, color=text_col,
                      aes(x = midpoint, y = y + 0.25, label = .labels), 
                      ylim=c(1.25, Inf), size=font_size) + 
      lims(x=xrange) + ggtheme
    
    # create annotations
    ann = lapply(seq_along(df.p[, id_col]), function(j) {
      list(x = df.p$midpoint[j], y = df.p$y[j] + 0.25, ax = 1, ay = 1,
           text = df.p$.labels[j], textangle = 0,
           font = list(color = text_col, size=font_size*4),
           arrowcolor = text_col, xanchor = "auto", yanchor = "auto",
           arrowwidth = 1, arrowhead = 0, arrowsize = 1.5
      )
    })
    plotly_location <- ggplotly(g.plot) %>%  layout(annotations = ann) %>%
      plotly::config(edits = list(annotationTail = TRUE),
                     toImageButtonOptions = list(format = "svg")) 
    
    return(list("gglocation"=gglocation, "plotly_location"=plotly_location))
  }
}
