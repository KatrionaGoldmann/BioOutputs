#' Gene Locations
#'
#' This function produces ggplot and plotly plots showing the genes in a 
#' genomic range
#' @param chromosome Which chromosome are you looking at 
#' @param xrange chromosome position
#' @param subset_gene List of genes to plot/subset to within xrange
#' @param ggtheme ggplot theme to give plots. 
#' @param output_plots Logical whether to output plots
#' @importFrom ggplot2 aes lims
#' @importFrom ggbio autoplot theme_null
#' @importFrom plotly ggplotly layout config
#' @importFrom gginnards extract_layers
#' @importFrom ensembldb genes
#' @importFrom ggrepel geom_text_repel
#' @import EnsDb.Hsapiens.v75
#' @export
#' @examples
#' bio_gene_locations()

bio_gene_locations <- function(chromosome = 6, 
                          xrange = c(33.5e6, 34e6), 
                          subset_genes=c(), 
                          ggtheme=theme_null(), 
                          output_plots=TRUE){
  
  # subset to genes of interest if input
  if(length(subset_genes) > 0) { 
    TX =  data.frame(genes(EnsDb.Hsapiens.v75, 
                           filter = ~(symbol %in% subset_genes)))
  } else{ TX = data.frame(genes(EnsDb.Hsapiens.v75)) }
  TX = TX[! is.na(TX$start), ]
  TX = TX[TX$seqnames == chromosome, ]
  TX = TX[TX$start > xrange[1], ]
  TX = TX[TX$end < xrange[2], ]
  TX$direction = TX$start - TX$end
  
  if(! output_plots){
    return(TX)
  } else {
  
  # Create two plots: one stacked and one reduced
  p_txdb_stacked = autoplot(EnsDb.Hsapiens.v75,
                                   ~ (gene_id %in% TX$gene_id &
                                        symbol %in% TX$symbol),
                                   names.expr="gene_name", heights=1)
  
  p_txdb_reduced = autoplot(EnsDb.Hsapiens.v75,
                                   ~ (gene_id %in% TX$gene_id &
                                        symbol %in% TX$symbol),
                                   names.expr="",
                                   heights=1, stat="reduce")
  
  # get the gene information from the stacked plot
  gt = extract_layers(attr(p_txdb_stacked, 'ggplot'), "GeomText")
  df.p = gt[[1]]$data
  
  # output for the reduced plot
  g.plot =  attr(p_txdb_reduced, 'ggplot') + ggtheme
  
  # create the stacked names
  df.p = df.p[order(df.p$stepping, df.p$midpoint), ]
  df.p$length = rep(c(-50/length(unique(df.p$stepping)),
                      -100/length(unique(df.p$stepping))),
                    ceiling(nrow(df.p)/2))[1:nrow(df.p)]
  
  df.p =df.p[! duplicated(df.p$.labels), ]
  
  # ggplot gene location
  gglocation = g.plot + 
    geom_text_repel(data=df.p, max.overlaps=35, 
                    aes(x = midpoint, y = 1.25, label = .labels), 
                    ylim=c(1.25, Inf), size=3) + 
    lims(y=c(0.75, 1.7), x=xrange) + ggtheme
  
  # create annotations
  ann = lapply(seq_along(df.p$tx_id), function(j) {
    list(x = df.p$midpoint[j], y = 1.27, ax = 1, ay = 1,
         text = df.p$.labels[j], textangle = 0,
         font = list(color = "black"),
         arrowcolor = "black", xanchor = "auto", yanchor = "auto",
         arrowwidth = 1, arrowhead = 0, arrowsize = 1.5
    )
  })
  plotly_location <- ggplotly(g.plot) %>%  layout( annotations = ann) %>%
    config(edits = list(annotationTail = TRUE),
                   toImageButtonOptions = list(format = "svg")) 
  
  return(list("gglocation"=gglocation, "plotly_location"=plotly_location, 
              "ggstacked"=attr(p_txdb_stacked, 'ggplot') + 
                ggtheme + lims(x=xrange)))
  }
}