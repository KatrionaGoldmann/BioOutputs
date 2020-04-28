#' Plot Drivers of Omic Variation
#'
#' This function was adapted from dswatsons function and visualizes the 
#' strength of associations between the principal
#' components of an omic data matrix and a set of biological and/or technical
#' features.
#'
#' @param pcs PCA on expression data 
#' @param clin Data frame or matrix with rows correponding to samples and
#'   columns to technical and/or biological features to test for associations
#'   with omic data.
#' @param parametric Compute \emph{p}-values using parametric association tests?
#'   If \code{FALSE}, rank-based alternatives are used instead. See Details.
#' @param block String specifying the name of the column in which to find the
#'   blocking variable, should one be accounted for. See Details.
#' @param unblock Column name(s) of one or more features for which the 
#'   \code{block} covariate should not be applied, if one was supplied. See 
#'   Details.
#' @param kernel The kernel generating function, if using KPCA. Options include
#'   \code{"rbfdot"}, \code{"polydot"}, \code{"tanhdot"}, \code{"vanilladot"}, 
#'   \code{"laplacedot"}, \code{"besseldot"}, \code{"anovadot"}, and 
#'   \code{"splinedot"}. To run normal PCA, set to \code{NULL}. See Details.
#' @param kpar A named list of arguments setting parameters for the kernel
#'   function. Only relevant if \code{kernel} is not \code{NULL}. See Details.
#' @param top Optional number (if > 1) or proportion (if < 1) of most variable
#'   probes to be used for PCA.
#' @param n_pc Number of principal components to include in the figure.
#' @param label Print association statistics over tiles?
#' @param alpha Optional significance threshold to impose on associations. 
#'   Those with \emph{p}-values (optionally adjusted) less than or equal to 
#'   \code{alpha} are outlined in black.
#' @param p_adj Optional \emph{p}-value adjustment for multiple testing. Options
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{
#'   "bonferroni"}, \code{"BH"}, \code{"BY"}, and \code{"fdr"}. See \code{
#'   \link[stats]{p.adjust}}.
#' @param title Optional plot title.
#' @param legend Legend position. Must be one of \code{"bottom"}, \code{"left"},
#'   \code{"top"}, \code{"right"}, \code{"bottomright"}, \code{"bottomleft"},
#'   \code{"topleft"}, or \code{"topright"}.
#' @param hover Show \emph{p}-values by hovering mouse over tiles? If 
#'   \code{TRUE}, the plot is rendered in HTML and will either open in your 
#'   browser's graphic display or appear in the RStudio viewer.
#'
#' @details
#' Strength of association is measured by -log \emph{p}-values, optionally
#' adjusted for multiple testing. When \code{parametric = TRUE}, significance
#' is computed from Pearson correlation tests (for continuous features) or 
#' ANOVA \emph{F}-tests (for categorical features). When \code{parametric =
#' FALSE}, significance is computed from rank-based alternatives, i.e. Spearman 
#' correlation tests (for continuous features) or Kruskal-Wallis tests (for 
#' categorical features). 
#'
#' An optional blocking variable may be provided if samples violate the
#' assumption of independence, e.g. for studies in which subjects are observed
#' at multiple time points. If a blocking variable is identified, it will be
#' regressed out prior to testing for all variables except those explicitly
#' exempted by the \code{unblock} argument. Significance is then computed from
#' partial correlation tests for continuous data (Pearson if \code{parametric = 
#' TRUE}, Spearman if \code{parametric = FALSE}) or repeated measures ANOVA 
#' \emph{F}-tests (under rank-transformation if \code{parametric = FALSE}).
#'
#' When supplying a blocking variable, be careful to consider potential
#' confounding effects. For instance, features like sex and age are usually
#' nested within subject, while subject may be nested within other variables
#' like batch or treatment group. The \code{block} and \code{unblock} arguments
#' are designed to help parse out these relationships.
#' 
#' Numeric and categorical features are tested differently. To protect against
#' potential mistakes (e.g., one-hot encoding a Boolean variable), 
#' \code{plot_drivers} automatically prints a data frame listing the class of
#' each feature.
#'
#' If \code{kernel} is non-\code{NULL}, then KPCA is used instead of PCA. See
#' \code{\link{plot_kpca}} for more info. Details on kernel functions and their
#' input parameters can be found in \code{kernlab::\link[kernlab]{dots}}.
#' #'
#' @examples
#' library(SummarizedExperiment)
#' library(edgeR)
#' library(dplyr)
#' data(airway)
#' cnts <- assay(airway)
#' keep <- rowSums(cpm(cnts) > 1) >= 4
#' mat <- cpm(cnts[keep, ], log = TRUE)
#' clin <- colData(airway) %>%
#'   as_tibble(.) %>%
#'   select(Run, cell, dex)
#' plot_drivers(mat, clin)
#'
#' @seealso
#' \code{\link{plot_pca}}, \code{\link{plot_kpca}}
#'
#' @export
#' @importFrom limma is.fullrank
#' @importFrom purrr map_chr
#' @importFrom kernlab eig
#' @importFrom kernlab rotated
#' @import dplyr
#' @import ggplot2
#'

plot_drivers <- function(pcs,
                         clin,
                         block = NULL,
                         unblock = NULL,
                         kernel = NULL,
                         kpar = NULL,
                         top = NULL,
                         n.pc = 5L,
                         label = FALSE,
                         alpha = 0.05,
                         p.adj = NULL,
                         title = 'Variation By Feature',
                         legend = 'right',
                         hover = FALSE) {
  
  
  sig <- function(j, pc) {  # p-val fn
    mod <- lm(pcs[, pc] ~ clin[[j]])
    if_else(clin[[j]] %>% is.numeric, summary(mod)$coef[2L, 4L], 
            anova(mod)[1L, 5L])
  }
  
  # remove columns all identical or all different
  keep.cols = which(apply(clin, 2, function(x) {
    length(unique(x[! is.na(x)])) < nrow(clin) &
      length(unique(x[! is.na(x)])) > 1
  }))
  clin = clin[, keep.cols]
  
  df <- expand.grid(Feature = colnames(clin), PC = colnames(pcs)) %>%
    rowwise(.) %>%
    dplyr::mutate(Association = sig(Feature, PC)) %>%   # Populate
    ungroup(.)
  if (!p.adj %>% is.null) {
    df <- df %>% mutate(Association = p.adjust(Association, method = p.adj))
  }
  df <- df %>% 
    mutate(Significant = if_else(Association <= alpha, TRUE, FALSE), 
           Association = -log(Association))
  
  # Build plot
  if (!p.adj %>% is.null && p.adj %in% c('fdr', 'BH', 'BY')) {
    leg_lab <- expression(~-log(italic(q)))
  } else {
    leg_lab <- expression(~-log(italic(p)))
  }
  p <- ggplot(df, aes(PC, Feature, fill = Association, text = Association,
                      color = Significant)) +
    geom_tile(size = 1L, width = 0.9, height = 0.9) +
    coord_equal() +
    scale_fill_gradientn(
      colors = c('white',  'dodgerblue1', 'dodgerblue3', 'dodgerblue4'), 
      name = leg_lab) +
    scale_color_manual(values = c('grey90', 'black')) +
    guides(color = FALSE) +
    labs(title = title, x = '', y='') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(angle = 315, hjust = 0))
  if (label) {
    p <- p + geom_text(aes(label = round(Association, 2L)))
  }
  return(p)
}
