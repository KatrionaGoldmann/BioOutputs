#' Calculate the module score for expression data
#'
#' This function calculates the mean expression (z-score) of sample types
#' @param exp Expression Data
#' @param module.list List of genes in each module
#' @return data frame with expression levels for modules
#' @keywords module, svd, expression
#' @export

bio_modscores <- function(exp, module.list){
  # check how many exp genes are in module.list
  exp = t(exp)
  checkgenes <- lapply(module.list, function(x) {
    colnames(exp)[colnames(exp) %in% x]
  })
  
  modulelengths <- sapply(checkgenes, length)
  # Remove those modules with 0 genes in exp
  module.list <- module.list[!(modulelengths==0)]
  namelist <- colnames(exp[,1:ncol(exp)-1])
  
  module.output3 <- sapply(module.list, function(x) {
    data <- as.matrix(exp[,colnames(exp) %in% x])  # Catagorises the data into gene groups
    colsd <- apply(data, 2, sd)
    if (sum(colsd==0)>0) print(colnames(data)[colsd==0])
    scaledata <- scale(data[, colsd!=0])
    
    #scaledata <- scaledata[,colSums(scaledata != 0) != 0]
    svd1 <- svd(t(scaledata), nu=0, nv=3)
    #svd1 <- svd(t(data), nu=0, nv=3)       ### I CHANGED THIS LINE FROM THE ABOVE LINE TO EXPRESS THE SCALED BLOOD BETTER.
    ret <- svd1$v[,1:3]
    sign1 <- sign(sum(cor(ret[,1],  data[,colsd!=0] )))
    if (sign1 != 0)  ret <- sign1* ret
    dimnames(ret) <- list(rownames(data), paste0("ALL", 1:3))
    ret[,1]
  })
  t(module.output3)
}

