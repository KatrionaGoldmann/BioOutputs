#' Longitudinal Sankey
#' @description Creates a sankey plot with heatmaps showing how individuals progress over time. 
#' @param samp.orders A list of named vectors for sample order at each timepoint. Vector names must correspond to matchable
#' ids. 
#' @import networkD3
#' @import htmlwidgets
#' @examples
#'library(survival)
#'data(kidney)
#'df <- kidney

#'library(lme4)
#'data(sleepstudy)
#'data = sleepstudy
#'data = data[data$Days !=1 | data$Subject != "308", ]
#'time.col="Days"
#'exp.cols = "Reaction"
#'sub.col = "Subject"
#'
#'row.order=list()
#'for(i in unique(data[, time.col])){
#'  hm = Heatmap(data[data[, time.col] == i, exp.cols], cluster_rows=cluster)
#'  row.order[[paste('timepoint', as.character(i))]] = setNames(unlist(row_order(hm)), data[data[, time.col] == i, sub.col])
#'}
#'
#'bio_sankey(samp.order=row.order)

bio_sankey = function(samp.order){

  all.subj = unique(unlist(lapply(samp.order, function(x) names(x))))
  
  # Start from zero for javascript
  for(i in 1:length(samp.order)){samp.order[[i]] = samp.order[[i]]-1  }
  
  # Order the later nodes accordingly
  for(i in  2:(length(samp.order))){samp.order[[i]] = samp.order[[i]] + max(samp.order[[i-1]]+1, na.rm=T)  }
  
  links = data.frame("source"=c(), "target"=c(), "value"=c())
  for(sub in all.subj){
    keep.order = samp.order[unlist(lapply(samp.order, function(x) sub %in% names(x)))]
    for(i in 1:(length(keep.order)-1)){
      link.sub = data.frame("source"=keep.order[[i]][sub], "target"=keep.order[[i+1]][sub], "value"=1)
      links = rbind(links, link.sub)
    }
  }

  # Add the node names
  nodes = data.frame("name"=c())
  for(i in samp.order){nodes = rbind(nodes, data.frame("name"=names(i)[order(i)]))}
  
  # create the sankey plot
  p = sankeyNetwork(Links = links, Nodes = nodes, nodePadding=0,
                     Source = "source", Target = "target",
                     NodeID = "name", Value="value", margin=0,
                     fontSize= 12, nodeWidth = 30, iterations = 0)
  
  p
}


bio_sankey2 = function(samp.order){
  
  
  all.subj = unique(unlist(lapply(samp.order, function(x) names(x))))
  init.order.w = init.order[match(all.subj, names(init.order))]
  sub.order.w = sub.order[match(all.subj, names(sub.order))]
  
  # Create order for the links
  sub.order.w = sub.order.w + length(init.order.w)
  links = data.frame(init.order.w-1, sub.order.w-1, 1, names(init.order.w), names(sub.order.w))
  names(links) = c("source", "target", "value")
  
  # Create the node names
  nodes = data.frame("name" = c(names(init.order.w)[order(as.numeric(links$source)+1)], 
                                names(sub.order.w)[order(as.numeric(links$target)+1-length(init.order.w))]))
  
  p = sankeyNetwork(Links = links, Nodes = nodes, nodePadding=0,
                    Source = "source", Target = "target",
                    NodeID = "name", Value="value", margin=0,
                    fontSize= 12, nodeWidth = 30, iterations = 0)
  
  p
}










