# glmnet predict function
p_glmnet = function(object,newdata){
  val = predict(object, newdata, 
                type="response", 
                s = object$tuneValue$lambda, 
                alpha = object$tuneValue$alpha
  )
  as.numeric(val)
}

# gbm predict function
p_gbm = function(object, newdata){
  gbm::predict.gbm(object, newdata=newdata, 
                   type="response", 
                   n.trees = object$n.trees)
}


# shap plot function
shap_plot <- function(object, newdata, pred_func){
  ex <- explain(object, X = newdata, 
                pred_wrapper = pred_func)
  
  plot = autoplot(ex, type = "contribution")
  plot$data = plot$data[plot$data$Shapley != 0, ]
  plot$mapping$colour="black"
  plot = plot + 
    scale_fill_distiller(palette = "PiYG") + 
    scale_color_manual(values = "black")+
    theme_bw() 
  
  o = lapply(colnames(newdata), #rownames(plot$data), 
             function(feature){
               shap_dep <- data.frame(
                 value = newdata[, feature, drop = TRUE],
                 scaled_value = scale(newdata[, feature, drop = TRUE]),
                 shap = ex[, feature, drop = TRUE],
                 gene=feature
               )
             })
  o = do.call("rbind", o)
  o = o[o$shap != 0, ]
  
  shap_plot = ggplot(o, aes(x = shap, y = gene, color = value)) + 
    geom_point() + 
    labs(x="SHAP value (impact on model)", color="Feature value", y="Feature")+ 
    scale_color_gradient(low = "royalblue", high = "#f7145c", 
                         labels=c("low", "high"), breaks=range(o$value)) +
    theme_bw()
  
  return(list("explain"=ex, "shap"=shap_plot, "contribution"=plot))
}


# create force plots for all new data, save as html convert to png 
path = "./"
lapply(1:nrow(newdata), function(j){
  p1 = force_plot(object = outputs$explain[j, ], feature_values=newdata[j, ], display = "html") 
  write(p1, paste0(path, "temp", j, ".html"))
  webshot(paste0(path, "temp", j, ".html"), 
          paste0(path, "temp", j, ".png"), delay=5, 
          vwidth=3000, vheight=1)
  file.remove(paste0(path, "temp", j, ".html"))
})

# read and combine png plots
files = list.files(path, pattern='temp', full.names = T)
files = files[grepl("png", files)]
plots <- lapply(files, function(x){
  img <- as.raster(readPNG(x))
  file.remove(x)
  p = as_ggplot(rasterGrob(img, interpolate = FALSE)) + 
    theme(plot.background = element_rect(colour = "black", fill=NA, size=1))
})

ggarrange(plotlist = plots, nrow=length(plots))
