# display results for a list of INLA objects
displayResults = function(InlaList, training, test, neighbors = T){
  pred_idx = as.integer(rownames(test))
  return(data.frame(
    k = sapply(InlaList, function(x){(ifelse(neighbors, (nrow(x$summary.hyperpar)-2)/2, nrow(x$summary.hyperpar)-2))}),
    dic = sapply(InlaList, function(x){x$dic$dic}),
    rmse = sapply(InlaList, function(x){sqrt(mean((test$prot - x$summary.fitted.values$mean[pred_idx])^2))}),
    rmsetrain = sapply(InlaList, function(x){sqrt(mean((training$prot[-pred_idx] - x$summary.fitted.values$mean[-pred_idx])^2))}),
    cpo = sapply(InlaList, function(x){sum(log(x$cpo$cpo[-pred_idx]))}),
    waic = sapply(InlaList, function(x){x$waic$waic})
  ))
}
