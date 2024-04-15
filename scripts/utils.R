displayResults = function(InlaList, training, test, neighbors = T){
  #pred_idx = as.integer(rownames(test))
  return(data.frame(
    k = sapply(InlaList, function(x){(ifelse(neighbors, (nrow(x$summary.hyperpar)-2)/2, nrow(x$summary.hyperpar)-2))}),
    dic = sapply(InlaList, function(x){x$dic$dic}),
    rmse = sapply(InlaList, function(x){sqrt(mean((test$prot - x$summary.fitted.values$mean[test$idx])^2))}),
    rmsetrain = sapply(InlaList, function(x){sqrt(mean((training$prot[-test$idx] - x$summary.fitted.values$mean[-test$idx])^2))}),
    cpo = sapply(InlaList, function(x){sum(log(x$cpo$cpo[-test$idx]))}),
    waic = sapply(InlaList, function(x){x$waic$waic})
  ))
}