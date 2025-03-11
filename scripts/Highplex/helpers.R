library(tidyverse)

## Returns protein | preds INLA model for highplex data (no fixed effects)
# df: dataframe with one row per spot assumes size columns are named size_rna & size_prot
# W: neighborhood matrix calculated based on df
# protein: character of length 1
# neighbors: boolean, false gives only spot to spot effects of the GMRF
# family: character of length 2. Specifies the likelihoods used for RNA and protein models
hpInla = function(df, W, protein, preds, neighbors = TRUE, family = c("poisson", "poisson")){
  k = length(preds)
  if(k == 0){
    # unconditional case
    m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
  } else if(k == 1){
    # conditional on 1 gene
    m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
    rnaform <- as.formula(paste(paste(preds, "~"), "1"))
    l.car <- inla(update(rnaform,.~. + f(idx, model = m)), data = df, family = family[1], offset = log(size_rna))
    if(neighbors){
      m <- inla.CCAR.model(W = W, alpha.min = 0, alpha.max = 1, phi = l.car$summary.random$idx$mean)
    } else {
      m <- inla.spotCCAR.model(W = W, alpha.min = 0, alpha.max = 1, phi = l.car$summary.random$idx$mean)
    }
  } else {
    # conditional on k genes
    X = kronecker(diag(rep(1,k)),rep(1, nrow(df)))
    mdat = data.frame("rna" = unname(unlist(as.vector(df[, names(df) %in% preds]))), 
                      "idx" = 1:(k*nrow(df)), 
                      "size" = rep(df$size_rna, k))
    mdat = cbind(mdat, X)
    names(mdat)[4:ncol(mdat)] = paste("x", 1:k, sep = "")
    rnaform = as.formula(paste("rna ~", paste(names(mdat)[4:(ncol(mdat))], collapse= "+"), "-1"))
    
    m <- inla.MCAR.model(W = W, k = k, alpha.min = 0, alpha.max = 1)
    m.car <- inla(update(rnaform, .~. + f(idx, model = m)), 
                  data = mdat, family = family[1], offset = log(size))
    if(neighbors){
      m <- inla.MCCAR.model(W = W, phi = matrix(m.car$summary.random$idx$mean, ncol = k), 
                            k = k, alpha.min = 0, alpha.max = 1)
    } else{
      m <- inla.spotMCCAR.model(W = W, phi = matrix(m.car$summary.random$idx$mean, ncol = k), 
                            k = k, alpha.min = 0, alpha.max = 1)
    }
  }
  protform <- as.formula(paste(paste(protein, "~"), "1"))
  mc.car <- inla(update(protform,.~. + f(idx, model = m)), 
                 data = df, family = family[2],
                 offset = log(size_prot),
                 control.compute = list(dic = TRUE,cpo=TRUE, waic=TRUE),
                 control.predictor = list(compute = TRUE, link = 1))
  
  return(mc.car)
}
