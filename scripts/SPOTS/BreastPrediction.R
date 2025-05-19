library(INLA)
library(scico)
library(foreach)

source("./scripts/SPOTS/helpers.R")
source("./scripts/utils.R")
source("./INLA/LCAR.R")
source("./INLA/MCAR.R")
source("./INLA/indepMCAR.R")
source("./INLA/CCAR.R")
source("./INLA/MCCAR.R")
source("./INLA/spotCCAR.R")
source("./INLA/spotMCCAR.R")

breast = readSpotscancerlist("./spots/cancer/") # set to data location
aar = c("fbh", "fbl", "mac2", "unknown", "lymph", "mac1")

# Podoplanin example
df = data.frame("spot" = dimnames(breast$Protein)[[2]], 
                "prot" = unname(breast$Protein[which(rownames(breast$Protein) == "Podoplanin"),]), 
                "size_prot" = unname(colSums(breast$Protein)) / median(colSums(breast$Protein)), 
                "size_rna" = unname(colSums(breast$RNA)) / median(colSums(breast$RNA)),
                "idx" = 1:ncol(breast$Protein)) %>%
  full_join(., breast$coords, by = "spot") %>%
  full_join(breast$AAR, by = "spot") %>%
  select(!AARs)


rna = data.frame(t(breast$RNA), row.names = c()) %>% select(all_of(c("Pdpn", "Sparc")))
names(rna) = paste("rna", 1:ncol(rna), sep = "")
df = cbind(df, rna)

set.seed(1); pred_idx = sample(1:nrow(df), 200)
df_pred = df[pred_idx,]
df$prot[pred_idx] = NA

W = matrix(0, nrow(df), nrow(df))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    cp = crossprod(c(df$imagerow[i], df$imagecol[i]) - c(df$imagerow[j], df$imagecol[j]))
    if(cp > 0 && cp < 55){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

my.cluster <- parallel::makeCluster(5)
parallel::clusterEvalQ(my.cluster, {
  source('./scripts/parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:(sum(str_detect(names(df), "^rna[0-9]*$"))+1)) %dopar% {
  if(i == 1){
    spotsInla(df, W, "prot", character(), aar = aar, neighbors = F)
  } else{
    spotsInla(df, W, "prot", paste("rna",1:(i-1),sep=""), aar = aar, neighbors = F)
  }
}
parallel::stopCluster(cl = my.cluster)

displayResults(models, df, df_pred, T)

## To get sigma^2 and pi on original scale (example)
marg.stdev <- inla.tmarginal(function(tau) exp(-tau),
                             models[[1]]$marginals.hyperpar[[2]])

marg.alpha <- inla.tmarginal(function(alpha) plogis(alpha),
                             models[[1]]$marginals.hyperpar[[1]])
inla.zmarginal(marg.stdev)
inla.zmarginal(marg.alpha)

