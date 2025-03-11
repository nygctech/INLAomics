library(INLA)
library(tidyverse)
library(foreach)

source("./scripts/Highplex/helpers.R")
source("./scripts/utils.R")
source("./scripts/utils.R")
source("./INLA/LCAR.R")
source("./INLA/MCAR.R")
source("./INLA/CCAR.R")
source("./INLA/MCCAR.R")

# downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213264
loc = ".../highplexdata/" 

## human skin
protein = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*skin_prot.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!c(X, unmapped))

## human thymus
protein = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*thymus_prot.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!c(X, unmapped))

## human tonsil
protein = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*humantonsil_prot.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!c(X, unmapped))

## mouse kdieny
protein = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*kidney_prot.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!c(X, unmapped))

names(protein) = str_replace(names(protein), "^(C{1}D{1}[0-9]+[A-z]*)\\.+.*$", "\\1") %>% 
  str_replace(., "[..]*[A-Z]{7,}", "")

protein$size = rowSums(protein[, 1:(ncol(protein)-2)]) / median(rowSums(protein[, 1:(ncol(protein)-2)]))

## human skin
rna = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*skin_RNA.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!X)

## human thymus
rna = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*thymus_RNA.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!X)

## human tonsil
rna = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*humantonsil_RNA.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!X)

## mouse kidney
rna = read.table(paste(loc, list.files(loc)[str_detect(list.files(loc), "^.*kidney_RNA.*$")], sep = ""), header = T,sep = "\t") %>%
  mutate(x = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 1)),
         y = as.numeric(str_extract(X, "^([0-9]+)x([0-9]+)", 2))) %>%
  select(!X)

names(rna) = str_replace(names(rna), "^CD(.*)$", "Cd\\1") # Avoid name clashses 
rna_top = names(sort(colSums(rna[,!colnames(rna) %in% c("x", "y", "size")]), T)[1:250])
rna$size = rowSums(rna[, 1:(ncol(rna)-2)]) / median(rowSums(rna[, 1:(ncol(rna)-2)]))

### EXAMPLE -- CD107a###

## Step 1: select genes
mdat = full_join(rna %>% select(all_of(rna_top), x, y), 
                 protein %>% select(x,y,CD107a,size), 
                 by = c("x","y"), suffix = c("rna", "prot")) %>%
  mutate(id = 1:nrow(.))

coords = cbind(mdat$x, mdat$y)
W = matrix(0, nrow(coords), nrow(coords))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    x = coords[i,]
    if(crossprod(coords[i,] - coords[j,]) == 1){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

protform = as.formula(paste("CD107a~ ", paste(names(mdat)[1:length(rna_top)], collapse= "+")))
m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
prot.car <- inla(update(protform, . ~. + f(id, model = m)), 
                 data = mdat, family = "poisson", offset = log(size))

top_preds = prot.car$summary.fixed %>% mutate(mean = abs(mean)) %>% arrange(desc(mean)) %>% rownames
top_preds = top_preds[which(!(top_preds %in% c("(Intercept)")))][1:5] # top 5 genes

# LAMP1 is pair, remove and take all top_preds if not to be included
df = full_join(rna %>% select(all_of(c("x", "y", c("LAMP1",top_preds[1:4]), "size"))), 
               protein %>% select(x,y,CD107a,size), 
               by = c("x","y"), suffix = c("rna", "prot")) %>%
  mutate(idx = 1:nrow(.))

# for ease of estimation in parallel
names(df) = c("x", "y", paste("rna", 1:5, sep = ""), "size_rna", "prot", "size_prot", "idx")

# nbhd matrix
coords = cbind(df$x, df$y)
W = matrix(0, nrow(coords), nrow(coords))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    x = coords[i,]
    if(crossprod(coords[i,] - coords[j,]) == 1){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

set.seed(12); pred_idx =sample(1:nrow(df), as.integer(nrow(df)*0.1))
pred_df = df[pred_idx,]
df$prot[pred_idx] = NA

cores = 6 # this will estimate models up to conditional G = 5, set to lower if convergence issue
my.cluster <- parallel::makeCluster(cores)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:cores) %dopar% {
  cat("iter = ", i, "\n", sep = "")
  if(i == 1){
    try(hpInla(df, W, "prot", character(), neighbors = T, family = c("poisson", "poisson")))
  } else{
    try(hpInla(df, W, "prot", paste("rna",1:(i-1),sep=""), neighbors = T, family = c("poisson", "poisson")))
  }
}
parallel::stopCluster(cl = my.cluster)

displayResults(models, df, df_pred)
