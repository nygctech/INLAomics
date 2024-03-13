library(INLA)
source("./scripts/SPOTS/helpers.R")
source("./INLA/LCAR.R")
source("./INLA/MCAR.R")
source("./INLA/CCAR.R")
source("./INLA/MCCAR.R")

loc = "~/Documents/postdoc/MCAR/data/spots/spleen/" # set to folder with data
df = SpotsProteinData(loc, fig_pairs)

## Create neighborhood matrix
coordsmat = cbind(df$imagerow, df$imagecol)
W = matrix(0, nrow(coordsmat), nrow(coordsmat))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    cp = crossprod(coordsmat[i,] - coordsmat[j,])
    if(cp > 0 && cp < 45){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

########################
### Recreate CD3 row ###
########################

# Treat mRNA as known and part of the fixed effects. Effect sizes are then used
# to sort the genes in the order they are introduced to the model
PROT = 1
protform = as.formula(paste(paste(names(fig_pairs)[PROT], "~"), paste(c("bf", "mz", "pals",unname(sapply(fig_pairs, c))), collapse = "+")))
m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
prot.car <- inla(update(protform, . ~. + f(idx, model = m)), data = df, family = "poisson", offset = log(size_prot))
top_preds = prot.car$summary.fixed %>% mutate(mean = abs(mean)) %>% arrange(desc(mean)) %>% rownames
top_preds = top_preds[which(!(top_preds %in% c("(Intercept)","bf","pals","mz", fig_pairs[PROT][[1]])))]

# Protein-gene pair
preds = fig_pairs[[PROT]]
ccar = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# Protein | 2 genes
preds = c(fig_pairs[[PROT]], top_preds[1])
c2car = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# continune if dic lowers
c2car$dic$dic < ccar$dic$dic

# Protein | 3 genes
preds = c(fig_pairs[[PROT]], top_preds[1:2])
c3car = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# continune if dic lowers
c3car$dic$dic < c2car$dic$dic

# Protein | 4 genes
preds = c(fig_pairs[[PROT]], top_preds[1:3])
c4car = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# continune if dic lowers
c4car$dic$dic < c3car$dic$dic

# Protein | 5 genes
preds = c(fig_pairs[[PROT]], top_preds[1:4])
c4car = spotsInla(df, W, names(fig_pairs)[PROT], preds)

# DIC goes up after 5 so we stop.