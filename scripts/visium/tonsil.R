library(INLA)
library(tidyverse)
library(foreach)
library(Seurat)

source("./scripts/SPOTS/helpers.R")
source("./scripts/utils.R")
source("./scripts/utils.R")
source("./INLA/LCAR.R")
source("./INLA/MCAR.R")
source("./INLA/CCAR.R")
source("./INLA/MCCAR.R")

data_dir = "./data/visium"  # set to location of data

data <- Read10X(data.dir = paste(data_dir, "raw_feature_bc_matrix", sep = ""))
tonsil = CreateSeuratObject(counts = data$`Gene Expression`)
tonsil[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)
seruat_img <- Read10X_Image(paste(data_dir, "spatial", sep = ""))
tonsil@images <- list(tonsil = seruat_img)

# extract data
coords = GetTissueCoordinates(tonsil, scale = "lowres") %>% 
  mutate(spot = rownames(.))
#prot = as.matrix(tonsil@assays$Protein@counts)
prot = RNA = GetAssayData(tonsil, assay = "Protein", layer = "counts")
#RNA = as.matrix(tonsil@assays$RNA@counts) # Does not work in V5.x.x
RNA = GetAssayData(tonsil, assay = "RNA", layer = "counts")

# compute sizes
rna_size = unname(colSums(RNA)) / median(colSums(RNA))
protein_size = unname(colSums(prot)) / median(colSums(prot))

# name fix
rownames(prot) = str_replace(rownames(prot), "^(.*)\\.1$", "\\1")
colnames(tpm) = str_replace(rownames(prot), "^(.*)\\.1$", "\\1")
rownames(RNA) = str_replace(rownames(RNA), "^(.*)\\.1$", "\\1")

## Find top genes, CD3E example
df = as.data.frame(t(prot)) %>%
  select("CD3E") %>%
  mutate(spot = rownames(.),
         rna = unname(RNA[which(rownames(RNA) == "CD3E"),]),
         size_prot = protein_size,
         size_rna = rna_size) %>%
  right_join(coords, by = "spot") %>%
  mutate(idx = 1:nrow(.))

## var selection (top 250 genes)
top_rna = names(sort(apply(RNA, 1, sum), decreasing = T)[1:250])
mdat = as.data.frame(t(RNA)) %>%
  select(all_of(top_rna)) %>%
  mutate(spot = rownames(.),
         prot = unname(prot[which(rownames(prot) == "CD3E"),]),
         size = protein_size) %>%
  right_join(coords, by = "spot") %>%
  mutate(idx = 1:nrow(.))
names(mdat) = str_replace(names(mdat), "-", "_")


## Compute neighborhood matrix
W = matrix(0, nrow(mdat), nrow(mdat))
for(i in 1:nrow(W)){
  for(j in i:nrow(W)){
    cp = crossprod(c(mdat$x[i], mdat$y[i]) - c(mdat$x[j], mdat$y[j]))
    if(cp <= 35 && cp > 0){
      W[i,j] = 1
      W[j,i] = 1
    }
  }
}

## find the top genes from the regression coefficents
protform = as.formula(paste("prot ~ ", paste(names(mdat)[!(names(mdat) %in% c("spot", "prot", "size", "idx", "x", "y"))], collapse= "+")))
m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
prot.car <- inla(update(protform, . ~. + f(idx, model = m)), data = mdat, family = "poisson", offset = log(size))
top_preds = prot.car$summary.fixed %>% mutate(mean = abs(mean)) %>% arrange(desc(mean)) %>% rownames
top_preds = top_preds[which(!(top_preds %in% c("(Intercept)")))]

## above procedure for reported proteins
varselection = list("CD163" = c("CD163", "SH3BGRL3", "TRAC", "BCL11A", "SLC38A2"),
                    "VIM" = c("VIM", "IGFBP7", "SPINK5", "KRT6A", "CHCHD2"),
                    "CD19" = c("CD19","BANK1", "MS4A1", "SMIM14", "TRAC"),
                    "CD4" = c("CD4", "IFNAR1", "FXYD5", "DRAM2","HERPUD1"),
                    "CD8A" = c("CD8A", "CD22", "STK17A", "IL7R", "LRMP"),
                    "MS4A1" = c("MS4A1","EVI2B","MS4A1","TMSB10","BANK1"),
                    "ACTA2" = c("ACTA2", "IGFBP7", "DSTN", "C3", "REL"),
                    "CD3E" = c("CD3E","TRBC2","TRBC1","CD22","BCL11A"))


idx = 8 # CD3E
df = data.frame("spot" = dimnames(prot)[[2]], 
                "prot" = unname(prot[which(rownames(prot) == names(varselection)[idx]),]), 
                "rna1" = unname(RNA[which(rownames(RNA) == varselection[[idx]][1]),]),
                "rna2" = unname(RNA[which(rownames(RNA) == varselection[[idx]][2]),]),
                "rna3" = unname(RNA[which(rownames(RNA) == varselection[[idx]][3]),]),
                "rna4" = unname(RNA[which(rownames(RNA) == varselection[[idx]][4]),]),
                "rna5" = unname(RNA[which(rownames(RNA) == varselection[[idx]][5]),]),
                "size_prot" = protein_size, 
                "size_rna" = rna_size) %>%
  right_join(., coords, by = "spot") %>%
  mutate("idx" = 1:nrow(.))


## prediction square
pred_idx = which(df$x >= 200 & df$x <= 300 & df$y <= 400 & df$y >= 300)
df_pred = df[pred_idx,]
df$prot[pred_idx] = NA


## estimates the models prot | G genes, G = 0,1,...,5
cores = 6
my.cluster <- parallel::makeCluster(cores)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:cores) %dopar% {
  cat("iter = ", i, "\n", sep = "")
  if(i == 1){
    try(hpInla(df, W, "prot", character(), neighbors = F, family = c("poisson", "poisson")))
  } else{
    try(hpInla(df, W, "prot", paste("rna",1:(i-1),sep=""), neighbors = F, family = c("poisson", "poisson")))
  }
}
parallel::stopCluster(cl = my.cluster)

displayResults(models, df, df_pred, T)