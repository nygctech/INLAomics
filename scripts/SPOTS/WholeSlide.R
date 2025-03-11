library(INLA)
library(tidyverse)
library(foreach)

source("./scripts/SPOTS/helpers.R")
source("./scripts/highplex/helpers.R")
source("./scripts/utils.R")
source("./INLA/LCAR.R")
source("./INLA/MCAR.R")
source("./INLA/CCAR.R")
source("./INLA/MCCAR.R")

## genes are obtained via 
# dat = predData("F4-80", NULL, spots, aar, geneindex = 1:200, genepair = "Adgre1")
# dat$top_preds
# cd4
genes = c("Cd4", "Ms4a1", "Ccl21a","Sh3bgrl3","H2-Q7" )
# f480
genes = c("Adgre1","Vcam1","Ltb","Rps13","Igfbp7")

# load both sections
spotsdata = readSpotsSpleen("~/Documents/postdoc/MCAR/data/spots/spleen/", 2)

# calculate sizes jointly (i.e. not slidewise)
prot = cbind(spotsdata$protein1, spotsdata$protein2)
RNA = cbind(spotsdata$RNA1, spotsdata$RNA2)
rna_size = unname(colSums(RNA)) / median(colSums(RNA))
protein_size = unname(colSums(prot)) / median(colSums(prot))

# Create a dataframe for slide 1 with the genes matching the proteins
rna_1 = as.data.frame(t(spotsdata$RNA1)) %>%
  select(all_of(genes)) %>%
  mutate(spot = colnames(spotsdata$RNA1),
         size_rna = rna_size[1:nrow(.)])

prot_1 = as.data.frame(t(spotsdata$protein1)) %>%
  select("F4-80") %>%
  mutate(spot = colnames(spotsdata$protein1),
         size_prot = protein_size[1:nrow(.)])

df_1 = full_join(rna_1, prot_1, by = "spot") %>%
  full_join(spotsdata$coords1, by = "spot") %>%
  full_join(spotsdata$AAR1, by = "spot")

# Create a dataframe for slide 2 with the genes matching the proteins
rna_2 = as.data.frame(t(spotsdata$RNA2)) %>%
  select(all_of(genes)) %>%
  mutate(spot = colnames(spotsdata$RNA2),
         size_rna = rna_size[(nrow(rna_1)+1):length(rna_size)])

prot_2 = as.data.frame(t(spotsdata$protein2)) %>%
  select("F4-80") %>%
  mutate(spot = colnames(spotsdata$protein2),
         size_prot = protein_size[(nrow(rna_1)+1):length(protein_size)])

df_2 = full_join(rna_2, prot_2, by = "spot") %>%
  full_join(spotsdata$coords2, by = "spot") %>%
  right_join(spotsdata$AAR2, by = "spot")

# move replicate 1 to the right of replicate 2
df_2$imagecol = df_2$imagecol + (max(df_1$imagecol)-min(df_2$imagecol)) + 5

df = rbind(df_1, df_2) %>% mutate(idx = 1:nrow(.)) %>% select(!AARs)
names(df)[1:5] = paste("rna", 1:5, sep = "")
names(df)[8] = "prot"

pred_idx = which(df$imagecol > 520)
df_pred = df[pred_idx,]
df$prot[pred_idx] = NA

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
# in parallel
my.cluster <- parallel::makeCluster(6)
parallel::clusterEvalQ(my.cluster, {
  source('parlibs.R')
})
doParallel::registerDoParallel(cl = my.cluster)
models = foreach(i = 1:6) %dopar% {
  if(i == 1){
    try(spotsInla(df, W, "prot", character(), aar = aar, neighbors = T, family = c("poisson", "poisson")))
    #try(hpInla(df, W, "prot", character(),  neighbors = T, family = c("poisson", "poisson"))) # without annotations
  } else{
    try(spotsInla(df, W, "prot", paste("rna",1:(i-1),sep=""), aar = aar, neighbors = T, family = c("poisson", "poisson")))
    #try(hpInla(df, W, "prot", paste("rna",1:(i-1),sep=""), neighbors = T, family = c("poisson", "poisson"))) # without annotations
  }
}
parallel::stopCluster(cl = my.cluster)

displayResults(m, df, df_pred)

# sequential
m0 = try(spotsInla(df, W, "prot", character(), aar = aar, neighbors = T, family = c("poisson", "poisson")))
m4 = try(spotsInla(df, W, "prot", paste("rna",1:4,sep=""), aar = aar, neighbors = T, family = c("poisson", "poisson")))

# visualize (example)
df_pred$lambda0 = m0$summary.fitted.values$mean[pred_idx]
df_pred$lambda1 = m4$summary.fitted.values$mean[pred_idx]

color_range = range(df_pred %>% pivot_longer(cols = c(prot, lambda0, lambda1)) %>% pull(value))
colour_breaks <- c(0, 10, 30, 60, 90, 140, 180) #f480
colours <- c("#112a47","#1e497a", "#3487c2", "#ffed9e", "#fc4e03", "#c4271b","#800b03")

dat_text = df_pred %>% pivot_longer(cols = all_of(c("lambda0", "lambda1"))) %>%
  mutate(name = factor(name, levels = c("lambda0", "lambda1"))) %>%
  group_by(name) %>%
  summarise(rmse = sqrt(mean((prot-value)^2))) %>%
  ungroup() %>%
  mutate(label = paste("RMSE = ", round(rmse, 1), sep = ""))

df_pred %>%
  pivot_longer(cols = c(prot, lambda0, lambda1)) %>%
  mutate(name = factor(name, levels = c("prot", "lambda1", "lambda0"))) %>%
  ggplot() + 
  geom_point(aes(x = imagecol, y = imagerow, color = value)) +
  facet_grid(. ~ name, labeller = labeller(name = c("prot" = "CD4", 
                                                    "lambda1" = "CD4 | 4 genes\n(RMSE = 15.6)", 
                                                    "lambda0" = "CD4 | 0 genes\n(RMSE = 20.5)"))) +
  theme_void() +
  scale_color_gradientn(
    limits  = color_range,
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = color_range), 1),
  ) +
  #geom_text(data = dat_text, mapping = aes(x = Inf, y = Inf, label = label),hjust= 1.02, vjust = 1.25, size = 3.5) +
  theme(text = element_text(size = 14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(0.8, "lines"),
        legend.position="bottom") +
  labs(color = "Expression")
