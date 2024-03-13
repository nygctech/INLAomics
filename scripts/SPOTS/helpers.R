library(Seurat)
library(tidyverse)

fig_pairs = list("CD3" = "Cd3e",
                 "F480" = "Adgre1",
                 "CD163" = "Cd163",
                 "CD29" = "Itgb1",
                 "CD68" = "Cd68",
                 "IgM" = "Ighm",
                 "CD38" = "Cd38",
                 "MadCAM1" = "Madcam1",
                 "EpCAM" = "Epcam",
                 "CD11b" = "Itgam",
                 "CD105" = "Eng",
                 "CD31" = "Pecam1",
                 "CD20" = "Ms4a1",
                 "CD169" = "Siglec1",
                 "IgD" = "Ighd",
                 "CD4" = "Cd4",
                 "CD8" = "Cd8a",
                 "CD19" = "Cd19",
                 "B220" = "Ptprc")

## loc: location of the required files
## nreplicates: whether replicate 1 or 1&2 should be used
# Required files can be found at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198353
readSpots = function(loc, nreplicates = 1){
  spleen_data <- Read10X_h5(paste(loc, "GSE198353_spleen_rep_1_filtered_feature_bc_matrix.h5", sep = ""))
  spleen_img <- Read10X_Image(paste(loc,'spatial',sep = ""))
  spleen <- CreateSeuratObject(spleen_data$`Gene Expression`, assay = "RNA", project = "spleen")
  spleen_cite <- CreateSeuratObject(spleen_data$`Antibody Capture`, assay = "CITE", project = "spleen")
  
  spleen@assays$CITE <- spleen_cite@assays$CITE
  spleen$nCount_CITE <- spleen_cite$nCount_CITE
  spleen$nFeature_CITE <- spleen_cite$nFeature_CITE
  spleen_img@assay <- c("RNA", "CITE")
  spleen_img@key <- "spleen"
  spleen@images <- list(spleen = spleen_img)
  SpatialDimPlot(spleen)
  coords1 = GetTissueCoordinates(spleen)
  
  # This is a manual fix of a spot that does not align with that of Figure 1 in Spots paper
  coords1[which((coords1$imagerow > 250) & (coords1$imagecol > 500)),1:2] = c(235.7749, 379.6919 + (379.6919-373.2788))
  coords1$spot = rownames(coords1)
  
  prot1 = as.matrix(spleen@assays$CITE@counts)
  RNA1 = as.matrix(spleen@assays$RNA@counts)
  
  aar1 <- read.csv(paste(loc, 'GSE198353_spleen_rep_1.csv', sep = ""), header = T)
  aar1[nrow(aar1)+1, ] = c("ATCATGGACTACCGAC-1", "Red pulp")
  aar1 = aar1 %>%
    mutate(spot = Barcode,
           pulp = ifelse(AARs == "Red pulp", 1, 0),
           bf = ifelse(AARs == "B follicle", 1, 0),
           mz = ifelse(AARs == "Marginal zone", 1, 0),
           pals = ifelse(AARs == "PALS", 1, 0)) %>% 
    select(!Barcode)
  
  if(nreplicates == 1){
    return(list("RNA" = RNA1, "protein" = prot1, "AAR" = aar1, "coords" = coords1))
  } else {
    spleen_data <- Read10X_h5(paste(loc, 'GSE198353_spleen_rep_2_filtered_feature_bc_matrix.h5', sep = ""))
    spleen_img <- Read10X_Image(paste(loc,'spatial2', sep=""))
    spleen <- CreateSeuratObject(spleen_data$`Gene Expression`, assay = "RNA", project = "spleen")
    spleen_cite <- CreateSeuratObject(spleen_data$`Antibody Capture`, assay = "CITE", project = "spleen")
    
    spleen@assays$CITE <- spleen_cite@assays$CITE
    spleen$nCount_CITE <- spleen_cite$nCount_CITE
    spleen$nFeature_CITE <- spleen_cite$nFeature_CITE
    spleen_img@assay <- c("RNA", "CITE")
    spleen_img@key <- "spleen"
    spleen@images <- list(spleen = spleen_img)
    SpatialDimPlot(spleen)
    coords2 = GetTissueCoordinates(spleen)
    coords2$spot = rownames(coords2)
    
    prot2 = as.matrix(spleen@assays$CITE@counts)
    RNA2 = as.matrix(spleen@assays$RNA@counts)
    
    aar2 <- read.csv(paste(loc, 'GSE198353_spleen_rep_2.csv', sep = ""), header = T)
    aar2 = aar2 %>%
      mutate(spot = Barcode,
             pulp = ifelse(AARs == "Red pulp", 1, 0),
             bf = ifelse(AARs == "B follicle", 1, 0),
             mz = ifelse(AARs == "Marginal zone", 1, 0),
             pals = ifelse(AARs == "PALS", 1, 0)) %>% 
      select(!Barcode)
    
    return(list("RNA1" = RNA1, "RNA2"=RNA2, "protein1" = prot1, "protein2" = prot2, 
                "AAR1" = aar1, "AAR2" = aar2, "coords1" = coords1, "coords2" = coords2))
  }
}

SpotsProteinData = function(loc, genepairs){
  spotsdata = readSpots("~/Documents/postdoc/MCAR/data/spots/spleen/", 2)
  
  # calculate sizes jointly (i.e. not slidewise)
  prot = cbind(spotsdata$protein1, spotsdata$protein2)
  RNA = cbind(spotsdata$RNA1, spotsdata$RNA2)
  rna_size = unname(colSums(RNA)) / median(colSums(RNA))
  protein_size = unname(colSums(prot)) / median(colSums(prot))
  
  # Create a dataframe for slide 1 with the genes matching the proteins
  rna_1 = as.data.frame(t(spotsdata$RNA1)) %>%
    select(unname(sapply(genepairs, c))) %>%
    mutate(spot = colnames(spotsdata$RNA1),
           size_rna = rna_size[1:nrow(.)])
  
  prot_1 = as.data.frame(t(spotsdata$protein1)) %>%
    mutate(spot = colnames(spotsdata$protein1),
           size_prot = protein_size[1:nrow(.)])
  
  df_1 = full_join(rna_1, prot_1, by = "spot") %>%
    full_join(spotsdata$coords1, by = "spot") %>%
    full_join(spotsdata$AAR1, by = "spot")
  
  df_1 = full_join(rna_1, prot_1, by = "spot") %>%
    full_join(spotsdata$coords1, by = "spot") %>%
    full_join(spotsdata$AAR1, by = "spot")
  
  # Create a dataframe for slide 2 with the genes matching the proteins
  rna_2 = as.data.frame(t(spotsdata$RNA2)) %>%
    select(unname(sapply(genepairs, c))) %>%
    mutate(spot = colnames(spotsdata$RNA2),
           size_rna = rna_size[(nrow(rna_1)+1):length(rna_size)])
  
  prot_2 = as.data.frame(t(spotsdata$protein2)) %>%
    mutate(spot = colnames(spotsdata$protein2),
           size_prot = protein_size[(nrow(rna_1)+1):length(protein_size)])
  
  df_2 = full_join(rna_2, prot_2, by = "spot") %>%
    full_join(spotsdata$coords2, by = "spot") %>%
    right_join(spotsdata$AAR2, by = "spot")
  
  # move replicate 1 to the right of replicate 2
  df_2$imagecol = df_2$imagecol + (max(df_1$imagecol)-min(df_2$imagecol)) + 5
  
  df = rbind(df_1, df_2) %>% mutate(idx = 1:nrow(.)) %>% select(!AARs)
  names(df)[str_detect(names(df), "-")] = c("F480", "CD169", "B220")
  
  return(df)
}

# Returns protein | preds INLA model
spotsInla = function(df, W, protein, preds, aar = c("pulp", "bf", "mz", "pals")){
  k = length(preds)
  if(k == 1){
    m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
    rnaform <- as.formula(paste(paste(preds, "~") , paste(aar[2:4],collapse = "+")))
    l.car <- inla(update(rnaform,.~. + f(idx, model = m)), data = df, 
                  family = "poisson", offset = log(size_rna))
    
    m <- inla.CCAR.model(W = W, alpha.min = 0, alpha.max = 1, phi = l.car$summary.random$idx$mean)
  } else {
    X = kronecker(diag(rep(1,k)), as.matrix(df[, names(df) %in% aar]))
    mdat = data.frame("rna" = unname(unlist(as.vector(df[, names(df) %in% preds]))), 
                      "idx" = 1:(k*nrow(df)), 
                      "size" = rep(df$size_rna, k))
    mdat = cbind(mdat, X)
    names(mdat)[4:ncol(mdat)] = paste(aar, rep(paste("_", 1:k, sep = ""), each = 4), sep = "")
    rnaform = as.formula(paste("rna ~", paste(names(mdat)[4:(ncol(mdat))], collapse= "+"), "-1"))
    
    m <- inla.MCAR.model(W = W, k = k, alpha.min = 0, alpha.max = 1)
    m.car <- inla(update(rnaform, .~. + f(idx, model = m)), 
                  data = mdat, family = "poisson", offset = log(size))
    
    m <- inla.MCCAR.model(W = W, phi = matrix(m.car$summary.random$idx$mean, ncol = k), 
                          k = k, alpha.min = 0, alpha.max = 1)
  }
  
  protform <- as.formula(paste(paste(protein, "~"), paste(aar[2:4],collapse = "+")))
  mc.car <- inla(update(protform,.~. + f(idx, model = m)), 
                 data = df, family = "poisson",
                 offset = log(size_prot),
                 control.compute = list(dic = TRUE))
  
  return(mc.car)
}