library(Seurat)
library(tidyverse)

## loc: location of the required files
## nreplicates: whether replicate 1 or 1&2 should be used
# Required files can be found at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198353
readSpotsSpleen = function(loc, nreplicates = 1){
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
  spotsdata = readSpotsSpleen("~/Documents/postdoc/MCAR/data/spots/spleen/", 2)
  
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

## Returns protein | preds INLA model
# df: dataframe with one row per spot assumes size columns are named size_rna & size_prot
# W: neighborhood matrix calculated based on df
# protein: character of length 1
# aar: character of length k that specifies the names of one-hot-encoded AARs. First value treated as reference
# neighbors: boolean, false gives only spot to spot effects of the GMRF
# family: character of length 2. Specifies the likelihoods used for RNA and protein models
spotsInla = function(df, W, protein, preds, aar, neighbors = TRUE, family = c("poisson", "poisson")){
  k = length(preds)
  naars = length(aar)
  if(k == 0){
    # unconditional case
    m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
  } else if(k == 1){
    # conditional on 1 gene
    m <- inla.LCAR.model(W = W, alpha.min = 0, alpha.max = 1)
    rnaform <- as.formula(paste(paste(preds, "~") , paste(aar[-1],collapse = "+")))
    l.car <- inla(update(rnaform,.~. + f(idx, model = m)), data = df, 
                  family = family[1], offset = log(size_rna))
    if(neighbors){
      m <- inla.CCAR.model(W = W, alpha.min = 0, alpha.max = 1, phi = l.car$summary.random$idx$mean)
    } else {
      m <- inla.spotCCAR.model(W = W, alpha.min = 0, alpha.max = 1, phi = l.car$summary.random$idx$mean)
    }
  } else {
    # conditional on k genes
    X = kronecker(diag(rep(1,k)), as.matrix(df[, names(df) %in% aar]))
    mdat = data.frame("rna" = unname(unlist(as.vector(df[, names(df) %in% preds]))), 
                      "idx" = 1:(k*nrow(df)), 
                      "size" = rep(df$size_rna, k))
    mdat = cbind(mdat, X)
    names(mdat)[4:ncol(mdat)] = paste(aar, rep(paste("_", 1:k, sep = ""), each = naars), sep = "")
    rnaform = as.formula(paste("rna ~", paste(names(mdat)[naars:(ncol(mdat))], collapse= "+"), "-1"))
    
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
  
  protform <- as.formula(paste(paste(protein, "~"), paste(aar[-1],collapse = "+")))
  mc.car <- inla(update(protform,.~. + f(idx, model = m)), 
                 data = df, family = family[2],
                 offset = log(size_prot),
                 control.compute = list(dic = TRUE),
                 control.predictor = list(compute = TRUE, link = 1))
  
  return(mc.car)
}

readSpotsBreast = function(loc){
  mmtv_gex <- Read10X_h5(paste(loc,'GSE198353_mmtv_pymt_GEX_filtered_feature_bc_matrix.h5', sep = ""))
  mmtv_adt <- read.csv(paste(loc, 'GSE198353_mmtv_pymt_ADT.csv.gz', sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
  mmtv_image <- Read10X_Image(paste(loc, 'spatial', sep = ""))
  mmtv <- CreateSeuratObject(mmtv_gex, assay = "RNA", project = "MMTV")
  mmtv_adt <- CreateSeuratObject(mmtv_adt, assay = "CITE", project = "MMTV")
  mmtv@assays$CITE <- mmtv_adt@assays$CITE
  mmtv$nCount_CITE <- mmtv_adt$nCount_CITE
  mmtv$nFeature_CITE <- mmtv_adt$nFeature_CITE
  mmtv_image@key <- "A"
  mmtv@images <- list(A = mmtv_image)
  coords = GetTissueCoordinates(mmtv)
  coords$spot = rownames(coords)
  rownames(coords) = NULL

  # spot missing, either Fibroblast-high or low so take the one with highest prevalence
  aar <- read.csv(paste(loc, 'GSE198353_mmtv_pymt.csv', sep = ""), header = T) %>%
    add_row(Barcode = "CCAGTTCGGTAACTCA-1", AARs = "Fibroblast-high") %>%
    mutate(spot = Barcode,
           fbh = ifelse(AARs == "Fibroblast-high", 1, 0),
           fbl = ifelse(AARs == "Fibroblast-low", 1, 0),
           mac2 = ifelse(AARs == "Mac2-enriched", 1, 0),
           unknown = ifelse(AARs == "Unknown", 1, 0),
           lymph = ifelse(AARs == "Lymphocyte-enriched", 1, 0),
           mac1 = ifelse(AARs == "Mac1-enriched", 1, 0)
           ) %>%
    select(!Barcode)

    return(list("RNA" = as.matrix(mmtv@assays$RNA@counts), 
                "Protein" = as.matrix(mmtv@assays$CITE@counts), 
                "AAR" = aar,
                "coords" = coords))
}

# returns a dataframe with one spot per row with all proteins and genes specified by genes argument
SpotsCancerData = function(loc, genes){
  breast = readSpotsBreast("~/Documents/postdoc/MCAR/data/spots/cancer/")
  rna_size = unname(colSums(breast$RNA)) / median(colSums(breast$RNA))
  protein_size = unname(colSums(breast$Protein)) / median(colSums(breast$Protein))
  prot = as.data.frame(t(breast$Protein)) %>%
    mutate(spot = rownames(.))
  names(prot) = str_replace_all(names(prot), "[\\-|\\s|/|\\-|\\.|\\(|\\)]*", "")
  rna = as.data.frame(t(breast$RNA[rownames(breast$RNA)[which(rownames(breast$RNA) %in% genes)],])) %>%
    mutate(spot = rownames(.))
  names(rna) = str_replace_all(names(rna), "-", "_")
  df = full_join(breast$AAR, breast$coords, by = "spot") %>%
    full_join(prot, by = "spot") %>%
    full_join(rna, by = "spot") %>%
    mutate(size_prot = protein_size,
           size_rna = rna_size)
  return(df)
}
