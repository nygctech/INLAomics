# SPOTS
## Breast cancer
The script `BreastPrediction.R` contains a specific example

```
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
```
which produces the data frame (first 6 rows)

|spot|prot|size_prot|size_rna|idx|imagerow|imagecol|fbh|fbl|mac2|unknown|lymph|mac1|rna1|rna2|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|AAACAAGTATCTCCCA-1|4858|1.00966284153675|0.597958382410679|1|347.98578496|426.71801312|0|0|1|0|0|0|0|10|
|AAACACCAATAACTGC-1|10214|2.02897627766279|1.44706190289229|2|405.47986128|141.20260784|0|1|0|0|0|0|1|22|
|AAACAGGGTCTATATT-1|2700|1.53248502810674|0.248135060855909|3|333.6789128|119.60900576|0|0|0|0|1|0|0|1|
|AAACAGTGTTCCTGGG-1|4985|1.11128801146314|1.81520743358199|4|488.47749232|225.00000192|1|0|0|0|0|0|1|33|
|AAACATGGTGAGAGGA-1|6983|2.51966247412832|1.35689045936396|5|424.22986144|75.88862624|0|0|0|1|0|0|0|19|
|AAACATTTCCCGGATT-1|3745|0.757743132524218|0.631069231775946|6|414.18839216|410.27843952|1|0|0|0|0|0|0|11|

In the dataframe, prot are the measurements of Podoplanin, likewise rna1 and rna2 are the measurements of Pdpn and Sparc respectively. 
To produce different dataframes, changes should be made to `rownames(breast$Protein) == "Podoplanin")` and `select(all_of(c("Pdpn", "Sparc"))`
accordingly. For example changing `"Podoplanin"` to `"CD117"` and `c("Pdpn", "Sparc")` to `c("B2m", "Rps16", "Rpl37a")` produces


|spot|prot|size_prot|size_rna|idx|imagerow|imagecol|fbh|fbl|mac2|unknown|lymph|mac1|rna1|rna2|rna3|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|AAACAAGTATCTCCCA-1|57|1.00966284153675|0.597958382410679|1|347.98578496|426.71801312|0|0|1|0|0|0|13|5|9|
|AAACACCAATAACTGC-1|64|2.02897627766279|1.44706190289229|2|405.47986128|141.20260784|0|1|0|0|0|0|26|16|8|
|AAACAGGGTCTATATT-1|23|1.53248502810674|0.248135060855909|3|333.6789128|119.60900576|0|0|0|0|1|0|4|3|0|
|AAACAGTGTTCCTGGG-1|67|1.11128801146314|1.81520743358199|4|488.47749232|225.00000192|1|0|0|0|0|0|16|28|19|
|AAACATGGTGAGAGGA-1|58|2.51966247412832|1.35689045936396|5|424.22986144|75.88862624|0|0|0|1|0|0|23|16|13|
|AAACATTTCCCGGATT-1|24|0.757743132524218|0.631069231775946|6|414.18839216|410.27843952|1|0|0|0|0|0|5|7|10|