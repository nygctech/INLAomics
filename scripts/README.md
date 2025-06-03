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
accordingly. For example changing `"Podoplanin"` to `"CD117"` and `c("Pdpn", "Sparc")` to `c("B2m", "Rps16", "Rpl37a")` produces the dataframe


|spot|prot|size_prot|size_rna|idx|imagerow|imagecol|fbh|fbl|mac2|unknown|lymph|mac1|rna1|rna2|rna3|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|AAACAAGTATCTCCCA-1|57|1.00966284153675|0.597958382410679|1|347.98578496|426.71801312|0|0|1|0|0|0|13|5|9|
|AAACACCAATAACTGC-1|64|2.02897627766279|1.44706190289229|2|405.47986128|141.20260784|0|1|0|0|0|0|26|16|8|
|AAACAGGGTCTATATT-1|23|1.53248502810674|0.248135060855909|3|333.6789128|119.60900576|0|0|0|0|1|0|4|3|0|
|AAACAGTGTTCCTGGG-1|67|1.11128801146314|1.81520743358199|4|488.47749232|225.00000192|1|0|0|0|0|0|16|28|19|
|AAACATGGTGAGAGGA-1|58|2.51966247412832|1.35689045936396|5|424.22986144|75.88862624|0|0|0|1|0|0|23|16|13|
|AAACATTTCCCGGATT-1|24|0.757743132524218|0.631069231775946|6|414.18839216|410.27843952|1|0|0|0|0|0|5|7|10|

## Spleen
### Prediction
The script `SpleenPred.R` contains the specific example of the CD4 protein which produces the dataframe used for prediction in the paper. The return object is a list where the
second element contains the top five predictors (`npreds = 5` default) among the top 200 most highly expressed genes and the genepair Cd4. 
```
spots = readSpotsSpleen("./spots/spleen/") # set to data location
aar = c("pulp", "bf", "mz", "pals")
names(spots)[2] = "Protein"
dat = predData("CD4", NULL, spots, aar, geneindex = 1:200, genepair = "Cd4")
dat$top_preds
#[1] "Cd4"      "Ms4a1"    "Ccl21a"   "Sh3bgrl3" "H2_Q7"   
```
|spot|prot|size_prot|size_rna|idx|imagerow|imagecol|pulp|bf|mz|pals|rna1|rna2|rna3|rna4|rna5|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|AAACACCAATAACTGC-1|57|1.57103825136612|1.1387048803569|1|383.46431769|175.03929504|1|0|0|0|0|7|14|7|5|
|AAACAGAGCGACTCCT-1|25|0.62568306010929|0.572529403812356|2|137.12668911|411.1914474|1|0|0|0|1|3|12|2|1|
|AAACAGCTTTCAGAAG-1|74|1.55874316939891|0.950385291334325|3|295.75605024|143.53976673|1|0|0|0|1|6|17|5|2|
|AAACAGGGTCTATATT-1|69|1.1974043715847|1.02757874814114|4|317.63596212|156.17730204|1|0|0|0|1|6|21|6|2|
|AAACCACTACACAGAT-1|84|1.59426229508197|0.707448965796945|5|76.95693144|483.62150052|1|0|0|0|0|2|16|5|1|
|AAACCGGGTAGGTACC-1|72|1.36543715846995|0.883466270109504|6|290.28607227|203.33228454|1|0|0|0|1|4|20|5|1|

If we instead want to consider the protein CD19, and choose from the top 300 highly expressed genes and gene pair Cd19 including top 10 predictors the above code is modified as follow
```
dat = predData("CD19", NULL, spots, aar, geneindex = 1:300, genepair = "Cd19", npreds = 10)
dat$top_preds
# [1] "E07Rik" "Cst3"   "mt_Co1" "Rplp2"  "Ms4a1"  "Fcer2a" "Cd19"   "Lcp1"   "Clu"    "Limd2" 
```
|spot|prot|size_prot|size_rna|idx|imagerow|imagecol|pulp|bf|mz|pals|rna1|rna2|rna3|rna4|rna5|rna6|rna7|rna8|rna9|rna10|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|AAACACCAATAACTGC-1|66|1.57103825136612|1.1387048803569|1|383.46431769|175.03929504|1|0|0|0|0|5|6|6|7|3|11|10|12|3|
|AAACAGAGCGACTCCT-1|39|0.62568306010929|0.572529403812356|2|137.12668911|411.1914474|1|0|0|0|0|2|3|1|3|3|11|5|7|1|
|AAACAGCTTTCAGAAG-1|71|1.55874316939891|0.950385291334325|3|295.75605024|143.53976673|1|0|0|0|0|5|1|1|6|4|10|6|8|1|
|AAACAGGGTCTATATT-1|124|1.1974043715847|1.02757874814114|4|317.63596212|156.17730204|1|0|0|0|0|6|5|1|6|3|7|7|11|1|
|AAACCACTACACAGAT-1|109|1.59426229508197|0.707448965796945|5|76.95693144|483.62150052|1|0|0|0|0|1|3|0|2|1|9|3|8|6|
|AAACCGGGTAGGTACC-1|100|1.36543715846995|0.883466270109504|6|290.28607227|203.33228454|1|0|0|0|0|5|8|1|4|1|10|4|10|1|