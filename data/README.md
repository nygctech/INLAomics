## Breast cancer
To read the data
```
source(url("https://raw.githubusercontent.com/nygctech/INLAomics/refs/heads/main/scripts/SPOTS/helpers.R"))
# loc points to a folder as specified in the main `README`.
breast = readSpotscancerlist(loc)
df = SpotsCancerData(breast, c(""))
```