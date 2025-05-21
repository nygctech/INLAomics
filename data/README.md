## Breast cancer
To read the data
```
source(url("https://raw.githubusercontent.com/nygctech/INLAomics/refs/heads/main/scripts/SPOTS/helpers.R"))
# loc points to a folder as specified in the main `README`.
breast = readSpotscancerlist(loc)
df = SpotsCancerData(breast, c(""))

df %>% 
  ggplot() + 
  geom_point(aes(x = imagerow, y = imagecol, color = AARs), size = 0.5) +
  labs(color = "") +
  theme_bw() +
  theme(text = element_text(size = 14), 
        legend.title.align=0.2,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("") + ylab("") + 
  coord_flip()
```
![github-small](https://github.com/nygctech/INLAomics/blob/main/data/breastcancer.png)
