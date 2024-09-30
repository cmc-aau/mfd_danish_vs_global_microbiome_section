### Setup env
library(tidyverse)

setwd("/mfd_danish_vs_global_microbiome_section")

### Create individual maps for each MFDO1 category for 16S fragment data

## Load data
data <- data.table::fread('output/2024-05-10_MFD-FL-metadata-subset.csv')

## Subset based on lacking information and geographical location
metadata <- data %>%
  filter(!is.na(mfd_hab1),
         !is.na(longitude),
         latitude <= 58.5 & latitude >= 54,
         longitude <= 19 & longitude >= 8)

## Reduce size and create "complex" correponding to full MFDO1 string
metadata.sub <- ampvis.sub$metadata %>%
  select(longitude,latitude,fieldsample_barcode,mfd_sampletype:mfd_hab3) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "),
         across(complex, ~factor(.))) %>%
  rename(long = longitude,
         lat = latitude) %>%
  mutate(group = 1)

## Create group summary based on MFDO1 categories
group.summary <- metadata.sub %>%
  group_by(complex) %>%
  reframe(samples = n())

## Extract MFDO1 categories
levels.complex <- metadata.sub %>%
  pull(complex) %>% 
  droplevels() %>%
  levels()

## Create empty data list
var.list = list()

## Change data to list format based on the MFDO1 categories
for (i in 1:length(levels.complex)) {
  
  var.list[[i]] = metadata.sub[metadata.sub$complex == levels.complex[i],]
  
}

## Create empty plot list
plot.list = list()

## Create plot for each MFDO1 category
for (i in 1:length(var.list)) {
  p = mapDK::mapDK() +
    theme_bw() + 
    ggtitle(str_c(levels.complex[i], "\n", var_list[[i]] %>% nrow(), sep = ", samples = ")) +
    geom_point(data = var_list[[i]], fill = "red", pch = 21) + 
    theme(legend.position = "right", 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size =20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20, face = "bold"), 
        plot.title = element_text(size=20))
  
  plot.list[[i]] = p
}


### Save as multi-page pdf, with each page representing af MFDO1 category
pdf("output/map_MFDO1_groups-FL16S.pdf")
for (i in 1:length(var_list)) {
  print(plot.list[[i]])
}
dev.off()
