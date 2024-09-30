### Setup env
library(tidyverse)
library(vegan)
library(ampvis2)
library(iNEXT)
library(patchwork)
library(rstatix)

### Import color palettes
source('scripts/MFD_colors.R')

setwd("/mfd_danish_vs_global_microbiome_section")


### Format data
## Import ASV data aggregated to OTUs (98.7% similarity ~ Species representatives)
fl.asv <- data.table::fread("data/2024-06-11_MFD_FL16S_aggregated_OTU.csv", sep = ",", header = TRUE) %>%
  rename(OTU = V1) %>%
  select(-MFD10339) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Summarise read counts per sample
read.count <- fl.asv %>%
  select(where(is.numeric)) %>%
  colSums() %>%
  data.frame("reads" = .) %>%
  rownames_to_column(var = "fieldsample_barcode")

## Import metdata file and create "complex" correponding to full MFDO1 string
metadata <- readxl::read_excel('data/2024-02-13_mfd_db.xlsx') %>%
  filter(fieldsample_barcode %in% colnames(fl.asv)) %>%
  select(-c(project_id, habitat_typenumber, sitename)) %>%
  relocate(coords_reliable, .after = "longitude") %>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  left_join(read.count) %>%
  mutate(across(mfd_sampletype, ~factor(., levels = c("Soil", "Sediment", "Water"))),
         across(mfd_areatype, ~factor(., levels = c("Natural", "Subterranean", "Agriculture", "Urban")))) %>%
  arrange(mfd_sampletype, mfd_areatype)


### Filter data
## Create summary of read count per MFDO1 level
read.count.summary <- metadata %>%
  group_by(complex) %>%
  summarise(effort = mean(reads))

## Select samples based on minimum read count (<10,000 reads)
min.reads <- read.count %>%
  filter(reads < 10000) %>%
  pull(reads) %>%
  max()

## Create summary on MFDO1 before filtering
groups.summary <- metadata %>%
  group_by(complex) %>%
  summarise(size = n())

## Create summary on MFDO1 after filtering
groups.summary.subset <- metadata %>%
  filter(reads >= min.reads) %>%
  group_by(complex) %>%
  summarise(size = n())

## Remove groups with less than 10 representative samples
groups.summary.filt <- groups.summary.subset %>%
  filter(size >= 10)

## Create metadata subset based on read counts
metadata.subset <- metadata %>%
  filter(reads >= min.reads) %>%
  mutate(across(complex, ~factor(., levels = names(mfdo1.palette)))) %>%
  filter(!is.na(complex))

## Filter metadata based on representative samples
metadata.filt <- metadata.subset %>%
  filter(complex %in% c(groups.summary.filt %>% pull(complex)))

### Write metadata files to output
data.table::fwrite(metadata.subset, "output/2024-05-10_MFD-FL-metadata-subset.csv", col.names = FALSE)
data.table::fwrite(metadata.filt, "output/2024-05-10_MFD-FL-metadata-filtered.csv", col.names = FALSE)

## Pull fieldsample_barcode from metadata subset
samples.subset <- metadata.subset %>%
  pull(fieldsample_barcode)

## Pull fieldsample_barcode from filtered metadata
samples.filt <- metadata.filt %>%
  pull(fieldsample_barcode)

## Filter observational table based on subsetting
fl.asv.subset <- fl.asv %>%
  select(any_of(samples.subset), OTU, Kingdom:Species) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Filter observational table based on filtering
fl.asv.filt <- fl.asv %>%
  select(any_of(samples.filt), OTU, Kingdom:Species) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Evaluate differences based on subsetting and filtering
nrow(fl.asv)-nrow(fl.asv.subset) # difference of 294 species
nrow(fl.asv)-nrow(fl.asv.filt) # difference of 9,872 species

### Write metadata files to output
data.table::fwrite(fl.asv.subset, "output/2024-05-10_MFD_FL16S_OTU_subset.csv", col.names = FALSE)
data.table::fwrite(fl.asv.filt, "output/2024-05-10_MFD_FL16S_OTU_filtered.csv", col.names = FALSE)

### Perform random subsample without replacement
## Random subsample without replacement for subsetted data
set.seed(123)
fl.asv.subset.ra <- fl.asv.subset %>%
  column_to_rownames(var = "OTU") %>%
  select(where(is.numeric)) %>%
  t() %>%
  rrarefy(., sample = min.reads) %>%
  t() %>%
  data.frame() %>%
  filter(rowSums(across(where(is.numeric)))!=0)

## Random subsample without replacement for filtered data
set.seed(123)
fl.asv.filt.ra <- fl.asv.filt %>%
  column_to_rownames(var = "OTU") %>%
  select(where(is.numeric)) %>%
  t() %>%
  rrarefy(., sample = min.reads) %>%
  t() %>%
  data.frame() %>%
  filter(rowSums(across(where(is.numeric)))!=0)


### Reduce color palettes
## Pull levels from different metadata files
levels <- metadata %>%
  pull(complex) %>%
  levels()

levels.subset <- metadata.subset %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

levels.filt <- metadata.filt %>%
  pull(complex) %>%
  droplevels() %>%
  levels()

## Filter color palettes
mfdo1.palette.subset <- mfdo1.palette[levels.subset]
mfdo1.palette.filt <- mfdo1.palette[levels.filt]

## Evaluate dropped levels in the different sets
setdiff(levels, levels.subset)
setdiff(levels.subset, levels)

setdiff(levels, levels.filt)
setdiff(levels.filt, levels)


### Write metadata files to output
data.table::fwrite(fl.asv.subset.ra, "output/2024-05-10_MFD_FL16S_OTU_subset-ra.csv", col.names = FALSE)
data.table::fwrite(fl.asv.filt.ra, "output/2024-05-10_MFD_FL16S_OTU_filtered-ra.csv", col.names = FALSE)


### Estimate DK species richness

## Format data for iNEXT (change to presence/absence format)
data.total <- fl.asv %>%
  column_to_rownames(var = "OTU") %>%
  select(-c(Kingdom:Species)) %>%
  mutate(across(everything(), ~+as.logical(.x)))

list <- lst()

list[["data"]] <- data.total

## Run iNEXT for all data
## Endpoint of extrapolation is by default 2x sample size (n)
set.seed(123)
object.iNEXT.total <- iNEXT(list, q = 0, datatype = "incidence_raw")

## Inspect iNEXT summary
## Observed and estimates of total richness is given with 95% CI 
## for total as well as Shannon and Simpson diversity
object.iNEXT.total

data.table::fwrite(as.data.frame(object.iNEXT.total), "output/2024-05-10_MFD-FL-iNEXT-richness-total.csv")

## Create function to format data for iNEXT by MFDO1 categories
format.iNEXT <- function(df, levels, metadata) {
  data.iNEXT <- list()
  
  abund <- df %>%
    mutate(across(everything(), ~+as.logical(.x)))
  
  for (i in 1:length(levels)) {
    filter <- metadata %>%
      filter(complex == levels[i]) %>%
      pull(fieldsample_barcode)
    
    data.filt <- abund %>%
      select(any_of(filter)) %>%
      filter(rowSums(across(where(is.numeric)))!=0) %>%
      as.matrix()
    
    data.iNEXT[[levels[i]]] <- data.filt
  }
  return(data.iNEXT)
}

## Run formatting function
data.iNEXT <- format.iNEXT(fl.asv.subset.ra, levels.subset, metadata.subset)

## Run iNEXT for MFDO1 categories
## Endpoint for extrapolation fixed at n=100
set.seed(123)
object.iNEXT <- iNEXT(data.iNEXT, q = 0, datatype = "incidence_raw", endpoint = 100)

## Inspect iNEXT summary
## Observed and estimates of total richness is given with 95% CI 
## for total as well as Shannon and Simpson diversity
object.iNEXT

data.table::fwrite(as.data.frame(object.iNEXT), "output/2024-05-10_MFD-FL-iNEXT-richness-MFDO1.csv")

## Change grouping value to factor
object.iNEXT$DataInfo$Assemblage <- factor(object.iNEXT$DataInfo$Assemblage, levels = levels.subset)


### Visualise results from iNEXT using internal iNEXT ggplot function
## Species accumulation curve per MFDO1 category
plot.iNEXT <- ggiNEXT(object.iNEXT, type = 1, color.var = "Assemblage") +
  scale_color_manual(values = mfdo1.palette.subset, breaks = levels.subset,
                     name = "MFDO1") +
  scale_fill_manual(values = mfdo1.palette.subset,
                    breaks = levels.subset,
                    name = "MFDO1") +
  ggtitle("Species accumulation curves",
          subtitle = "Interpolation and Extrapolation") +
  labs(x = "Number of sampling sites",
       y = "Number of species") +
  theme_bw(base_size = 12) +
  theme(title = element_text(face = "bold"),
        legend.position = "right",
        legend.box = "vetical")

## Remove point layer (curve end-point)
plot.iNEXT$layers[[1]] <- NULL

## Re-plot
plot.iNEXT

### save plots
png(file = 'output/MFD_FL_estimated_richness.png',
    width = 1900,
    height = 1000) 
plot.iNEXT
dev.off()

pdf(file = 'output/MFD_FL_estimated_richness.pdf',
    width = 1900,
    height = 1200)
plot.iNEXT
dev.off()

tiff(file = 'output/MFD_FL_estimated_richness.tiff',
     width = 1900,
     height = 1200)
plot.iNEXT
dev.off()

ggsave("output/MFD_FL_estimated_richness.svg", plot = plot.iNEXT, width = 19, height = 12, units = "in", dpi = "retina")


### Species abundance distribution (SAD)
## Compare the abundance distribution of V1-V8 OTUs to the log-normal distribution

## Calculate rank relative abundance of OTUs
sum.abund <- fl.asv.subset.ra %>%
  rownames_to_column(var = "OTU") %>%
  # slice_head(n = 5) %>%
  group_by(OTU) %>%
  rowwise() %>%
  reframe(SRA = sum(c_across(everything())))

## Plot rank-abundance of OTUs using ggplot
plot.rank <- sum.abund %>%
  # mutate(across(MA, ~log10(.))) %>%
  arrange(desc(SRA)) %>%
  mutate(Rank = seq(1:nrow(.))) %>%
  ggplot(aes(x = Rank, y = SRA, color = "Rank-abundance")) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("black"), limits = c("Rank-abundance")) +
  guides(color = guide_legend(title = "SAD", position = "top", reverse = TRUE)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000),
                labels = c("1", "10", "100", "1000", "10000", "100000")) +
  scale_x_continuous(breaks = seq(0, 150000, 25000)) +
  xlab('Species rank') +
  ylab('Log10(Abundance)') +
  theme_minimal(base_size = 32) +
  theme(aspect.ratio = 1, 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.x = element_line())

plot.rank

## Calculate mean relative abundance of OTUs and transform using log10
mean.abund <- fl.asv.subset %>%
  column_to_rownames(var = "OTU") %>%
  select(where(is.numeric)) %>%
  mutate(across(where(is.numeric), ~./sum(.)*100)) %>%
  rownames_to_column(var = "OTU") %>%
  # slice_head(n = 5) %>%
  group_by(OTU) %>%
  rowwise() %>%
  reframe(MRA = mean(c_across(everything()))) %>%
  mutate(across(MRA, ~log10(.)))

## Plot species abundance distribution as histogram, density curve and 
## the fitted log-normal distribution
plot.sad <- mean.abund %>%
  ggplot(aes(x = MRA)) +
  geom_histogram(aes(y = after_stat(count)), 
                 color = "black", binwidth = 0.22, alpha = 0) +
  geom_density(aes(y=after_stat(count)*0.22, 
               color = "Density"), fill = "grey", alpha = 0.25, linewidth = 1, adjust = 4) +
  stat_function(fun=function(x,mean,sd)0.22*nrow(plot.sad)*dnorm(x,mean,sd),
                args=MASS::fitdistr(plot.sad$MRA,"normal")$estimate,
                aes(color = "Fitted log-normal"), linewidth = 1, linetype = 1) +
  scale_color_manual(values = c("#2C789D", "black"), 
                     limits = c("Fitted log-normal", "Density")) +
  guides(color = guide_legend(title = "", position = "top", reverse = TRUE)) +
  scale_x_continuous(limits = c(-7,0)) +
  scale_y_continuous(limits = c(0,24000), breaks = seq(0, 24000, 4000)) +
  xlab('Log10(MRA)') +
  ylab('Number of species') +
  theme_minimal(base_size = 32) +
  theme(aspect.ratio = 1, 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.x = element_line())

plot.sad

## Arrange plots
ggarranged <- plot.rank + plot.sad + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c('A')) &
  theme(plot.tag = element_text(face = 'bold', size = 32),
        text = element_text(face = "bold"),
        legend.position = "bottom")

## Render arranged plot
ggarranged


### save plots
png(file = 'output/MFD_FL_distributions.png',
    width = 1900,
    height = 1000) 
ggarranged
dev.off()

pdf(file = 'output/MFD_FL_distributions.pdf',
    width = 1900,
    height = 1200)
ggarranged
dev.off()

tiff(file = 'output/MFD_FL_distributions.tiff',
     width = 1900,
     height = 1200)
ggarranged
dev.off()

ggsave("output/MFD_FL_distributions.svg", plot = ggarranged, width = 19, height = 12, units = "in", dpi = "retina")


### Save image
save.image(file = 'output/2024-03-19_estimated_richness.RData')