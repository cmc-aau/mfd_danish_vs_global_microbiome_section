# mfd_danish_vs_global_microbiome_section
This is a repo for the Danish vs global microbiome section of the MFD paper. 

The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 
The scripts are used to generate maps of each category of the [MFD Ontology](https://github.com/cmc-aau/mfd_wiki/wiki/Ontology) of both 16S fragments derived from metagenomic sequencing as well as based on FL16S sequences. As a continuation of this, the repo contains scripts for mapping to the 10 and 1 km reference grid of Denmark with subsequent spatial thinning. The last script uses the 10 km representative set and produces a list of files and read patterns to extract the 16S fragments, to be used in subsequent taxonomic classification. 

Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

## Scripts
### Amplicon 16S data 
`scripts/maps_FL.R` generates maps of each MFDO1 category across the country based on 16S amplicon sequencing (V1-V8). The script outputs a multi-paged PDF. 


`scripts/richness_estimates_FL` estimates the total richness, as well as Shannon and Simpson diversity - [not to confused with their corresponding indexes](https://johnsonhsieh.github.io/iNEXT/). Estimates are made for the total data and categoty-specific estimates based on MFDO1. 


### Metagenomic 16S data 
`scripts/maps_FL.R` generates maps of each MFDO1 category across the country based on 16S fragments derived from the metagenomic sequencing. The script outputs a multi-paged PDF. 


`scripts/grid_spatial_thinning.R` thins the data by mapping to the 10 and 1 km reference grids of Denmark. If more than one sample of the same MFDO1 category is present within the same grid cell, the samples with the overall highest (mean) similarity to the other samples is chosen. 


`scripts/grid_reads_for_classification.R` uses the spatially thinned 10 km reference dataset to identify names of the corresponding sequencing files and identified 16S fragments. Then it generates a list of patterns to search for based on a file-read key, from which unambigious reads has previosly been removed. The cpmmand used was:


`find /dir | grep -F -f DATE_libs-classification.txt -exec cat {} \; | grep -F -f DATE_forward-reads-classification.txt > DATE_grid-reads-for-classification.fq`


## Data
The scripts rely on data files available from the MFD Zenodo [repo](https://zenodo.org/records/12605769) and the MFD [github](https://github.com/cmc-aau/mfd_metadata), from where the original output files are also available. 

