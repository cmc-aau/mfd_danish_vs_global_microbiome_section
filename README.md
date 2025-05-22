# mfd_danish_vs_global_microbiome_section
This is a repo for the Danish vs global microbiome section of the MFD paper. 

## Background
The Microflora Danica near full-length 16S rRNA dataset (V1-V8) contains 21.3 million sequences representing 141,252 species-level (98.7%) OTUs. 
The near full-length 18S rRNA dataset (V4-V9) contains 13.4 million sequences representing 12,447 species-level (99%) OTUs. 

Using only the near full-length 16S rRNA gene UMI data generated on the Nanopore platform (5.8 million reads and 101,423 98.7) we investigated how well the data captured the collective Danish terrestrial species pool. Pan-habitat rarefactions of the 16S rRNA sequences indicates that the 16S rRNA data captures Denmark’s dominant species in the investigated habitats. The taxonomic diversity measured as the number of species representatives in the Danish habitats under investigation was quantified using rarefaction (interpolation) and prediction (extrapolation) with [Hill numbers of order <em>q</em>](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/13-0133.1). Hill numbers, <em><sup>q</sup></em>Δ, differs by the sensitivity to the relative incidence of the species representative OTUs. The first three Hill numbers correspond to the species richness (<em>q</em> = 0), the exponential of Shannon entropy (<em>q</em> = 1), and the inverse Simpson concentration (<em>q</em> = 2), with the later two being referred to as Shannon and Simpson diversity. The last three scritps are needed to make the combined figure panel used in the manuscript. 


The scripts in this repository are part of the [Microflora Danica project](https://github.com/cmc-aau/mfd_wiki/wiki). 
The scripts are used to generate maps of each category of the [MFD Ontology](https://github.com/cmc-aau/mfd_wiki/wiki/Ontology) of both 16S fragments derived from metagenomic sequencing as well as based on near full-length 16S UMI sequences. As a continuation of this, the repo contains scripts for mapping to the 10 and 1 km reference grid of Denmark with subsequent spatial thinning. The last script uses the 10 km representative set and produces a list of files and read patterns to extract the 16S fragments, to be used in subsequent taxonomic classification. 

Be advised, that for the metagenomic-derived data, the term "OTU" is only used due to format requirements by ampvis2, and they do not represent classical OTUs. 
The generated profiles can be thought of as taxonomic bins. 

## Scripts
### UMI amplicon 16S data 
`scripts/maps_FL.R` generates maps of each MFDO1 category across the country based on the 16S UMI amplicon sequencing. The script outputs a multi-paged PDF. 


`scripts/richness_estimates_FL15S.R` estimates the total richness based on the 16S UMI amplicon sequencing, as well as Shannon and Simpson diversity - [not to confused with their corresponding indexes](https://johnsonhsieh.github.io/iNEXT/). Estimates are made for the total data and categoty-specific estimates based on MFDO1. 


`scripts/250502_Figure2a_Sequence_novelty.Rmd` estimates the novelty of the seqeunces from both the prokaryotic (98.7% OTUs) data and the eukaryotic (99% OTUs) data.


`scripts/250502_Figure2b_Rarefraction_curves.Rmd` creates the habitat-specific and pan-habitat rarefaction curves. 


`scripts/250501_Figure2cd_Database_coverage.Rmd` creates the summaries of the database evaluations using both the GPC (Global Prokaryotic Census) amplicons and the spatial thinned metagenomic-derived 16S fragments from this study. 


## Data
The scripts rely on data files available from the MFD Zenodo [repo](https://zenodo.org/records/12605769) and the MFD [github](https://github.com/cmc-aau/mfd_metadata), from where the original output files are also available. 

