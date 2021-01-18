# New Phytologist Supporting Information

## Article title: <br/> Pollen sterols are associated with phylogeny and environment but not with pollinators

### Authors:

Pengjuan Zu*, Hauke Koch, Orlando Schwery, Samuel Pironon, Charlotte Phillips, Ian Ondo, Iain W. Farrell, W. David Nes, Elynor Moore, Geraldine A. Wright, Dudley I. Farman, Philip C. Stevenson

**Article acceptance date:** 14 January 2021

The repository contains:

> **Inputs**: input files for each analysis (part of Table S1 in Supporting Information)

-> **NewestSterolTaxa_nov0419.csv** Initial list of taxa to get phylogeny for

—> **PollenSterolPercentPCA.csv** Relative sterol amount (percentage) for each plant species

—> **PollenTaxoFam.csv** Family names for monophyly-testing

—> **sterolCombo_new.csv** Sterol data for phylogenetic signal and plot


> **Scripts:** R scripts for the analyses

-> **MRM.R** Codes for MRM analysis (multiple regression on distance matrices).

-> **Niche.R** Clean records of environmental variables for each species, produce 3 main PCs and calculate alpha-shaped 3D volume

—> **PollenSterolPCA.R** PCA/factor analysis

—> **PollenTree.R** Building phylogenetic tree of our plant species

—> **SterolTreeDating.R** Time-calibrating the phylogenetic tree

—> **Sterol_PhyloSig_Plot.R** Testing for phylogenetic signal and plotting Fig. 2
