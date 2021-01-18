# Script for Pollen Sterol analyses - tree dating
# Orlando Schwery 2018

library(adephylo)
library(phytools)
library(phylosignal)
library(geiger)
library(picante)
library(phylobase)
library(ape)

# loading tree
phy <- read.tree(file= "PollenSterol123.tre")

# randomly resolve polytomies
set.seed(1337)
phy <- multi2di(phy, random = TRUE)

# compute branch length using Grafen's method (https://rdrr.io/cran/ape/man/compute.brlen.html)
phy <- compute.brlen(phy, method="Grafen", power=1)

# manually looked up clades to be dated from Zanne et al.
Angiospermae <- getMRCA(phy, tip=c("Lamium_album", "Nymphaea_macrosperma"))
Monocotyledoneae <- getMRCA(phy, tip=c("Colchicum_speciosum", "Aloe_perdita"))
Eudicotyledoneae <- getMRCA(phy, tip=c("Lamium_album", "Helleborus_foetidus"))
Superrosidae <- getMRCA(phy, tip=c("Paeonia_peregrina", "Genista_tinctoria"))
Rosidae <- getMRCA(phy, tip=c("Alcea_rosea", "Genista_tinctoria"))
Superasteridae <- getMRCA(phy, tip=c("Drosera_regia", "Lamium_album"))
Asteridae <- getMRCA(phy, tip=c("Lamium_album", "Philadelphus"))
Nymphaea <- getMRCA(phy, tip=c("Nymphaea_macrosperma", "Nymphaea_alexii"))
Spermatophyta <- getMRCA(phy, tip=c("Lamium_album", "Pinus_ayacahuite"))

# manually looked up ages from Zanne et al.
ageAngiospermae <- 243
ageMonocotyledoneae <- 171
ageEudicotyledoneae <- 137
ageSuperrosidae <- 118
ageRosidae <- 117
ageSuperasteridae <- 117
ageAsteridae <- 108
ageNymphaea <- 78.07  # from DateLife
ageSpermatophyta <- 327  # from Smith, Beaulieu & Donoghue 2009 An uncorrelated relaxed-clock analysis suggests and earlier origin for flowering plants, PNAS

# vector of all node names to be calibrated
datenodes <- c(Angiospermae, Monocotyledoneae, Eudicotyledoneae, Superrosidae, Rosidae, Superasteridae, Asteridae, Nymphaea, Spermatophyta)
# vector of all ages (same order as nodes above)
dateages <- c(ageAngiospermae, ageMonocotyledoneae, ageEudicotyledoneae, ageSuperrosidae, ageRosidae, ageSuperasteridae, ageAsteridae, ageNymphaea, ageSpermatophyta)

# use penalised likelihood to date tree
phydat <- chronopl(phy, 0.5, age.min=dateages, node=datenodes)
is.rooted(phydat)
is.binary(phydat)
is.ultrametric(phydat)

# plot to see where calibrated nodes are
plot(phydat, cex=0.3)
nodelabels(node=datenodes, cex=0.3)

# make vector of node ages, in my from present
H<-nodeHeights(phydat)
H<-max(H)-H
## display root age
H[phydat$edge==125][1]


# compare ages to constraints
datedtimes <- c()
for (i in 1:length(datednodes)) {
  datedtimes <- c(datedtimes, H[phydat$edge==datednodes[i]][1])
}
cbind(dateages, datedtimes)

# drop another taxon that has issues
phy <- drop.tip(phydat, tip=c("Alcea_nudiflora"))

# ladderize tree
phylad <- ladderize(phy)

#write dated tree to file
write.tree(phylad, file= "PollenSterol122DatedLadder.tre")
