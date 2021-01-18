# Edited script to get synthetic tree of sterol taxa
# Pengjuan Zu & Orlando Schwery 21. August 2019

library(rotl)
library(ape)

mu <- read.csv("NewestSterolTaxa_nov0419.csv", stringsAsFactors = FALSE)  ##load species data
mu
taxon_search <- tnrs_match_names(names = mu$species)  # check Names
# fixed NAs
  # fix Petrocosmea cryptica to just Petrocosmea
  # Euphorbia milii splendens to without splendens
  # Hieracium umbellatum agg. to without agg.
  # cistus sp. to cistus
  # remove weird signs after Nymphaceae
  # remove trailing space after Helminthotheca echioides
  # add proper space in Parodia schumanniana

# fix fuchsia & NAs
taxon_search[48, 2] <- "Fuchsia"
taxon_search[39, 2] <- "Digitalis purpurea"
mu <- mu[-98, ] # remove duplicate Petrocosmea in data (so mu and taxon_search have same dimensions)


#### add two more columns in dataset
mu$ott_name <- taxon_search$unique_name  # add OTT names
mu$ott_id <- taxon_search$ott_id  # add OTT IDs

ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]  # get id's
tr <- tol_induced_subtree(ott_ids = ott_in_tree)  # get published and taxonomy synthetic tree subset

tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores = FALSE)  # clean tip labels from ott id's and underscores
plot(tr, cex=0.45)

write.tree(tr, "PollenSterol123.tre")
write.csv(mu, "PollenSterol123_ROTL_names.csv")
