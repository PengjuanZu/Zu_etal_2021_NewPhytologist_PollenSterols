
# load library ecodist
> library(ecodist)

# read data with absolute pollen sterol contents
> new_sterol_matrix_forMRM <- read.csv(“~/Desktop/pollen_sterol/SterolRedux122Ladder_absolute_for_distance.csv”)

# data formatting
> bray_matrix <- new_sterol_matrix_forMRM[,2:ncol(new_sterol_matrix_forMRM)]

# calculate Bray-Curtis distance matrix
> bray_curtis_sterol <- distance(bray_matrix, method=“bray-curtis”)

# load environmental principal component data (PC1-PC3)
>PC1 <- read.csv(“~/Desktop/pollen_sterol/PC1.csv”)
>PC2 <- read.csv(“~/Desktop/pollen_sterol/PC2.csv”)
>PC3 <- read.csv(“~/Desktop/pollen_sterol/PC3.csv”)

# read phylogenetic distance matrix
> new_phyl_dist <- read.csv(“~/Desktop/pollen_sterol/PollenSterol_PhyloDist_original.csv”)

# format as distance matrix
> phyl_dist <- as.dist(new_phyl_dist)

# run MRM analysis
> sterol.mrm <- MRM(bray_curtis_sterol~phyl_dist+dist(PC1)+dist(PC2)+dist(PC3))
> sterol.mrm



