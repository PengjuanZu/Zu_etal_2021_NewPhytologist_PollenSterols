## load library
library(ggbiplot)
library(factoextra)
library(stats)
library(pracma)


##DATA_UPLOAD#########################################################################################################
# Uploads the .csv file containing information on the % (i.e. relative abundance) of sterols in each pollen sample

sterolRelAb<-read.table("~/PollenSterolPercentPCA.csv",sep=",", header=TRUE)
head(sterolRelAb)

pollSterol.pca<-prcomp(sterolRelAb[,c(2:26)], center=TRUE, scale.=TRUE)
summary(pollSterol.pca)


# Eigenvalues
eigenValues <- get_eigenvalue(pollSterol.pca)
 
# Keep components with EV > 1 (12)
nComponents<-length(which(eigenValues$eigenvalue>=1))

# Varimax rotation with Kaiser normalisation
rawLoadings     <- pollSterol.pca$rotation[,1:nComponents] %*% diag(pollSterol.pca$sdev, nComponents, nComponents)
rotatedLoadings <- varimax(rawLoadings, normalize = TRUE)$loadings

print(rotatedLoadings)

write.csv(rotatedLoadings, 'pollSterolFactorTable.csv')


##########################################################
