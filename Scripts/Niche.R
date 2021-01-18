
library(raster)
library(alphashape3d)
library(factoextra)


# PCA with all climate, soil and topography layers across world

#Load all .tif layers from folder
setwd("C:/...")
grids <- list.files(pattern = "tif", full.names = T)
predictors <- stack(grids)

#Reduce the dimensionality of the environmental space to 3D using PCA
pca <- princomp(na.omit(values(predictors)), cor = TRUE)

#Get contributions of the variables for the 3 axes
var <- get_pca_var(pca)
var$coord[,1:3]

#Variance explained by 3 axes
eig <- (pca$sdev)^2
variance <- eig*100/sum(eig)
sum(variance[1:3])

#Map the 3 first axes of the PCA
PCs <- predict(predictors, pca, index = 1:3)
#Convert to Points
PCs.points <- rasterToPoints(PCs)



  
  tabFin<-c()
  
  #Beforehand, create raster files (same grid as environmental variables) for each species with 1=presence and 0=absence
  #Create a list of species occurrence raster files found in folder
  setwd(paste("C:/..."))
  splist <- list.files(pattern = "tif")
 
  #Loop for each species
  for (i in splist){
    
    #Load species occurrence raster
    sp<-raster(paste("C:/...",i,sep=""))
    #Transform to point data
    PA.data <- rasterToPoints(sp)
    #Extract presences only
    occ.data <- PA.data[which(PA.data[,3]>0),]
    
    #Extract positions on PCA axes and for each environmental variable at each occurrence point
    if (length(na.omit(occ.data))<4) {
      occur.PCs <- extract(PCs, rbind(occ.data[1:2],occ.data[1:2]))
      occur.Vars <- extract(predictors, rbind(occ.data[1:2],occ.data[1:2]))
    }else{
      occur.PCs <- extract(PCs, occ.data[,1:2])
      occur.Vars <- extract(predictors, occ.data[,1:2])
    }
    
    
    #Build Alpha-Shape 3D for species with five points or more
    #If less than 5 points, no niche size is calculated but positions are
    if (length(na.omit(occur.PCs[,1]))>=5) {
      
      #Alpha-Shape 3D with alpha defined as 2
      ashape3d.occ <- ashape3d(unique(as.matrix(na.omit(occur.PCs))), alpha = 2)
      #Extract niche characteristics for the species
      N<-length(na.omit(occur.PCs)[,1]) # number of occurrences
      vol<-volume_ashape3d(ashape3d.occ) # niche size
      pos1<-mean(na.omit(occur.PCs[,1])) # niche position axis 1
      pos2<-mean(na.omit(occur.PCs[,2])) # niche position axis 2
      pos3<-mean(na.omit(occur.PCs[,3])) # niche position axis 3
      tab<-c(N,vol,pos1,pos2,pos3)
      pos<-c() # niche position for each variable
      for (k in 1:nlayers(predictors)){
        pos<-mean(na.omit(occur.Vars[,k]))
        tab<-c(tab,pos)
      }
      tabFin<-rbind(tabFin,tab)
      
    } else {
      
      N<-length(unique(na.omit(occur.PCs)[,1])) # number of occurrences
      vol<-NA
      pos1<-mean(occur.PCs[,1]) # niche position axis 1
      pos2<-mean(occur.PCs[,2]) # niche position axis 2
      pos3<-mean(occur.PCs[,3]) # niche position axis 3
      tab<-c(N,vol,pos1,pos2,pos3)
      pos<-c() # niche position for each variable
      for (k in 1:nlayers(predictors)){
        pos<-mean(occur.Vars[,k])
        tab<-c(tab,pos)
      }
      tabFin<-rbind(tabFin,tab)
      next
      
    }
    
    
  }
  
  #format and save final table for all species
  colnames(tabFin)<-c("N", "Vol", "Pos1", "Pos2", "Pos3", names(predictors))
  rownames(tabFin)<-splist
  setwd(paste("C:/..."))
  write.csv(tabFin, file="niche.csv", row.names = TRUE, col.names = TRUE)



