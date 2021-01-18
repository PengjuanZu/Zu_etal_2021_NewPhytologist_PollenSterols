# Script for Pollen Sterol analyses
# Orlando Schwery 2018

library(ape)
library(geiger)
# library(plotrix)
library(phytools)

# load phylogeny and plot to check
phy <- read.tree(file= "PollenSterol122DatedLadder.tre")
plot(phy, cex=0.4)

# correct some tip labels
phy$tip.label[10] <- "Digitalis_purpurea"
phy$tip.label[34] <- "Campanula_fragilis"
phy$tip.label[84] <- "Fuchsia"
plot(phy, cex=0.4)

# read sterol trait data
sterolCombo <- read.csv(file= "sterolCombo_new.csv", stringsAsFactors=FALSE, row.names=1)
# row names to same format as
rownames(sterolCombo) <- gsub(" ", "_", rownames(sterolCombo))
# remove taxa from tree or sterol data that are not in respective other set
steroltrdat <- treedata(phy, sterolCombo)
# save cleaned up sterol data
write.csv(steroltrdat$dat, file = "SterolRedux122Ladder_2020_12_13.csv")

# divide up data by different type
discretes <- steroltrdat$dat[, c(4:5)]
class(discretes) <- "character"
sterolH <- steroltrdat$dat[, 6]
class(sterolH) <- "numeric"
sterolperc <- steroltrdat$dat[, 7:31]  # relative (this is used for most of the figure)
class(sterolperc) <- "numeric"
sterolabs1 <- steroltrdat$dat[, 33:57]  # absolute per mg pollen (the sum of this is used for the figure (see sterolSum below))
class(sterolabs1) <- "numeric"
sterolabs2 <- steroltrdat$dat[, 59:83]  # absolute per flower
class(sterolabs2) <- "numeric"
sterolCgrp <- steroltrdat$dat[, 85:87]  # % in C 24 substitution groups
class(sterolCgrp) <- "numeric"
sterolDeltagrp <- steroltrdat$dat[, 89:92]  # % in delta group B ring moiety
class(sterolDeltagrp) <- "numeric"
sterolSum <- steroltrdat$dat[, 58]  # sum of absolute sterol amount
class(sterolSum) <- "numeric"
# check by dimensions
dim(discretes)  # 122 2
length(sterolH)  # 122
dim(sterolperc)  # 122 25
dim(sterolabs1)  # 122 25
dim(sterolabs2)  # 122 25
dim(sterolCgrp)  # 122 3
dim(sterolDeltagrp)  # 122 4
length(sterolSum) # 122

# re-assign tree as crosschecked with data (though should be same in this case)
phy <- steroltrdat$phy
# get distance matrix
distmat <- cophenetic(phy)
write.csv(distmat, file = "PollenSterol_PhyloDist_2020_12_13.csv")

# calculate commonness (relative presence of each sterol across all plant species)
commonness <- c()
for (i in 1:ncol(sterolperc)) {
  commonness <- c(commonness, sum(sterolperc[,i] > 0))
}
commonness <- commonness/nrow(sterolperc)

# calculate abundance (average propoertion of each sterol in each plant species)
abundance <- c()
for (i in 1:ncol(sterolperc)) {
  abundance <- c(abundance, mean(sterolperc[,i]))
}

abundance <- c()
for (i in 1:ncol(sterolperc)) {
  abundance <- c(abundance, mean(sterolperc[,i]))
}

# combine both and add column names as in remaining data
richab <- rbind(abundance, commonness)
colnames(richab) <- colnames(sterolperc)
# Save Abundance and Commonness table
write.csv(richab, file = "AbundanceCommonness_2020_12_13.csv")

# convert sterol data to phylo4d objects
objdiscr <- phylobase::phylo4d(phy, discretes)
objH <- phylobase::phylo4d(phy, sterolH)
objperc <- phylobase::phylo4d(phy, sterolperc)
objabs1 <- phylobase::phylo4d(phy, sterolabs1)
objabs2 <- phylobase::phylo4d(phy, sterolabs2)
objCgrp <- phylobase::phylo4d(phy, sterolCgrp)
objDeltagrp <- phylobase::phylo4d(phy, sterolDeltagrp)
objsterolSum <- phylobase::phylo4d(phy, sterolSum)

# Phylosignal: Moran's I, Abouheif's Cmean, Pagel's Lambda, Blimberg's K, and K*
set.seed(1337)
AllPhyloSig_H <- phylosignal::phyloSignal(objH, methods = c("all"), reps = 999, W = NULL)
AllPhyloSig_perc <- phylosignal::phyloSignal(objperc, methods = c("all"), reps = 999, W = NULL)
AllPhyloSig_abs1 <- phylosignal::phyloSignal(objabs1, methods = c("all"), reps = 999, W = NULL)
AllPhyloSig_abs2 <- phylosignal::phyloSignal(objabs2, methods = c("all"), reps = 999, W = NULL)
AllPhyloSig_C <- phylosignal::phyloSignal(objCgrp, methods = c("all"), reps = 999, W = NULL)
AllPhyloSig_Delta <- phylosignal::phyloSignal(objDeltagrp, methods = c("all"), reps = 999, W = NULL)
AllPhyloSig_Sum  <- phylosignal::phyloSignal(objsterolSum, methods = c("all"), reps = 999, W = NULL)

# check correct format (list of two data frames, stat and pvalue, both 5 variables per observation)
str(AllPhyloSig_H)
str(AllPhyloSig_perc)
str(AllPhyloSig_abs1)
str(AllPhyloSig_abs2)
str(AllPhyloSig_C)
str(AllPhyloSig_Delta)
str(AllPhyloSig_Sum)

# combine all to full result dataset
Fulldata <- rbind(
  cbind(AllPhyloSig_H$stat, AllPhyloSig_H$pvalue),
  cbind(AllPhyloSig_perc$stat, AllPhyloSig_perc$pvalue),
  cbind(AllPhyloSig_abs1$stat, AllPhyloSig_abs1$pvalue),
  cbind(AllPhyloSig_abs2$stat, AllPhyloSig_abs2$pvalue),
  cbind(AllPhyloSig_C$stat, AllPhyloSig_C$pvalue),
  cbind(AllPhyloSig_Delta$stat, AllPhyloSig_Delta$pvalue),
  cbind(AllPhyloSig_Sum$stat, AllPhyloSig_Sum$pvalue)
)
# add more descriptive column names for p-values
colnames(Fulldata)[6:10] <- c("p-Cmean", "p-I", "p-K", "p-K.star", "p-Lambda")
# save full output table
write.csv(Fulldata, file = "PollenSterol_PhylosigResults_122_2020_12_13.csv")

#######################
# Prepare for plotting
#######################
# test monophyly of Families
library(MonoPhy)
# read taxon file of families
fams <- read.csv("PollenTaxoFam.csv", header=FALSE)
# run monophyly test
solution <- AssessMonophyly(phy, fams)
# summary of results, and detailed result
GetSummaryMonophyly(solution)
GetResultMonophyly(solution)
# output monophyly vs taxonomy and intruder plots
PlotMonoPhyly(solutionn, phy, plot.type="monoVStax", PDF=TRUE, PDF_filename="PollenFamMonoVsTax.pdf")
PlotMonoPhyly(solutionn, phy, plot.type="intruders", PDF=TRUE, PDF_filename="PollenFamIntruder.pdf")
# save monophyly results
write.csv(solution$Taxlevel_1$result, file="PollenFamSolution.csv")

# add function for coloured boxes aroudn clades from Liam Revell's blog http://blog.phytools.org/2017/05/simple-method-to-draw-boxes-around.html
cladebox<-function(tree,node,color=NULL,...){
   if(is.null(color)) color<-make.transparent("yellow",0.2)
     obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
   h<-max(nodeHeights(tree))
   parent<-tree$edge[which(tree$edge[,2]==node),1]
   x0<-max(c(obj$xx[node]+obj$xx[parent])/2,obj$xx[node]-0.05*h)
   x1<-obj$x.lim[2]
   dd<-getDescendants(tree,node)
   y0<-min(range(obj$yy[dd]))-0.5
   y1<-max(range(obj$yy[dd]))+0.5
   polygon(c(x0,x1,x1,x0),c(y0,y0,y1,y1),col=color,
       border=0)
}

# vectors with names in correct tree order

# fams with ntips >1
# mono_famnums <- c(19, 10, 4, 6, 15, 13, 20, 23, 18, 14, 2, 16, 7, 24, 17, 9, 11, 22, 3, 5, 25, 12, 1, 8, 21)
# mono_famnames <- c("Asteraceae", "Malvaceae", "Asphodelaceae",  "Amaryllidaceae", "Ericaceae", "Fabaceae", "Solanaceae", "Scrophulariaceae", "Campanulaceae", "Cactaceae", "Colchicaceae", "Primulaceae", "Papaveraceae", "Plantaginaceae", "Caprifoliaceae", "Onagraceae", "Euphorbiaceae", "Gesneriaceae", "Iridaceae", "Asparagaceae", "Lamiaceae", "Rosaceae", "Nymphaeaceae", "Paeoniaceae", "Boraginaceae")
# mono_famnodes <- c(156, 204, 237, 227, 175, 187, 152, 145, 164, 179, 242, 172, 216, 144, 167, 211, 201, 146, 240, 235, 138, 195, 243, 214, 147)

# all fams
mono_famnames <- c("Adoxaceae", "Amaryllidaceae", "Apiaceae", "Apocynaceae", "Araliaceae", "Asparagaceae", "Asphodelaceae", "Asteraceae", "Bignoniaceae", "Boraginaceae", "Cactaceae", "Campanulaceae", "Cannaceae", "Caprifoliaceae", "Caryophyllaceae", "Cistaceae", "Colchicaceae", "Convolvulaceae", "Cucurbitaceae", "Droseraceae", "Ericaceae", "Euphorbiaceae", "Fabaceae", "Fagaceae", "Geraniaceae", "Gesneriaceae", "Hydrangeaceae", "Iridaceae", "Lamiaceae", "Linaceae", "Malvaceae", "Menyanthaceae", "Myrtaceae", "Nyctaginaceae", "Nymphaeaceae", "Onagraceae", "Oxalidaceae", "Paeoniaceae", "Papaveraceae", "Phrymaceae", "Pinaceae", "Plantaginaceae", "Polemoniaceae", "Primulaceae", "Ranunculaceae", "Rosaceae", "Salicaceae", "Scrophulariaceae", "Solanaceae", "Strelitziaceae", "Theaceae")
mono_famnums <- c(37, 9, 36, 42, 35, 8, 7, 41, 49, 45, 29, 39, 4, 38, 27, 16, 3, 43, 23, 26, 33, 21, 25, 22, 13, 46, 30, 6, 51, 20, 17, 40, 14, 28, 2, 15, 18, 12, 11, 50, 1, 48, 31, 34, 10, 24, 19, 47, 44, 5, 32)
mono_famnodes <- c(38, 227, 39, 22, 40, 235, 237, 156, 8, 147, 179, 164, 117, 167, 55, 81, 242, 21, 68, 56, 175, 201, 187, 69, 87, 146, 48, 240, 138, 72, 204, 32, 86, 54, 243, 211, 74, 214, 216, 7, 122, 144, 47, 172, 97, 195, 73, 145, 152, 116, 46)


mono_fams <- cbind(mono_famnums, mono_famnames, mono_famnodes)
mono_fams_sorted <- as.data.frame(mono_fams[order(mono_famnums), ])

monocols <- rep(c("gray75", "gray35"), length=length(mono_famnodes))

# scaling pareameter for C1 and C2 dots
k=0.8
# ordering dot data by phylogeny
x <- picante::match.phylo.data(phy, cbind(sterolH, sterolCgrp, sterolDeltagrp, sterolSum))$data
# Colours for sterol heatmap
colfunc <- colorRampPalette(c( "white","red"))
# ordering sterol heatmap data by phylogeny
heat <- picante::match.phylo.data(phy, sterolperc)$data
#check
head(heat)
# Sterol names in desired order
sterolorder <- c("Cycloartenol", "X31Norcycloartanol", "X2425Dehydropollinastanol", "Pollinastanol", "Lathosterol", "Cholesterol", "X31Norcycloartenol", "X14MethylCholest8enol", "Desmosterol", "X24Methylenecholesterol", "X2428MethyleneCycloartanol", "Cycloeucalenol", "Obtusifoliol", "IsoObtusifoliol", "X2428MethyleneLophenol", "Episterol", "Epifungisterol", "Campesterol", "Campestanol", "Citrostradienol", "Schottenol", "Sitosterol", "Sitostanol", "Isofucosterol", "Stigmasterol")
# order heatmap data by sterols
heat <- heat[, sterolorder]
# check
head(heat)
# order abundance and commonness by sterols
richabheat <- richab[, sterolorder]
# check
head(richabheat)
# add column names to be printed
colnames(heat) <- c("1 Cycloartenol", "2 X31Norcycloartanol", "3 X2425Dehydropollinastanol", "4 Pollinastanol", "5 Lathosterol", "6 Cholesterol", "8 X31Norcycloartenol", "9 X14MethylCholest8enol", "10 Desmosterol", "12 X24Methylenecholesterol", "13 X2428MethyleneCycloartanol", "14 Cycloeucalenol", "15 Obtusifoliol", "16 IsoObtusifoliol", "17 X2428MethyleneLophenol", "18 Episterol", "19 Epifungisterol", "20 Campesterol", "21 Campestanol", "22 Avenasterol", "23 Schottenol", "24 Sitosterol", "25 Sitostanol", "26 Isofucosterol", "27 Stigmasterol")
# check
head(heat)
# numbers in heatmap (some skipped)
heatnums <- c(1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27)

#richabheat <- t(picante::match.phylo.data(PathTree, t(richab))$data)
richabheat <- richab[, sterolorder]

# set tip and edge colours of phylogeny (remnant from when branch colours were different)
phylotipcol <- c(rep("black", length(phy$tip.label)))
phyloedgecol <- rep("black", Nedge(phy))
phyloedgecol[which(phy$edge[,2] %in% 1:length(phy$tip.label))] <- phylotipcol  # with this only, it solely colors the

library(scales)
### Full combo figure ##########################################################
brackets = TRUE
################################################################################
pdf(file = "PollenSterolPlot_Final_2020_12_15.pdf", width = 10, height = 14)
plot.new()
#### Tree ######################################################################
par(fig=c(0, 0.35, 0.1, 0.85), new=TRUE)
#par(fig=c(0, 0.45, 0.1, 1), new=TRUE)
  plot(phy, cex=0.5, adj=1, no.margin=TRUE, align.tip.label=TRUE, tip.color=phylotipcol, edge.color=phyloedgecol, edge.width=2)
  obj <- get("last_plot.phylo",envir=.PlotPhyloEnv)
  y.tip <- obj$yy[1:obj$Ntip]
  #y.nodes <- obj$yy[as.numeric(mono_fams_sorted$mono_famnodes)]
  y.nodes <- c(122, 120.5, 118.5, 117, 116, 114.5, 111.5, 108, 102, 97, 93, 88.5, 87, 86, 83.5, 81, 77.5, 74, 73, 72, 70.5, 69, 68, 66, 60.5, 56, 55, 54, 51, 48, 47, 46, 44.5, 42, 40, 39, 38, 36, 33.5, 32, 27, 22, 21, 19.5, 16.5, 13.5, 11.5, 9.5, 8, 7, 3.5)

  for (famis in 1:length(mono_famnodes)) {
    cladebox(phy, as.numeric(mono_fams_sorted$mono_famnodes[famis]), make.transparent(monocols[famis], 0.5))
  }
  if (brackets==TRUE) {
      cladelabels(text=mono_fams_sorted$mono_famnums, node=as.numeric(mono_fams_sorted$mono_famnodes), cex=0.5, offset=20, orientation="horizontal")
  } else {
    text(x=rep(580, length(mono_famnodes)), y=y.nodes, labels=mono_fams_sorted$mono_famnums, srt=0, cex=0.5, adj=0.5, font=2)
  }
  par(fig=c(0, 0.35, 0.1, 0.85), new=TRUE)  # plot tree over the colured boxes again
  plot(phy, cex=0.5, adj=1, no.margin=TRUE, align.tip.label=TRUE, tip.color=phylotipcol, edge.color=phyloedgecol, edge.width=2)

### Main Heatmap ###############################################################
par(fig=c(0.32, 0.66, 0.125, 0.825), new=TRUE)
#par(fig=c(0.4, 0.75, 0.126, 0.9744), new=TRUE)
  image(t(heat), col=colfunc(500), axes=FALSE)
  grid(25, 122, lty=1, lwd=0.5)
### Abundance&Commonness Heatmap ###############################################
par(fig=c(0.32, 0.66, 0.1, 0.115), new=TRUE)
#  plot.new()
  image(t(richabheat), col=colfunc(500), axes=FALSE)
  grid(25, 2, lty=1, lwd=0.5)
### compound numbers between the heatmaps ######################################
par(fig=c(0.31, 0.663, 0.115, 0.126), new=TRUE)
plot.new()
  text(x=seq(from=0, to=1, length.out=length(heatnums)), y=rep(0.45, length(heatnums)), labels=heatnums, srt=0, cex=0.5, adj=0)
### names of Abundance and Commonnness Heatmap #################################
par(fig=c(0.2, 0.4, 0.1, 0.115), new=TRUE)
  plot.new()
  text(x= 0.575, y=0.77, labels="Commonness", col="black", cex=0.55, font=1, adj=1)
  text(x= 0.575, y=0.23, labels="Abundance", col="black", cex=0.55, font=1, adj=1)

#par(fig=c(0.4, 0.1, 0.84, 1), new=TRUE)
### Compound names above heatmap ###############################################
par(fig=c(0.305, 0.91, 0.825, 1), new=TRUE)
  plot.new()
#  plot(PathTree, no.margin=TRUE, cex=0.4, edge.width=1.5, direction="downwards", show.tip.label=TRUE, adj=0.5, align.tip.label=TRUE, label.offset=8, tip.color=sterollabelcol, edge.color=steroledgecol, lwd=2, srt=180)
  text(x=seq(from=0, to=0.58, length.out=length(colnames(heat))), y=rep(0, length(colnames(heat))), labels=colnames(heat), srt=90, cex=0.5, adj=0)
### dot plots for sums and groups ##############################################
par(fig=c(0.59, 0.95, 0.098, 0.85), new=TRUE)
  #dotTree(phy, sterolgrp[, 4:5],length=10, fsize=0.4)
  plot(rep(0, length(y.tip)), y.tip, pch=21, bg="gray", col="black", cex=1.5*k*x[,1], lwd=.4, xlim=c(-0.25, 0.75), xlab=NULL, ylab=NULL, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE, no.margin=TRUE)
  points(rep(0.075, length(y.tip)), y.tip, pch=21, bg=alpha("#808080", 0.65), col="black", cex=k*x[,2], lwd=.4)
  points(rep(0.15, length(y.tip)), y.tip, pch=21, bg=alpha("#808080", 0.75), col="black", cex=k*x[,3], lwd=.4)
  points(rep(0.225, length(y.tip)), y.tip, pch=21, bg=alpha("#808080", 0.90), col="black", cex=k*x[,4], lwd=.4)
  points(rep(0.3, length(y.tip)), y.tip, pch=21, bg="#eeeeee", col="black", cex=k*x[,5], lwd=.4)
  points(rep(0.375, length(y.tip)), y.tip, pch=21, bg="#e8f2a1", col="black", cex=k*x[,6], lwd=.4)
  points(rep(0.45, length(y.tip)), y.tip, pch=21, bg="#afd095", col="black", cex=k*x[,7], lwd=.4)
  points(rep(0.525, length(y.tip)), y.tip, pch=21, bg="#b3cac7", col="black", cex=k*x[,8], lwd=.4)
  points(rep(0.6, length(y.tip)), y.tip, pch=21, bg=alpha("gray", 0.5), col="black", cex=k*log(x[,9]+1), lwd=0.4)
  text(x= 0, y=126, labels="H", col="darkgray", cex=0.55, font=2)
  text(x= 0.075, y=126, labels="C0", col=alpha("#808080", 0.65), cex=0.55, font=2)
  text(x= 0.15, y=126, labels="C1", col=alpha("#808080", 0.75), cex=0.55, font=2)
  text(x= 0.225, y=126, labels="C2", col=alpha("#808080", 0.90), cex=0.55, font=2)
  text(x= 0.3, y=126, labels="D0", col="#eeeeee", cex=0.55, font=2)
  text(x= 0.375, y=126, labels="D5", col="#e8f2a1", cex=0.55, font=2)
  text(x= 0.45, y=126, labels="D7", col="#afd095", cex=0.55, font=2)
  text(x= 0.525, y=126, labels="D8", col="#b3cac7", cex=0.55, font=2)
  text(x=0.6, y=126, labels="Sum", col="darkgray", cex=0.55, font=2)
### legends ####################################################################
### heatmap colour legend ######################################################
par(fig=c(0.3, 0.5, 0.05, 0.08), new=TRUE)
  legend_image <- as.raster((matrix(colfunc(20), nrow=1, ncol=20)))
  plot.new()
  rasterImage(legend_image, 0, 0.5, 1, 1)
  text(x=seq(0.,1,l=5), y = 0.2, labels = seq(0,1,l=5), cex=0.5)
### bubble legend ##############################################################
par(fig=c(0.65, 0.89, 0.0475, 0.0975), new=TRUE)
  plot.new()
  plot(seq(1,4,l=9), rep(0.7, 9), pch = 21, cex=k*c(0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5), col="black", bg=alpha("white", 0.4), lwd=0.4, xlim=c(0,5.2), ylim=c(-0.3,1.2), axes = F, xlab = '', ylab = '')
#  points(1:5, rep(0.35, 5), pch = 21, cex=k*seq(0,1,l=5), col="black", bg="skyblue",lwd=0.4)
  text(x=seq(1,4,l=9), y = -0.1, labels = c(0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5), cex=0.5)

  dev.off()
