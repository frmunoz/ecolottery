################################################################## 
############## PLOTS OF THE DATA PAPER & EXPLORATION #############
############## Date: 25/03/2019 | Revision: 26/07/2019
############## Author: Alienor Jeliazkov
##################################################################
## Associated publication:
## Jeliazkov A., [78 co-authors] & J. Chase. A global database for metacommunity ecology:
## species, traits, environment and space. Nature's Scientific Data.

library(ggplot2)
library(gridExtra)
library(MASS)
library(plyr)
library(plotrix)
library(devtools)
library(easyGgplot2)
library(scales)
library(treemap)
library(cowplot)
library(doBy)

####### Set directory and load the meta-data ####### 
setwd("rCESTES") # path to rCESTES folder 
Metadat <- read.csv("Metadat.csv")
Listdat <- read.csv("ListDat.csv")
load("CESTES.RData")
names(LSmatraw)

####### Pie charts of the meta-data #######
explo <- table(Metadat$Group)
slices <- explo
lbls <- names(explo)
pie3D(slices,labels=lbls, explode=0, labelcex=1, mar=c(0,0,0,0), 
      border="white", shade=0.6)
explo <- table(Metadat$Ecosystem_linked)
slices <- explo
lbls <- names(explo)
pie3D(slices,labels=lbls, explode=0, labelcex=1, mar=c(0,0,0,0), 
      border="white", shade=0.6, col=c("lightblue", "darkblue", "brown"))
explo <- table(Metadat$Hemeroby_4)
slices <- explo
lbls <- names(explo)
pie3D(slices,labels=lbls, explode=0, labelcex=1, mar=c(0,0,0,0), 
      border="white", shade=0.6, col=c("grey", "#E6E600FF", "darkgreen", "lightgreen"))


####### Barplots of the meta-data #######
df1 <- as.data.frame(table(Metadat$Group))
df2 <- as.data.frame(table(Metadat$Ecosystem_linked))
df3 <- as.data.frame(table(Metadat$Hemeroby_4))
df1s <- df1[order(df1$Freq, decreasing=TRUE),]
df1s$Var1 <- factor(df1s$Var1, 
                    levels=df1s$Var1[order(df1s$Freq, decreasing=TRUE)])
df2s <- df2[order(df2$Freq, decreasing=TRUE),]
df2s$Var1 <- factor(df2s$Var1, 
                    levels=df2s$Var1[order(df2s$Freq, decreasing=TRUE)])
df3s <- df3[order(df3$Freq, decreasing=TRUE),]
df3s$Var1 <- factor(df3s$Var1, 
                    levels=df3s$Var1[order(df3s$Freq, decreasing=TRUE)])

piebar1 <- ggplot(df1s, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") + theme_minimal() +  
  xlab("Taxonomic group") + ylab("Number of datasets") + 
  theme(legend.position = "none", axis.text.y=element_text(size=17), 
        axis.text.x=element_text(size=17, angle=45, hjust=1), 
        axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20))
piebar2 <- ggplot(df2s, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") + theme_minimal() +  
  xlab("Ecosystem") + ylab("") + 
  theme(legend.position = "none", axis.text.y=element_text(size=17), 
        axis.text.x=element_text(size=17, angle=45, hjust=1), 
        axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20)) + 
  scale_fill_manual(values=c("brown4", "skyblue", "blue"))
piebar3 <- ggplot(df3s, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity") + theme_minimal() +  
  xlab("Level of human disturbance") + ylab("Number of datasets") + 
  theme(legend.position = "none", axis.text.y=element_text(size=17), 
        axis.text.x=element_text(size=17, angle=45, hjust=1), 
        axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20)) + 
  scale_fill_manual(values=c("black", "grey30", "grey60", "grey"))

### Histogram Extent
histo <- ggplot(data = Metadat, aes(x=Extent_km2)) + 
  geom_histogram(bins=10, fill="springgreen4") +
  scale_x_log10(labels = comma) + xlab("Extent (km2)") + 
  theme(axis.text.y=element_text(size=17), axis.text.x=element_text(size=15, angle=30, hjust=1), 
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), 
        plot.margin = margin(0.5, 0.8, 0, 0, "cm")) + ylab("")

grid.arrange(piebar1, piebar2, piebar3, histo, ncol=2)


####### Plots of the sampling properties #######
## Data reshaping
nbETS <- reshape(Metadat[c("dat", "nbEnv","nbTra","nbSpe", "nbSit")], idvar="data", 
                  varying=list(2:5), times=names(Metadat[(8:11)]), v.names="nb", direction="long")
nbETS$type <- recodeVar(nbETS$time, src=unique(nbETS$time), 
                         tgt=c("Envir & traits", "Envir & traits", "Sp & sites", "Sp & sites"))
nbETS$leg <- recodeVar(nbETS$time, src=unique(nbETS$time), 
                        tgt=c("Environment", "Traits", "Species", "Sites"))

#### Violin plots (not in the paper)
### Plots of mean/variation nb of sites, nb of species and nb of traits per dataset
ggplot2.violinplot(data=nbETS, xName='time', yName='nb', 
                   groupName='leg', groupColors=c('aquamarine3','chartreuse1','goldenrod1', 'red'), 
                   addDot=TRUE, dotSize=0.5, faceting=TRUE, facetingVarNames="type",
                   facetingDirection="vertical", facetingScales="free", 
                   addMean=TRUE, meanPointShape=21, meanPointSize=4,
                   meanPointColor="black", meanPointFill="black", 
                   xtitle="Data component", xShowTitle=FALSE, ytitle="Number of variables per dataset",
                   legendTitle="Data component", xShowTickLabel=FALSE) # yScale='log10'

#### Histograms 
gsamp <- ggplot(data=nbETS, aes(nb)) +
  geom_histogram(bins=19, aes(fill=leg)) + facet_wrap(leg ~ ., scales = "free") + 
  theme(legend.position = "none", axis.text.y=element_text(size=17), 
        axis.text.x=element_text(size=17, angle=45, hjust=1), 
        axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), strip.text = element_text(size=20), 
        panel.grid.major = element_line("lightgray",0.5)) + 
  xlab("Number of variables") + ylab("Number of datasets")
gsamp


####### Plots of the data request success #######
explo <- table(Listdat$RequestSuccess)
explo2 <- as.data.frame(explo)
### Tree map 
piepercent <- round(100*explo2$Freq/sum(explo2$Freq), 1)
explo2$Var2 <- paste(explo2$Var1, paste(piepercent, "%", sep=""), sep='\n')
treemap(explo2, index="Var2", vSize="Freq", type="index", title="", palette="Set2")
explo2$perc <- round(100*explo2$Freq/sum(explo2$Freq), 1)
### Barplot (after revision 26/07/2019)
explo2o <- explo2[order(explo2$perc, decreasing=TRUE),]
explo2o$Var1 <- factor(explo2o$Var1, 
                       levels=explo2o$Var1[order(explo2o$perc, decreasing=TRUE)])
req <- ggplot(explo2o, aes(x=Var1, y=perc, fill=Var1)) +
  geom_bar(stat="identity") + theme_minimal() +  
  xlab("Outputs from the data collection") + ylab("% of eligible datasets") + 
  theme(legend.position = "none", axis.text.y=element_text(size=17), 
        axis.text.x=element_text(size=17, angle=45, hjust=1), 
        axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20))


###### Power zone plot based on Dray et Legendre 2008  ###### 
plot(Metadat$nbSit ~ Metadat$nbSpe, pch=16, xlab="Number of species", cex=1.5,
     ylab="Number of sites", cex.axis=1.5, cex.lab=1.5)
points(Metadat$nbSpe, Metadat$nbSit, pch=16, cex=0.7,col="lightblue")
rect(0, 10, 100, 100, lty=3)
spref.95 <- c(18,20,30,40,50,60,70,80,90,100)
traref.95 <- c(100,70,45,35,28,25,20,19,18,17)
spref.90 <- c(15,20,30,40,50,60,70,80,90,100)
traref.90 <- c(100,60,35,29,22,20,18,17,15,14)
spref.70 <- c(10,20,30,40,50,60,70,80,90,100)
traref.70 <- c(60,26,19,15,12, NA, NA, NA, NA, NA)
spref.100 <- c(36,40,54,50,49,50,60,70,75,75, 80,90,100)
traref.100 <- c(100, 90,70, 65,60,59,58,58,50,40,39,38, 45)
lines(traref.100 ~ spref.100, col="grey", lwd=3, lty=1)
lines(traref.95 ~ spref.95, col="yellow", lwd=3, lty=1)
lines(spref.90, traref.90, col="orange", lwd=3, lty=1)
lines(spref.70, traref.70, col="red", lwd=3, lty=1)
