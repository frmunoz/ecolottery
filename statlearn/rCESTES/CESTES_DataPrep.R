##################################################################### 
############## DATA PREPARATION FOR THE META-STUDY ##################
############## Date: 09/03/2018 | Revision: 26/07/2019
############## Author: Alienor Jeliazkov
##################################################################### 
## Associated publication:
## Jeliazkov A., [78 co-authors] & J. Chase. A global database for metacommunity ecology:
## species, traits, environment and space. Nature's Scientific Data.


############ Load packages #####
library(doBy)
library(readxl)
library(plyr)
library(dplyr)
library(gdata)
library(cellranger)


############ Set working directory and load packages #####
setwd("Path_to_CESTES_data_dir") # directory where ONLY the DATA Excel files are stored
setwd("C:/AJ/POSTDOC_iDIV/PROJECTS/DB_CESTES_REV/xCESTES")
datfiles <- list.files() # list the files in the directory


############ Load the datasets from the directory #####
LSrawdat <- list()
LoadData <- function(datfiles){
  # datfiles = vector of file names of the datasets
  for (i in 1:length(datfiles)) {
    comm <- as.data.frame(read_excel(datfiles[i], sheet = "comm", na="NA"))
    traits <- as.data.frame(read_excel(datfiles[i], sheet = "traits", na="NA"))
    envir <- as.data.frame(read_excel(datfiles[i], sheet = "envir", na="NA"))
    coord <- as.data.frame(read_excel(datfiles[i], sheet="coord", na="NA"))
    if ("blo" %in% excel_sheets(datfiles[i])){
      blo <- as.data.frame(read_excel(datfiles[i], sheet="blo", na="NA"))
      LSrawdat[[i]] <- list(comm=comm, traits=traits, envir=envir, coord=coord, blo=blo)
      names(LSrawdat)[i] <- gsub("_AJ.xlsx", "", datfiles[i])
    } else {
      LSrawdat[[i]] <- list(comm=comm, traits=traits, envir=envir, coord=coord)
      names(LSrawdat)[i] <- gsub("_AJ.xlsx", "", datfiles[i])
    }
  }
  return(LSrawdat)
}
LSrawdat <- LoadData(datfiles)


############ DATA CHECKING 1 - Overall structure #####
### Check NAs (but then, remove manually in the original files to make sensible choices)
### Check empty lines / empty columns in species matrix (and remove manually)
### 'DataCheck' function returns information on where are NAs, ghost species, empty sites, etc.

Single <- function(x) {
  if(length(which(x!=0))<2) 
    return(x)
}

DataCheck <- function(LSrawdat){
  for (i in 1:length(LSrawdat)){
    if(nrow(LSrawdat[[i]]$comm) != nrow(LSrawdat[[i]]$env))
      cat(paste("Warning:", names(LSrawdat[i]), " - comm and env have different number of rows.", sep="\t"), "\n")
    if(nrow(LSrawdat[[i]]$comm) != nrow(LSrawdat[[i]]$coord))
      cat(paste("Warning:", names(LSrawdat[i]), "- comm and coord have different number of rows.", sep="\t"), "\n")
    if(nrow(LSrawdat[[i]]$env) != nrow(LSrawdat[[i]]$coord))
      cat(paste("Warning:", names(LSrawdat[i]), "- env and coord have different number of rows.", sep="\t"), "\n")
    if(ncol(LSrawdat[[i]]$comm[-1]) != nrow(LSrawdat[[i]]$traits))
      cat(paste("Warning:", names(LSrawdat[i]), "- comm and traits have different number of species.", sep="\t"), "\n")
    if(("blo" %in% ls(LSrawdat[[i]])) && (ncol(LSrawdat[[i]]$traits[-1]) != nrow(LSrawdat[[i]]$blo)))
      cat(paste("Warning:", names(LSrawdat[i]), "- traits and blo have different number of lines.", sep="\t"), "\n")
    nam <- names(LSrawdat[[i]])
    for(j in 1:length(nam)){
      df <- LSrawdat[[i]][[nam[j]]][-1]
      if(sum(!complete.cases(df))>0){
        noNA <- sum(!complete.cases(df))
        cat(paste("There are:", noNA, "NAs in the rows of:", names(LSrawdat[i]), nam[j], sep="\t"), "\n")
      }
    }
    if(any(colSums(LSrawdat[[i]]$comm[-1])==0)){
      ghost <- length(which(colSums(LSrawdat[[i]]$comm[-1])==0)==TRUE)
      cat(paste("There are:", ghost, "ghost species in the comm table of:", names(LSrawdat[i]), sep="\t"), "\n")
    }
    # if(any(colSums(LSrawdat[[i]]$comm[-1])==1)){
    #   singlesp <- length(which(colSums(LSrawdat[[i]]$comm[-1])==1))
    #   cat(paste("There are:", singlesp,"singleton species in the comm table of:", names(LSrawdat[i]), sep="\t"), "\n")
    # }
    if(any(rowSums(LSrawdat[[i]]$comm[-1])==0)){
      emptysi <- length(which(rowSums(LSrawdat[[i]]$comm[-1])==0))
      cat(paste("There are", emptysi, "empty sites in the comm table of:", names(LSrawdat[i]), sep="\t"), "\n")
    }
    # if(any(rowSums(LSrawdat[[i]]$comm[-1])==1)){
    #   singlesi <- length(which(rowSums(LSrawdat[[i]]$comm[-1])==1))
    #   cat(paste("There are:", singlesi, "singleton sites in the comm table of:", names(LSrawdat[i]), sep="\t"), "\n")
    # }
    # singleval <- apply(LSrawdat[[i]]$traits[-1], 2, Single)
    # nosingleval <- sapply(singleval, is.null)
    # if(any(nosingleval==FALSE)){
    #   singletra <- sum(nosingleval==FALSE)
    #   cat(paste("There are:", singletra, "singleton values in the trait table of:",  names(LSrawdat[i]), sep="\t"), "\n")
    #   LSrawdat[[i]]$traits <- LSrawdat[[i]]$traits[which(nosingleval)]
    # }
    if(any(!(names(LSrawdat[[i]]$comm[-1]) %in% unique(LSrawdat[[i]]$traits$Sp)))){
      cat(paste("There is a mismatch in species labels between comm and traits.", names(LSrawdat[i]), sep="\t"), "\n")
    }
    if(setequal(LSrawdat[[i]]$comm$Sites, LSrawdat[[i]]$coord$Sites)==FALSE){
      cat(paste("There is a mismatch in site names between comm and coord.", names(LSrawdat[i]), sep="\t"), "\n")
    }
    if(setequal(LSrawdat[[i]]$comm$Sites, LSrawdat[[i]]$envir$Sites)==FALSE){
      cat(paste("There is a mismatch in site names between comm and envir.", names(LSrawdat[i]), sep="\t"), "\n")
    }
  }
}

DataCheck(LSrawdat)


############ DATA CHECKING 2 - Variable types #####
### TRAITS and ENVIR variables are of mixed types:
## Makes sure that: 
# numerics are numerics (and not factors), 
# factors are factors and not characters if stringsAsFactors=TRUE - all characters if FALSE
# ordinal factors are ordinal factors and not numerics,
# integers are integers and not factors,
# binary are numerics (0-1) (although R already converts most of them automatically) and not factors.
setwd("Path_to_Working_dir") # here, use  the path of the WORKING directory
setwd("C:/AJ/POSTDOC_iDIV/PROJECTS/DB_CESTES_REV")
FindBin <- function(x) {if(is.character(x)==TRUE && length(unique(x))==2) return("Bin")}

VarCheck <- function(LSrawdat, str2fac){
  require(dplyr)
  require(doBy)
  capture.output(
    for(j in c("traits", "envir")){
    for(i in 1:length(LSrawdat)) {
      ## binary character to binary 0-1
      binvar <- sapply(LSrawdat[[i]][[j]], FindBin) 
      if(any(binvar=="Bin")){
        nam <- names(which(binvar=="Bin"))
        temp <- as.data.frame(LSrawdat[[i]][[j]][, which(binvar=="Bin")])
        names(temp) <- nam
        LSrawdat[[i]][[j]][,which(binvar=="Bin")] <- as.data.frame(apply(temp, 2, function(x) 
          as.numeric(recodeVar(x, src=unique(x), tgt=c(0,1)))))
      }
      ## character to factor or factors to character
      charvar <- sapply(LSrawdat[[i]][[j]], is.character)
      facvar <- sapply(LSrawdat[[i]][[j]], is.factor)
      if(any(charvar==TRUE) && str2fac==TRUE){
        nam <- names(which(charvar==TRUE))
        temp <- as.data.frame(LSrawdat[[i]][[j]][, which(charvar==TRUE)])
        names(temp) <- nam
        LSrawdat[[i]][[j]][,which(charvar==TRUE)] <- as.data.frame(apply(temp, 2, factor), 
                                                                   stringsAsFactors=TRUE) 
      } else if (any(facvar==TRUE) && str2fac==FALSE) {
        nam <- names(which(facvar==TRUE))
        temp <- as.data.frame(LSrawdat[[i]][[j]][, which(facvar==TRUE)])
        names(temp) <- nam
        LSrawdat[[i]][[j]][,which(facvar==TRUE)] <- as.data.frame(apply(temp, 2, as.character), 
                                                                  stringsAsFactors=FALSE) 
      }
      ## false numeric to ordinal
      h <- names(LSrawdat[[i]][[j]])
      LSrawdat[[i]][[j]][grep("_ORDINAL", h)] <- 
        as.data.frame(sapply(LSrawdat[[i]][[j]][grep("_ORDINAL", h)], ordered))
      ## species are characters, not factors
      LSrawdat[[i]][["traits"]][, "Sp"] <- sapply(LSrawdat[[i]][["traits"]][,"Sp"], as.character)
      ## For a final visual check that numerics are numerics, factors are factors, etc.:
      cat(paste(names(LSrawdat[i]), j, sep=" @@@@@@@ "))
      cat(str(LSrawdat[[i]][[j]]))
    }
  }, file="CheckVarTypes.txt") ###/!\ writes a file into the working directory for visual check /!\
  return(LSrawdat)
} 

LSrawdat <- VarCheck(LSrawdat=LSrawdat, str2fac=TRUE)
#LSrawdat <- VarCheck(LSrawdat=LSrawdat, str2fac=FALSE)


############ DATA CHECKING 3 - Community matrix #####
### Community matrix can include different types of value: 
# counts (integers), densities (decimal), presences/absences (0-1). 
# 'CommCheck' function checks and transforms wherever necessary. 

CommCheck <- function(LSrawdat){
  capture.output(
  for(i in 1:length(LSrawdat)){
    if(length(unique(apply(LSrawdat[[i]]$comm[-1], 2, class)))>1){
      cat(paste("There are heterogeneous types of variables in Comm from:", 
                names(LSrawdat[i]), sep="\t"), "\n")
    }
    cat(paste(names(LSrawdat[i]), "MIN=", min(LSrawdat[[i]]$comm[-1]), "MAX=", 
              max(LSrawdat[[i]]$comm[-1]), "Excerpt:", sep="\t"), "\n")
    if(nrow(LSrawdat[[i]]$comm[-1])>9 & ncol(LSrawdat[[i]]$comm[-1])>9) {
      print(head(LSrawdat[[i]]$comm[-1], 10)[c(1:10)])
    } else {
      print(LSrawdat[[i]]$comm[-1])
    }
  }, file="CheckComm.txt")
}

CommCheck(LSrawdat)
# 12 datasets are PA data: Carvalho2015, Drew2017, Krasnov2015, Purschke2012abcde, Urban2014ab, 
# Westgate2012, Eallonardo2013


############ DATA CHECKING 4 - Spatial coordinates matrix #####
### simply makes sure that coordinates are numerics and that only X and Y are kept in the final list
SpaceCheck <- function(LSrawdat){
  for(i in 1:length(LSrawdat)){
    if(is.numeric(LSrawdat[[i]]$coord$X)==FALSE || is.numeric(LSrawdat[[i]]$coord$Y)==FALSE){
      cat(paste("Coord variables are not numeric in:", names(LSrawdat[i]), sep="\t"), "\n")
      # Convert accordingly
      LSrawdat[[i]]$coord$X <- as.numeric(LSrawdat[[i]]$coord$X)
      LSrawdat[[i]]$coord$Y <- as.numeric(LSrawdat[[i]]$coord$Y)
    }
  }
}
SpaceCheck(LSrawdat)


############ DATA CHECKING 5 - Cleaning process applied and other infos #####
### It can be useful to gather in one file the information of which cleaning procedures
## were applied in each dataset. For this purpose, one can access the 
## "Notes" sheet of each data file. 
## Note that because the "Notes" sheets include heterogeneous information, the information will appear 
## relatively messy in the output file. This is only for quick checking purpose.
setwd("Path_to_CESTES_data_dir") # directory where ONLY the DATA files are stored
LSnotes <- list()
LoadNotes <- function(datfiles){
  # datfiles = vector of file names of the datasets
  for (i in 1:length(datfiles)) {
    if ("Notes" %in% excel_sheets(datfiles[i])){
      LSnotes[[i]] <- read_excel(datfiles[i], sheet = "Notes", na="") # range="Notes!A1:Z200"
      names(LSnotes)[i] <- gsub("_AJ.xlsx", "", datfiles[i])
    } else {
      LSnotes[[i]] <- c("NO NOTES")
      names(LSnotes)[i] <- gsub("_AJ.xlsx", "", datfiles[i])
    }
  }
  return(LSnotes)
}
LSnotes <- LoadNotes(datfiles)

setwd("Path_to_Working_dir") # here, go back to the WORKING directory
NotesCheck <- function(LSnotes){
  capture.output(
    for(i in 1:length(LSnotes)){
      print(LSnotes[[i]])
      cat(paste(names(LSnotes)[i], "#####################################################", sep=""))
    }, file="CheckNotes.txt")
}
NotesCheck(LSnotes)


############ DATA CHECKING 6 - Extraction of DataKey information #####
### Extract the variable names of some traits and environmental variables 
### (as for the Metadata table of the data paper) 
setwd("Path_to_CESTES_data_dir") # directory where ONLY the DATA files are stored
LSkey <- list()
LoadKey <- function(datfiles){
  for (i in 1:length(datfiles)) {
    if ("DataKey" %in% excel_sheets(datfiles[i])){
      LSkey[[i]] <- read_excel(datfiles[i], sheet="DataKey", na="", range=cell_cols("A:C"))
      names(LSkey)[i] <- gsub("_AJ.xlsx", "", datfiles[i])
    } else {
      LSkey[[i]] <- c("NO DATAKEY")
      names(LSkey)[i] <- gsub("_AJ.xlsx", "", datfiles[i])
    }
  }
  return(LSkey)
}
LSkey <- LoadKey(datfiles)

setwd("Path_to_Working_dir") # here, go back to the WORKING directory
KeyCheck <- function(LSkey){
  MetaDataKey <- data.frame(dat=names(LSkey), variables=rep(NA,length(LSkey)))
  capture.output(
    for(i in 1:length(LSkey)){
      toprint <- paste(na.exclude(LSkey[[i]]$Variable), collapse = "|")
      cat(paste(names(LSkey)[i], "#####################################################", sep=""))
      print(toprint)
      MetaDataKey[i,2] <- toprint
    }, file="CheckDataKey.txt")
  return(MetaDataKey)
}
MetaDataKey <- KeyCheck(LSkey)
write.csv(MetaDataKey, file="MetaDataKey.csv") ### writes the information in a csv table


############ DATA PREPARATION FULL VERSION for raw data and trait-based analyses #####
### When "blo" sheet exists, then use the info as weights for the dummy coding
LSmat <- list()

DataPrep <- function(LSrawdat){
  # For some unknown reason, model.matrix does not work within the function's 
  # environment, so one has to run the code from the "for" loop.
  require(vegan)
  require(FD)
  require(stats)
  require(picante)
  require(ape)
  require(ade4)
  require(adespatial)
  nberror = 0
  
  for (i in 1:length(LSrawdat)){
    
    # Transfo community matrix - Hellinger ok for both abund & PA, see Legendre et Gallagher 2001
    speraw <- as.data.frame(LSrawdat[[i]]$comm[-1], row.names=as.integer(LSrawdat[[i]]$comm$Sites))
    matspe <- as.matrix(speraw)
    spehell <- as.matrix(decostand(matspe, "hellinger"))
    
    # Model-matrix of environmental variables
    envir <- as.data.frame(LSrawdat[[i]]$envir[-1], row.names=as.integer(LSrawdat[[i]]$envir$Sites))
    if(ncol(envir) == 1 && is.numeric(envir[,1])==TRUE){
      env <- as.data.frame(model.matrix(~., data=envir)[,-1])
      names(env) <- names(envir[1])
    } else {
      env <- model.matrix(~., data=envir)[,-1] 
    }
      
    # Transfo the environmental data
    ## scale & centre
    matenv <- scale(env, scale=TRUE, center=TRUE) # ok for MRM (Lichstein, 2007) 
    # will return some minor warning messages when the variance is null
    
    # Compute orthonormal envir variables (mix PCA / PCA)
    if(ncol(envir) > 1){
      test <- dudi.hillsmith(envir, scannf=FALSE)
      noaxes <- min(which(inertia(test)$TOT$`cum(%)` > 99))
      dud <- dudi.hillsmith(envir, scannf=FALSE, nf=noaxes)
      envortho <- as.matrix(dud$li)
    } else if (ncol(envir)==1 && is.numeric(envir[,1])==FALSE){
      test <- dudi.pca(env, scannf=FALSE)
      noaxes <- min(which(inertia(test)$TOT$`cum(%)` > 99))
      dud <- dudi.pca(env, scannf=FALSE, nf=noaxes)
      envortho <- as.matrix(dud$li)
    }else{
      envortho <- env
    }

    # Dummy recoding of traits
    tra <- as.data.frame(LSrawdat[[i]]$traits[-1], row.names=LSrawdat[[i]]$traits$Sp)
    if(ncol(tra) == 1 && is.numeric(tra[,1])==TRUE){
      mattra <- as.data.frame(model.matrix(~., data=tra)[,-1])
      names(mattra) <- names(tra[1])
    } else {
      mattra <- model.matrix(~., data=tra)[,-1] # needs to be improved for applying appropriate blo weights
    }
    
    ## Space data
    if((dim(LSrawdat[[i]]$coord)[1]==0)==FALSE){
      spa <- as.data.frame(cbind(X=as.numeric(LSrawdat[[i]]$coord$X), 
              Y=as.numeric(LSrawdat[[i]]$coord$Y)), 
              row.names=as.integer(LSrawdat[[i]]$coord$Sites))
      matspa <- as.matrix(spa)
      matspa <- matspa[order(match(rownames(matspa), rownames(speraw))),,drop=FALSE]
    } else {
      matspa <- as.matrix(NA)
    }
    
    ## Space processing / MEM generation and selection ## /!\ generates errors for several datasets /!\
    # XYpcnm <- pcnm(dist(matspa))
    # pca.species <- dudi.pca(spehell, scale=F, scannf=F)
    # rda.sp <- pcaiv(pca.species, XYpcnm$vectors, scannf=F)
    # R2 <- sum(rda.sp$eig)/sum(pca.species$eig)
    # # matmem <- as.matrix(XYpcnm$vectors[,which(XYpcnm$values>0)]) # keep only MEMs of positive value
    # fw <- try(forward.sel(matspe, XYpcnm$vectors, alpha = 0.05,
    #     adjR2thresh=RsquareAdj(R2, nrow(XYpcnm$vectors), ncol(XYpcnm$vectors)))) 
    # if(class(fw) != "try-error"){ 
    #   fw <- forward.sel(matspe, XYpcnm$vectors, alpha = 0.05,
    #      adjR2thresh=RsquareAdj(R2, nrow(XYpcnm$vectors), ncol(XYpcnm$vectors)))
    #   matmem <- as.matrix(XYpcnm$vectors[,fw$order])
    #   colnames(matmem) <- fw$variables
    # } else {
    #   nberror <- nberror+1
    #   matmem <- matspa
    # }
    
    ## Check and match order of sites and species across matrices
    speraw <- speraw[order(as.numeric(rownames(speraw))),,drop=FALSE]
    matspe <- matspe[order(match(rownames(matspe), rownames(speraw))),,drop=FALSE]
    tra <- tra[order(match(rownames(tra), colnames(matspe))),,drop=FALSE]
    spehell <- spehell[order(match(rownames(spehell), rownames(matspe))),,drop=FALSE]
    matenv <- matenv[order(match(rownames(matenv), rownames(matspe))),,drop=FALSE]
    envortho <- envortho[order(match(rownames(envortho), rownames(matspe))),,drop=FALSE]
    env <- env[order(match(rownames(env), rownames(matspe))),,drop=FALSE]
    #matmem <- matmem[order(match(rownames(matmem), rownames(matspe))),,drop=FALSE]
    
    ## Compute Gower distances (accounting for potential blocks) # USE SYNTHETIC or RAW TRAITS 
    if("blo" %in% names(LSrawdat[[i]])) {
      w <- unlist(with(LSrawdat[[i]]$blo, by(LSrawdat[[i]]$blo, as.factor(Blocks), function(x) rep(1/nrow(x),nrow(x)))))
      gowdist <- gowdis(tra, w=w) # here the raw traits are used, acc. to Dray's & ter Braak's works (dbCCA)
      # blo <- LSrawdat[[i]]$blo
    } else {
      gowdist <- gowdis(tra)
      # blo <- NULL
    }
    
    ## Compute FDis (Functional Dispersion, Laliberte & Legendre 2010) close to Rao's QE
    # fundisp <- data.frame(FDis=fdisp(gowdist, matspe)$FDis)
    # w <- prep.fuzzy.var(as.data.frame(mattra), col.blocks=as.vector(blo))
    # pcaThv <- dudi.fca(w, scannf=FALSE)
 
    ## Compute Community Weighted Means (Laliberte's method)
    # First log transfo and scale traits: apply(log(traits), 2, scale) 
    # Then apply a PCA / PCoA / dimension reduction for correcting trait covariation 
    # Use the synthetic axes as traits to compute the CWM
    ## With raw traits
    CMW <- functcomp(tra, matspe, CWM.type="all")
    ## With orthonormalized traits
    if(ncol(tra) > 1){
      test <- dudi.hillsmith(tra, scannf=FALSE)
      noaxes <- min(which(inertia(test)$TOT$`cum(%)` > 95))
      mixpca <- dudi.hillsmith(tra, scannf=FALSE, nf=noaxes)
      CMWortho <- functcomp(mixpca$li, matspe, CWM.type="all") # orthonormalized traits via mix PCA
      # pcoaT <- ape::pcoa(gowdist, correction="lingoes")
      # CMWortho <- functcomp(pcoaT$vectors, matspe, CWM.type="all") # orthonormalized traits via PCoA
    } else {
      CMWortho <- mattra
    }
    
    ## Compute dendrogram-based FD (e.g. Podani & Schmera 2006) 
    # dendro <- hclust(gowdist, method="average")
    # coph <- cophenetic(dendro)
    # MPD <- taxondive(comm=speraw, dis=gowdist) # or dis=coph or dis=gowdist
    # mat <- mpd(samp=LSmat[[1]]$speraw, dis=as.matrix(LSmat[[1]]$gowdist), abundance.weighted=TRUE) # Rao's index
    # matMPD <- dist(summary(MPD)[-nrow(summary(MPD)),"Delta"])
    # cor.test(bcdist(spehell), matMPD)
    
    LSmat[[i]] <- list(spe=spehell, speraw=matspe, cmw=CMW, tra=tra, env=matenv, 
                       envraw=envir, envortho=envortho, spa=matspa, dumtra=mattra,  
                       gowdist=gowdist, dumenv=env, cmwortho=CMWortho) 
    # to uncomment and restore in the list above if needed: MPD=MPD, fundisp=fundisp, mems=matmem, blo=blo
    names(LSmat)[i] <- names(LSrawdat)[i]
  }
  return(LSmat)
  cat(paste("nberror = ", nberror, sep="")) 
}

LSmat <- DataPrep(LSrawdat)


### Quick check of the gowdist matrices 
for(i in 1:length(LSmat)){
  cat(head(LSmat[[i]]$gowdist), names(LSmat)[i], sep="\n")
}

### Save the list
save(LSmat, file="LSmat80.RData") # save the list in an R object


############ DATA PREPARATION MINIMAL VERSION for raw data only #####
LSmatred <- lapply(LSmat, function(x) list(comm=x$speraw, traits=x$tra, envir=x$envraw, coord=x$spa))
save(LSmatred, file="CESTES.RData") # save the list in an R object


################## FINAL OUTPUT OF THE DATA PREPARATION FULL VERSION ##########
# a list with all the datasets: LSdat
# every element should include these matrices:
# 1. cmw = Community Weighted Means calculated from {FD}
# 2. cmwortho = Community Weighted Means calculated from {FD} based on orthonormalised traits (mix PCA)
# 3. dumenv = dummy recoded environmental variables (i.e. factor levels are splitted into binary variables)
# 4. dumtra = dummy recoded traits (i.e. factor levels are splitted into binary variables)
# 5. env = scaled environmental variables from dumenv (scaled and centered) 
# 6. envortho = orthonormalised environmental components
# 7. envraw = raw environmental variables
# 8. gowdist = Gower distance matrix computed on raw traits
# 9. spa = matrix of spatial coordinates
# 10. spe =  spehell = hellinger transformed species matrix
# 11. speraw = species original abundances/PA/densities...
# 12. tra = raw trait data
# (opt.) mems (spatial MEMs before or after forward selection)
# (opt.) fundisp = functional dispersion calculated from {FD}
# (opt.) blo = blocks of traits when relevant (otherwise, blo=NULL)


################## FINAL OUTPUT OF THE DATA PREPARATION MINIMAL VERSION ##########
# a list with all the datasets: LSdat
# every element should include these matrices:
# 1. comm = raw community matrix 
# 2. traits = raw trait matrix
# 3. envir = raw environment matrix
# 4. coord = raw spatial coordinates matrix


#### Other useful functions to load ####
# From: https://gist.github.com/pimentel/256fc8c9b5191da63819#file-head-list-r
head.list <- function(obj, n = 6L, ...)
{
  stopifnot(length(n) == 1L)
  origN <- n
  n <- if (n < 0L)
    max(length(obj) + n, 0L)
  else min(n, length(obj))
  lapply(obj[seq_len(n)], function(x)
  {
    tryCatch({
      head(x, origN, ...)
    }, error = function(e) {
      x
    })
  })
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y)) # abs() could be removed if we want the sign of the corr
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
} 

