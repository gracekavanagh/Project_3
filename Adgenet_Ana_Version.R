#===============================================================================
#                               Welcome!
# This adgenet script was provided by Levi Yant and further altered.
# From an initial VCF file (from GATK) it does a PCA, K-means clustering
# and also calculates a matrix of genetic distances and does an AMOVA.
#
# It was written for plants so it can deal with different ploidies (useful!)
#
# If a comment doesn't have name, it was made by Levi, otherwise it has the
# name of the person that has done it between hashtags
#
# Script by: Filip Kolar 2017, further edits by Sian Bray 2018 & Levi Yant 2022
# altered by: Ana C. da Silva (Matthew Gaskins helped correct for Nas)
# Date: May-July 2022
#===============================================================================

#Ana# install packages if you don't have them!
#install.packages("adegenet", dep=TRUE)
#install.packages("StAMPP")

#Ana# this setting should print warnings as they occur
options(warn=1)

#Ana# call the libraries needed:
library(adegenet)
library(StAMPP)
library(vcfR)
library(ggplot2)
library(MASS)
library(adegraphics) #not strictly necessary for all of this (homebrew R installs will interfere)
library(geosphere) #Ana# added to allow geographic distance in Km
library(dplyr)#Ana# added here to make your life easier
library(tidyverse)
#library(pegas) #Ana# Loaded automatically with StAMPP
#library(ape) #Ana# Loaded automatically with StAMPP
#library(ade4) #Ana# Loaded automatically with adegenet!

################################################################################
######################=========MODIFIED FUNCTIONS=========######################

# a function for conversion from vcfR object to genlight in tetraploids
##Levi##: note not all of this is necessary for LIFE4136 project, but some is helpful
vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0     #Ana# for diploids it's only lines with hashtag below
  x[x == "0|1"] <- 1     #diploid#
  x[x == "1|0"] <- 1     #diploid#
  x[x == "1|1"] <- 2     #diploid#
  x[x == "0/0"] <- 0     #diploid#
  x[x == "0/1"] <- 1     #diploid#
  x[x == "1/0"] <- 1     #diploid#
  x[x == "1/1"] <- 2     #diploid#
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "1/1/1/1/1/1"] <- 6
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}

# a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")

  # keep the original mean / var code, as it's used further down
  # and has some NA checks...
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  
  # convert to full data, try to keep the NA handling as similar to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)

  # all dot products at once using underlying BLAS to support thousands of samples,
  # this could be replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
################################################################################
  ## PERFORM THE ANALYSIS ## ---------------------------------------------------
  
  # eigen analysis

  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  # scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }

  # rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  
  
  ## GET LOADINGS ## -----------------------------------------------------------
  # need to decompose X^TDV into a sum of n matrices of dim p*r
  # but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }

  ## FORMAT OUTPUT ## ----------------------------------------------------------
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}

################################################################################
######################====================================######################
# IMPORT SNP data from VCF
vcfname<-readLines("Vcf_Name.txt")
vcf <- read.vcfR(vcfname)
#vcf <- read.vcfR("Cochlearia_112_dip_tet_4dg.purged.ann.vcf")  

# convert to genlight 	
##Levi## this uses the modified function vcfR2genlight.tetra (see Modified functions section)
aa.genlight <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")  # add real SNP.names

### Levi reverts code back to three letter code in original pipe
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3)               # add pop names: here pop names are first 3 chars of ind name



#check    =====VERY IMPORTANT===
aa.genlight$pop
indNames(aa.genlight)
ploidy(aa.genlight)


#Ana# use this to save file and change names using vcftools so that order keeps the same!
#write.table(indNames(aa.genlight), file="population_names.txt", sep="\t", quote=F, col.names=F)

##############################################################################
# save nice figs

# pdf ("PCA_populations_location.pdf", width=14, height=7)
# ADEgS(c(g1, g2), layout = c(1, 2))
# dev.off()

#ploidy - differentiated plots
# pdf ("Ana_v_PCA_all_ploidycol_SNPs_ax12_1K_less.pdf", width=14, height=7)
# g3 <- s.class(pca.1$scores, as.factor(as.vector(ploidy(aa.genlight))), xax=1, yax=2, col=transp(c("#1e90ff", "#ffa500", "#7cfc00")), 
#               ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plab.cex = 0 , plot = FALSE)
# g4 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
#                                                                                optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
# ADEgS(c(g3, g4), layout = c(1, 2))
# dev.off()

##############################################################################
#===============================================================================
#  distance-based analyses     -------------------------------------------------

# Calculate Nei's distances between individuals/pops
#Ana# Note here raw data is used, not the corrected for NAs
# ---

aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE) # Nei's 1972 distance between indivs
# export matrix - for SplitsTree
stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance_4ds.phy.dst")

aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)   # Nei's 1972 distance between pops
# export matrix - for SplitsTree
stamppPhylip(aa.D.pop, file="aa.pops_Neis_distance_4ds.phy.dst") 

#Ana# the stamppNeisD  is having problem identifying populations here.
# if you use in line 240 pop(aa.genlight) <-substr(indNames(aa.genlight),1,6)
# matrices aa.D.ind and aa.D.pop are the same! what you are calculating is INDIVIDUALS
# since "6" corresponds to full name including the sample nr (1-5)

