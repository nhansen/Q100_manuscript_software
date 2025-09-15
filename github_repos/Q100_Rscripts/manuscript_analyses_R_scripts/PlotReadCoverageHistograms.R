setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Supplement/coveragetests")

library(colorspace)
library(Hmisc)
library(plotly)
library(tidyverse)
library(stats)
library(ggplot2)
library(RColorBrewer)

################
### Supplementary Figure? ###
### Sequencing coverage of reads divided into bins so that we can assess
### deviations from Poisson behavior
################

outdir <- "."

# currently plotted read sets (alternatives below)
readsetnames <- c("ont_epi2me_q28", "hifi_revio_pbmay24", "element_ultraq_jun2024", "illumina_googlepcrfree")
platformlabels <- c("ONT Q28", "HiFi Revio", "Element", "Illumina")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
readplatformcolors <- c("#6699CC", "#01541F", "#332288", "#661100", "#5D5D5D")
platformlinetype <- c(1, 2, 3, 4, 5, 6)
platformpchvals <- c(0, 1, 2, 5, 6, 8)
filledplatformpchvals <- c(15, 16, 17, 18, 23, 25, 8)
multiplevector <- c(1.0, 1.0, 1.0, 1.4, 1.0, 1.0, 1.0, 1.4, 1.4, 1.4, 1.4, 1.9)

title <- ""

hificov_filename <- "hifi_revio_pbmay24.binnedcoverage.bed"
hificovdf <- read.table(hificov_filename, sep="\t")
names(hificovdf) <- c("chrom", "start", "end", "arrivalcount")
hifimean <- mean(hificovdf$arrivalcount)
hifivar <- var(hificovdf$arrivalcount)

hifikmer_filename <- "hifi_revio_pbmay24.extremekmercounts.bed"
hifikmerdf <- read.table(hifikmer_filename, sep="\t")
names(hifikmerdf) <- c("chrom", "start", "end", "countat", "countag", "countac", "countgc")
#hifimean <- mean(hificovdf$arrivalcount)
#hifivar <- var(hificovdf$arrivalcount)

ontcov_filename <- "ont_epi2me_q28.binnedcoverage.bed"
ontcovdf <- read.table(ontcov_filename, sep="\t")
names(ontcovdf) <- c("chrom", "start", "end", "arrivalcount")
ontmean <- mean(ontcovdf$arrivalcount)
ontvar <- var(ontcovdf$arrivalcount)

ontkmer_filename <- "ont_epi2me_q28.extremekmercounts.bed"
ontkmerdf <- read.table(ontkmer_filename, sep="\t")
names(ontkmerdf) <- c("chrom", "start", "end", "countat", "countag", "countac", "countgc")

illuminacov_filename <- "illumina_googlepcrfree.binnedcoverage.bed"
illuminacovdf <- read.table(illuminacov_filename, sep="\t")
names(illuminacovdf) <- c("chrom", "start", "end", "arrivalcount")
illuminamean <- mean(illuminacovdf$arrivalcount)
illuminavar <- var(illuminacovdf$arrivalcount)

illuminakmer_filename <- "illumina_googlepcrfree.extremekmercounts.bed"
illuminakmerdf <- read.table(illuminakmer_filename, sep="\t")
names(illuminakmerdf) <- c("chrom", "start", "end", "countat", "countag", "countac", "countgc")

elementcov_filename <- "element_ultraq_jun2024.binnedcoverage.bed"
elementcovdf <- read.table(elementcov_filename, sep="\t")
names(elementcovdf) <- c("chrom", "start", "end", "arrivalcount")
elementmean <- mean(elementcovdf$arrivalcount)
elementvar <- var(elementcovdf$arrivalcount)

elementkmer_filename <- "element_ultraq_jun2024.extremekmercounts.bed"
elementkmerdf <- read.table(elementkmer_filename, sep="\t")
names(elementkmerdf) <- c("chrom", "start", "end", "countat", "countag", "countac", "countgc")

plotcomparisonhist <- function(covdf, kmerdf, key="countat", ymax=0.01, lowcovlim=NA, xmax=500) {
  allvals <- kmerdf[,key]
  if (is.na(lowcovlim)) {
    meancov <- mean(covdf$arrivalcount)
    sdcov <- sd(covdf$arrivalcount)
    lowcovlim <- meancov - 3*sdcov
  }
  lowcovvals <- kmerdf[covdf$arrivalcount < lowcovlim, key]
  hist(allvals, prob=TRUE, ylim=c(0, ymax), xlim=c(0,xmax))
  lines(density(allvals), col=4, lwd=2)
  lines(density(lowcovvals), col=6, lwd=2)
}

plotcomparisonhist(hificovdf, hifikmerdf, xmax=100, key="countag", ymax=0.03)

# This function was moved to the file "ReadBenchComparisonPlotFunctions.R"
#plotcdfcurves(ontcovdf$arrivalcount, xmax=500, datacolor='darkgreen', title='Cumulative distribution of ONT read coverage')
#plotcdfcurves(hificovdf$arrivalcount, addedcurve=TRUE, xmax=500, datacolor='darkred', title='Cumulative distribution of HiFi Revio read coverage')
#plotcdfcurves(illuminacovdf$arrivalcount, addedcurve=TRUE, xmax=500, datacolor='lightblue', title='Cumulative distribution of Illumina read coverage')
#plotcdfcurves(elementcovdf$arrivalcount, addedcurve=TRUE, xmax=500, datacolor='pink', title='Cumulative distribution of Element read coverage')


plot2basecontentvscoverage <- function(covdf, kmerdf, key="countag", maxcount=3000, minx=NA, maxx=NA, miny=NA, maxy=NA, title="", hex=FALSE) {
  binsize <- covdf[1,"end"]-covdf[1,"start"]
  plotdf <- data.frame("arrivalcount"=covdf[covdf$arrivalcount < maxcount, "arrivalcount"], "twobasetype"=as.double(kmerdf[covdf$arrivalcount < maxcount,key]/binsize))
  d <- ggplot(data=plotdf, aes(x=arrivalcount, y=twobasetype))
  #if (is.na(minx)) {
    #minx <- min(plotdf$arrivalcount)
  #}
  #if (is.na(maxx)) {
    #maxx <- max(plotdf$arrivalcount)
  #}
  #if (is.na(miny)) {
    #miny <- min(plotdf$twobasetype)
  #}
  #if (is.na(maxy)) {
    #maxy <- max(plotdf$twobasetype)
  #}
  twobase <- toupper(substr(key, 6, 7))
  #d + geom_density2d()
  if (binsize > 1000000) {
    xlabel <- paste0("Reads per ", as.character(binsize/1000000), " Mb window")
    ylabel <- paste0("Fraction of window that is ", twobase, " kmers")
  }
  else if (binsize > 1000) {
    xlabel <- paste0("Reads per ", as.character(binsize/1000), " kb window")
    ylabel <- paste0("Fraction of window that are ", twobase, " kmers")
  }
  else {
    xlabel <- paste0("Reads per ", as.character(binsize), " base window")
    ylabel <- paste0("Fraction of window that are ", twobase, " kmers")
  }
  if(hex) {
    d + geom_hex() + ggtitle(title) + xlab(xlabel) + ylab(ylabel)
    #+ coord_cartesian(xlim = c(minx, maxx), ylim=c(miny, maxy))
  }
  else {
    d + geom_density2d() + ggtitle(title) + xlab(xlabel) + ylab(ylabel)
    #+ coord_cartesian(xlim = c(minx, maxx), ylim=c(miny, maxy))
  }
}

# ont:
plot2basecontentvscoverage(ontcovdf, ontkmerdf, "countag", title="ONT AG Content vs Coverage", color="blue")
plot2basecontentvscoverage(ontcovdf, ontkmerdf, "countac", title="ONT AC Content vs Coverage")
plot2basecontentvscoverage(ontcovdf, ontkmerdf, "countat", title="ONT AT Content vs Coverage")

# hifi:
plot2basecontentvscoverage(hificovdf, hifikmerdf, "countag", title="HiFi AG Content vs Coverage", color="red")
plot2basecontentvscoverage(hificovdf, hifikmerdf, "countac", title="HiFi AC Content vs Coverage")
plot2basecontentvscoverage(hificovdf, hifikmerdf, "countat", title="HiFi AT Content vs Coverage")

# illumina:
plot2basecontentvscoverage(illuminacovdf, illuminakmerdf, "countag", title="Illumina AG Content vs Coverage")
plot2basecontentvscoverage(illuminacovdf, illuminakmerdf, "countac", title="Illumina AC Content vs Coverage")
plot2basecontentvscoverage(illuminacovdf, illuminakmerdf, "countat", title="Illumina AT Content vs Coverage")

# element:
plot2basecontentvscoverage(elementcovdf, elementkmerdf, "countag", title="Element AG Content vs Coverage")
plot2basecontentvscoverage(elementcovdf, elementkmerdf, "countac", title="Element AC Content vs Coverage")
plot2basecontentvscoverage(elementcovdf, elementkmerdf, "countat", title="Element AT Content vs Coverage")

plot(hificovdf$arrivalcount, hifikmerdf$countag, xlim=c(0,1500), xlab="Arrival count per bin", ylab="Number of AG kmers in bin", main="AG kmer content vs. read coverage for HiFi reads")
plot(hificovdf$arrivalcount, hifikmerdf$countat, xlim=c(0,1500), xlab="Arrival count per bin", ylab="Number of AT kmers in bin", main="AT kmer content vs. read coverage for HiFi reads")
plot(hificovdf$arrivalcount, hifikmerdf$countac, xlim=c(0,1500), xlab="Arrival count per bin", ylab="Number of AC kmers in bin", main="AC kmer content vs. read coverage for HiFi reads")

plot(elementcovdf$arrivalcount, elementkmerdf$countag, xlim=c(0,2000), xlab="Arrival count per bin", ylab="Number of AG kmers in bin", main="AG kmer content vs. read coverage for Element reads")
plot(elementcovdf$arrivalcount, elementkmerdf$countac, xlim=c(0,2000), xlab="Arrival count per bin", ylab="Number of AC kmers in bin", main="AC kmer content vs. read coverage for Element reads")

ks.test()
