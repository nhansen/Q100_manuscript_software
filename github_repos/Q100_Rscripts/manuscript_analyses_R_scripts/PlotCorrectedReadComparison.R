setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Supplement/misc/CorrectionComparison")

library(colorspace)
library(Hmisc)

source("/Users/nhansen/OneDrive/HG002_diploid_benchmark/all_Rscripts/manuscript_analyses_R_scripts/ReadBenchPlotFunctions.R")

outdir <- "dorado_vs_herro"
readsetqvfiles <- c("dorado_vs_herro/herro_corrected.errorqvstats.txt", "dorado_vs_herro/dorado_corrected.errorqvstats.txt")
readsetnames <- c("hg002v1.1_dorado_corrected", "hg002v1.1_herro_corrected")
platformlabels <- c("Dorado corrected", "Herro corrected")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#readplatformcolors <- c("#6699CC", "#CC6677", "#01541F", "#332288", "#661100", "#5D5D5D")
readplatformcolors <- c("#6699CC", "#01541F", "#332288", "#661100", "#5D5D5D")
platformlinetype <- c(1, 2, 3, 4, 5, 6)
platformpchvals <- c(0, 1, 2, 5, 6, 8)
filledplatformpchvals <- c(15, 16, 17, 23, 25, 8)

title <- ""

# There are no read quality scores for the Dorado-corrected and Herro-corrected reads, so skipping QVAccuracy plot for this comparison:
#pdf("Figure4Output/QVAccuracy.pdf")
#read_qv_plot(readsetnames, platformlabels)
#dev.off()

# Plot mononucleotide accuracy

readmononuchist <- function(filename) {
  mndf <- read.table(filename, header=TRUE, sep="")
  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  mndf <- mndf[mndf$readlength != -1, ]
  mndf$reflength <- as.integer(mndf$reflength)
  mndf$readlength <- as.integer(mndf$readlength)
  mndf$numcorrect <- as.integer(mndf$numcorrect)
  mndf$numerror <- as.integer(mndf$numerror)
  
  return(mndf)
}

plotmononucaccuracy <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Accuracy of mononucleotide runs")
  }
  accarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numcorrect"])/sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})
  #  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  if (addtoplot) {
    points(seq(minlength, maxlength), accarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), ylim=c(0,1), col=color)
  }  
  else {
    plot(seq(minlength, maxlength), accarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), main=plottitle, ylim=c(0,1), col=color)
  }
}

plotmononucerrors <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Error rate in mononucleotide runs")
  }
  errorarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"])/sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})
  #  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  if (addtoplot) {
    points(seq(minlength, maxlength), errorarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Error rate"), ylim=c(0,1), col=color)
  }  
  else {
    plot(seq(minlength, maxlength), errorarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Error rate"), main=plottitle, ylim=c(0,1), col=color)
  }
}

plotmononucqvscores <- function(file, plottitle=NA, titlecex=1.0, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE, ymax=35, errorbars=FALSE){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Accuracy of mononucleotide runs")
  }
  errorarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"])})
  totalarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})
  qvarray <- sapply(seq(minlength, maxlength), function(x) {-10.0*log10(sum(mncounts[mncounts$reflength==x, "numerror"])/sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"]))})
  totalerrorrateconfints <- binconf(errorarray, totalarray, return.df=TRUE)
  qvhigharray <- as.numeric(-10.0*log10(totalerrorrateconfints$Lower+0.0000000001))
  qvlowarray <- as.numeric(-10.0*log10(totalerrorrateconfints$Upper))
  if (addtoplot) {
    points(seq(minlength, maxlength), qvarray, pch=pchval, col=color)
    if (errorbars) {
      arrows(x0=seq(minlength, maxlength), y0=qvlowarray, x1=seq(minlength, maxlength), y1=qvhigharray, code=3, angle=90, length=0.05, col=color)
    }
  }  
  else {
    plot(seq(minlength, maxlength), qvarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Phred QV Score"), main=plottitle, cex.main=titlecex, ylim=c(0,ymax), col=color)
    if (errorbars) {
      arrows(x0=seq(minlength, maxlength), y0=qvlowarray, x1=seq(minlength, maxlength), y1=qvhigharray, code=3, angle=90, length=0.05, col=color)
    }
  }
}

plotmononuccoveragecounts <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE, ymax=35){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Number of traversed mononucleotide runs assessed")
  }
  totalarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})/1000
  if (addtoplot) {
    points(seq(minlength, maxlength), totalarray, pch=pchval, col=color, log='y')
  }  
  else {
    plot(seq(minlength, maxlength), totalarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Thousands of reads"), main=plottitle, col=color, log='y', ylim=c(1,1500))
  }
}
read_mononucacc_plot <- function(readsetnames, platformlabels, maxlength=40, platformcolors=readplatformcolors) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".mononuchist.txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucaccuracy(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1])
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucaccuracy(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE)
    }
  }
  legend("bottomleft", platformlabels, col=platformcolors, pch=platformpchvals)
}

read_mononucerror_plot <- function(readsetnames, platformlabels, maxlength=40, platformcolors=readplatformcolors) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".mononuchist.txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucerrors(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1])
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucerrors(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE)
    }
  }
  legend(30, 0.4, platformlabels, col=platformcolors, pch=platformpchvals)
}

read_mononucqvscore_plot <- function(readsetnames, platformlabels, maxlength=40, ymax=45, platformcolors=readplatformcolors, errorbars=FALSE, plottitle=NA, titlecex=1.0, legendcex=1.0) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".mononuchist.txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucqvscores(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1], ymax=ymax, errorbars=errorbars, plottitle=plottitle, titlecex=titlecex)
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucqvscores(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax, errorbars=errorbars)
    }
  }
  legend(32, 42, platformlabels, col=platformcolors, pch=platformpchvals, cex=legendcex)
}

read_mononuccoverage_plot <- function(mnstatsfiles, platformlabels, maxlength=40, ymax=45, platformcolors=readplatformcolors) {
  
  if (length(mnstatsfiles)>0) {
    plotmononuccoveragecounts(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1], ymax=ymax)
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononuccoveragecounts(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax)
    }
  }
  legend(32, 500, platformlabels, col=platformcolors, pch=platformpchvals)
}
### Make mononucleotide accuracy plot ###


pdf("dorado_vs_herro/ReadMononucAccuracy.pdf")  
read_mononucacc_plot(readsetnames, platformlabels)
dev.off()

pdf("dorado_vs_herro/ReadMononucErrorRate.pdf")  
read_mononucerror_plot(readsetnames, platformlabels)
dev.off()

pdf("dorado_vs_herro/ReadMononucQVScores.pdf")  
read_mononucqvscore_plot(readsetnames, platformlabels)
dev.off()

# Plot substitution rate-by-type histogram
pdf("dorado_vs_herro/ReadSubstitutionRates.pdf")
read_substitutions_plot(readsetnames, platformlabels, outputdir="dorado_vs_herro", legend=TRUE, legendposx=4, legendposy=170)
dev.off()

### Plot indel error rate histogram ###
pdf("dorado_vs_herro/ReadIndelRates.pdf")  
read_indels_plot(readsetnames, platformlabels, legendposx=4, legendposy=500)
dev.off()

### Plot the full figure together:
pdf("dorado_vs_herro/ReadQualityMultiplot.pdf", width=11, height=11)
par(mfrow=c(2,2))

read_qv_plot(readsetnames, platformlabels)
read_mononucqvscore_plot(readsetnames, platformlabels, errorbars=TRUE)
read_substitutions_plot(readsetnames, platformlabels, outputdir="Figure4Output")
read_indels_plot(readsetnames, platformlabels)
dev.off()

# compare Illumina old and new
illuminareadsetnames <- c("illumina_2x250", "illumina_googlepcrplus_ds0.1", "illumina_googlepcrfree_ds0.1", "NIST_onso_2024Q1")
illuminaplatformlabels <- c("2x250 (2016)", "Google PCR Plus", "Google PCR Free", "NIST/PB Onso (2024)")

pdf("IlluminaSubsIndelComparison.pdf", width=20, height=11)
par(mfrow=c(1,2))
read_substitutions_plot(illuminareadsetnames, illuminaplatformlabels, outputdir="Figure4Output", legend=TRUE, legendposx=7)
read_indels_plot(illuminareadsetnames, illuminaplatformlabels, legendposx=7, legendposy=25)
dev.off()

pdf("IlluminaMononucAccuracy.pdf", width=11, height=11)
read_mononucqvscore_plot(illuminareadsetnames, illuminaplatformlabels, errorbars=TRUE, plottitle='Illumina Accuracy of Mononucleotide Runs')
dev.off()

# compare Illumina old and new
shortreadsetnames <- c("illumina_2x250", "illumina_googlepcrplus_ds0.1", "illumina_googlepcrfree_ds0.1", "NIST_onso_2024Q1", "element_ultraq_jun2024")
shortplatformlabels <- c("Illumina 2x250 (2016)", "Illumina PCR Plus (2020)", "Illumina PCR Free (2020)", "NIST/PB Onso (2024)", "Element Ultraq (2024)")
shortreadcolors <- c("#88CCEE", "#AA4499", "#661100", "#44AA99", "#332288")


pdf("ShortReadSubsIndelComparison.pdf", width=23, height=11)
par(mfrow=c(1,2))
read_substitutions_plot(shortreadsetnames, shortplatformlabels, platformcolors=shortreadcolors, outputdir="Figure4Output", legend=TRUE, legendcex=1.2, titlecex=1.3, legendposx=10, legendposy=3700)
read_indels_plot(shortreadsetnames, shortplatformlabels, platformcolors=shortreadcolors, legendcex=1.2, titlecex=1.3, legendposx=10, legendposy=28)
dev.off()

pdf("ShortReadMononucAccuracy.pdf", width=11, height=11)
#read_mononucqvscore_plot(shortreadsetnames, shortplatformlabels, errorbars=TRUE, plottitle='Short Read Accuracy of Mononucleotide Runs', titlecex=1.3, legendcex=1.2)
read_mononucqvscore_plot(shortreadsetnames, shortplatformlabels, errorbars=TRUE, plottitle='', legendcex=1.2)
dev.off()



