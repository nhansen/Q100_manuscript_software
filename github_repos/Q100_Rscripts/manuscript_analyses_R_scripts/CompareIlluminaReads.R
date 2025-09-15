setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/PaperFigures/Figures")

library(colorspace)
library(Hmisc)

source("/Users/nhansen/OneDrive/HG002_diploid_benchmark/all_Rscripts/manuscript_analyses_R_scripts/ReadBenchComparisonPlotFunctions.R")

################
### Comparison of new "advanced HP" Illumina reads with various old ones ###
################

outdir <- "Figure4Output"
#readsetnames <- c("illumina_advancedhp_rep1", "illumina_advancedhp_rep2", "illumina_advancedhp_rep3", "illumina_2x250", "illumina_googlepcrfree", "illumina_googlepcrplus")
#platformlabels <- c("AdvancedHP Rep1", "AdvancedHP Rep2", "AdvancedHP Rep3", "GIAB 2x250", "Google PCR free", "Google PCR plus")
readsetnames <- c("illumina_advancedhp_rep1", "illumina_advancedhp_rep2", "illumina_advancedhp_rep3", "illumina_2x250", "illumina_googlepcrfree", "element_ultraq_jun2024")
platformlabels <- c("AdvancedHP Rep1", "AdvancedHP Rep2", "AdvancedHP Rep3", "GIAB 2x250", "Google PCR free", "Element Ultraq Jun2024")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
readplatformcolors <- c("#01541F", "#332288", "#661100", "#5D5D5D", "#CC6677", "#999933")
platformlinetype <- c(1, 2, 3, 4, 5, 6)
platformpchvals <- c(0, 1, 2, 5, 6, 8)
filledplatformpchvals <- c(15, 16, 17, 23, 25, 8)

title <- ""

#pdf("Figure4Output/IlluminaQVAccuracy.pdf", width=11, height=11)
read_qv_plot(readsetnames, platformlabels)
#dev.off()


### Make mononucleotide accuracy plot ###


#pdf("Figure4Output/IlluminaReadMononucAccuracy.pdf")
pdf("Figure4Output/IlluminaReadMononucAccuracyWithElement.pdf")
#png("Figure4Output/IlluminaReadMononucAccuracy.png")
read_mononucacc_plot(readsetnames, platformlabels)
dev.off()

#pdf("Figure4Output/IlluminaReadMononucErrorRate.pdf")  
read_mononucerror_plot(readsetnames, platformlabels)
#dev.off()

#pdf("Figure4Output/IlluminaQVScores.pdf")  
read_mononucqvscore_plot(readsetnames, platformlabels)
#dev.off()

#pdf("Figure4Output/IlluminaReadSubstitutionRates.pdf", width=18, height=11)
pdf("Figure4Output/IlluminaReadSubstitutionRates.pdf", width=18, height=8)
read_substitutions_plot(readsetnames, platformlabels, outputdir="Figure4Output", legend=TRUE, ymax=4200, legendposx=4, legendposy=4200)
dev.off()

### Indel error plot ###

totalindelrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist$indelspermbaligned, na.rm=TRUE))  
}

totalinsertionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist[indellengthhist$indellength>0, "indelspermbaligned"], na.rm=TRUE))  
}

totaldeletionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist[indellengthhist$indellength<0, "indelspermbaligned"], na.rm=TRUE))  
}

read_indels_plot <- function(indelreadsetnames, platformlabels, platformcolors=readplatformcolors, xlabval="Read platform", ylabval="Indels per mb", titleval="Indel errors in reads", titlecex=1.0, ymax=NA, legendcex=1.0, legendposx=NA, legendposy=NA) {
  indellengthfiles <- sapply(indelreadsetnames, function(x) {paste(c(outdir, "/", x, ".indelerrorstats.txt"), sep="", collapse="")})
  totalinsertionrates <- sapply(indellengthfiles, totalinsertionrate)
  totaldeletionrates <- sapply(indellengthfiles, totaldeletionrate)
  
  platformlabelswithinsdel <- sapply(platformlabels, function(x) {label = paste(c(x, "\nIns/Del"), sep="", collapse=""); return(label)})
  
  barcolors <- sapply(platformcolors, function(x) {c(darken(x), lighten(x, 0.3))})
  out <- barplot(rbind(totalinsertionrates, totaldeletionrates), beside=TRUE, names.arg=platformlabelswithinsdel, col=barcolors, main=titleval, cex.main=titlecex, xlab=xlabval, ylab=ylabval)
  if (is.na(legendposx)) {
    legendposx <- out[2*length(totalinsertionrates)-2]
  }
  if (is.na(legendposy)) {
    legendposy <- max(rbind(totalinsertionrates, totaldeletionrates))-100
  }
  legend(legendposx, legendposy, platformlabels, col=platformcolors, pch=15, cex=legendcex)
  formattedrates <- as.integer((rbind(totalinsertionrates, totaldeletionrates)+0.5)*10)/10
  formattedrates <- ifelse(formattedrates>100, as.integer(formattedrates), formattedrates)
  text(out, rbind(totalinsertionrates, totaldeletionrates), formattedrates, pos=3, xpd=NA, cex=0.85)
}

pdf("Figure4Output/IlluminaReadIndelRates.pdf", width=18, height=8)  
read_indels_plot(readsetnames, platformlabels, legend=TRUE, legendposx=4, legendposy=25)
dev.off()

### Plot the full figure together:
pdf("PacBioHiFiComparisonMultiplot.pdf", width=13, height=11)
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



