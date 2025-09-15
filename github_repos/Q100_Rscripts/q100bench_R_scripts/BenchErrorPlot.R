setwd("/Users/nhansen/HG002_diploid_benchmark/plots/hetsite_plots")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#if (!requireNamespace("karyoploteR", quietly = TRUE))
  #BiocManager::install("karyoploteR")

library(stringr)
library(karyoploteR)
args = commandArgs(trailingOnly=TRUE)

#genomename <- ifelse(!( is.na(args[1])), args[1], "year1pat")
genomename <- ifelse(!( is.na(args[1])), args[1], "v2_trio.hap2")
outputdir <- ifelse(!( is.na(args[2])), args[2], ".")
resourcedir <- ifelse(! ( is.na(args[3])), args[3], ".")
#plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c(genomename, " Errors vs. hg002v1.0.1"), sep="", collapse=""))
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c("Verkko2.0 Trio Assembly Haplotype 2", " Errors vs. hg002v1.0.1"), sep="", collapse=""))
#plottitle <- c("Year 1 HPRC HG002 Paternal Assembly Alignment Coverage and Errors vs. hg002v1.0.1")
genomefile <- paste(c(outputdir, "/v1.0.1.genome.bed"), sep="", collapse="")
benchgenome <- toGRanges(genomefile)

maternalchroms <- paste0("chr", c(1:22, "X"), "_MATERNAL")
paternalchroms <- paste0("chr", c(1:22, "Y"), "_PATERNAL")
chroms <- c(maternalchroms, paternalchroms)

errorfile <- paste(c(outputdir, "/", genomename, ".errortype.v1.0.1.bed"), sep="", collapse="")
benchcoveredfile <- paste(c(outputdir, "/", genomename, ".benchcovered.v1.0.1.merged.bed"), sep="", collapse="")
benchcoveredpdf <- paste(c(outputdir, "/", genomename, ".benchcovered.v1.0.1.pdf"), sep="", collapse="")
#plotname <- paste(c(outputdir, "/", genomename, ".testcovered.v1.0.1.pdf"), sep="", collapse="")

benchcoveredranges <- toGRanges(benchcoveredfile)
benchcovereddf <- read.table(benchcoveredfile, header=FALSE, sep="\t")

nlocfile <- paste(c(resourcedir, "/v1.0.1.ncoords.bed"), sep="", collapse="")
nlocranges <- toGRanges(nlocfile)

errorsites <- read.table(errorfile, header=FALSE, sep="\t")
names(errorsites) <- c("chrom", "start", "end", "errorname", "score", "strand", "widestart", "wideend", "rgbcolor", "errortype", "assemblyerrorname")
#names(errorsites) <- c("chrom", "start", "end", "errorname", "errortype")

namefieldlengths <- sapply(seq(1, length(errorsites$chrom)), function(x) {length(strsplit(errorsites$errorname[x], split="_")[[1]])})
errorsites$ref <- sapply(seq(1,length(errorsites$chrom)), function(x) {strsplit(errorsites$errorname[x], split="_")[[1]][[namefieldlengths[x] - 1]]})
errorsites$alt <- sapply(seq(1,length(errorsites$chrom)), function(x) {strsplit(errorsites$errorname[x], split="_")[[1]][[namefieldlengths[x]]]})


plotheatmaps <- function(kp, phasingerrors, snperrors, indelerrors) {
  chromlengths <- width(ranges(benchgenome))
  names(chromlengths) <- seqnames(benchgenome)
  windowranges <- tileGenome(seqlengths=chromlengths, tilewidth=1000000)
  phasingerrorcounts <- sapply(seq(1, length(windowranges)), function(x) {sum(phasingerrorgranges %over% windowranges[x])})
  kpHeatmap(kp, data=windowgranges, y=phasingerrorcounts, r0=0.1, r1=0.3, colors=c("lightblue", "blue", "black"))
  snperrorcounts <- sapply(seq(1, length(windowranges)), function(x) {sum(snperrorgranges %over% windowranges[x])})
  kpHeatmap(kp, data=windowgranges, y=snperrorcounts, data.panel=1, r0=0.4, r1=0.6, colors=c("lightgreen", "darkgreen", "black"))
  indelerrorcounts <- sapply(seq(1, length(windowranges)), function(x) {sum(indelerrorgranges %over% windowranges[x])})
  kpHeatmap(kp, data=windowgranges, y=indelerrorcounts, data.panel=1, r0=0.7, r1=0.9, colors=c("white", "pink", "red"))
  
}

plotlinecurves <- function(kp, phasingerrorgranges, snperrorgranges, indelerrorgranges) {
  kp <- kpPlotDensity(kp, data=phasingerrorgranges, data.panel=1, r0=0.1, r1=0.9, border="blue", col=NA, window.size=100000, fill=NULL)
  maxphasingerrordensity <- 
  kp <- kpPlotDensity(kp, data=snperrorgranges, data.panel=1, r0=0.1, r1=0.9, col=NA, border="red", window.size=100000)
  kp <- kpPlotDensity(kp, data=indelerrorgranges, data.panel=1, r0=0.1, r1=0.9, col=NA, border="darkgreen", window.size=100000)
}

plotrectangles <- function(kp, phasingerrorgranges, snperrorgranges, indelerrorgranges) {
  kpRect(kp, data=phasingerrorgranges, border="blue", lwd=0.1, col="blue", data.panel=1, y0=rep(0, length(phasingerrorgranges)), y1=rep(0.2, length(phasingerrorgranges)))
  # then snps:
  kpRect(kp, data=snperrorgranges, border="red", lwd=0.1, col="red", data.panel=1, y0=rep(0.2, length(snperrorgranges)), y1=rep(0.4, length(snperrorgranges)))
  # then indels:
  kpRect(kp, data=indelerrorgranges, border="darkgreen", lwd=0.1, col="darkgreen", data.panel=1, y0=rep(0.4, length(indelerrorgranges)), y1=rep(0.6, length(indelerrorgranges)))
}

# calculate error profiles and make plot:

indelerrors <- errorsites[errorsites$errortype=="CONSENSUS" & (nchar(errorsites$alt)>1 | nchar(errorsites$ref)>1 | errorsites$alt=="*" | errorsites$ref=="*"),]
indelerrorgranges <- as(indelerrors, "GRanges")
snperrors <- errorsites[errorsites$errortype=="CONSENSUS" & nchar(errorsites$alt)==1 & nchar(errorsites$ref)==1 & errorsites$alt!="*" & errorsites$ref!="*",]
snperrorgranges <- as(snperrors, "GRanges")
phasingerrors <- errorsites[errorsites$errortype=="PHASING",]
phasingerrorgranges <- as(phasingerrors, "GRanges")

pdf(benchcoveredpdf, 8.5, 11.0)
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 400
pp$data1height <- 1200
#pp$data2height <- 1000
pp$topmargin <- 650

bordercol = 'black'
borderlwd = 1.0

aligncolor <- ifelse(str_detect(benchcovereddf$V1, "MAT"), "green", "blue")
kp <- plotKaryotype(genome=benchgenome, plot.type=1, chromosomes=chroms, main=plottitle, cex=0.4, plot.params=pp)
kpAddBaseNumbers(kp, tick.dist = 20000000, tick.len = 10, cex=0.3)
kpRect(kp, data=benchcoveredranges, border=bordercol, lwd=borderlwd, col=aligncolor, data.panel="ideogram", y0=rep(0, length(benchcoveredranges)), y1=rep(1, length(benchcoveredranges)))

kpRect(kp, data=nlocranges, border=bordercol, lwd=borderlwd, col="red", data.panel="ideogram", y0=rep(0, length(nlocranges)), y1=rep(1, length(nlocranges)))

plotlinecurves(kp, phasingerrorgranges, snperrorgranges, indelerrorgranges)
legend("bottomright", c("Phasing", "SingleNuc", "Indel"), col=c("blue", "red", "darkgreen"), lty=c(1,1,1,1), title="Error Type")
legend("topright", inset = c(0,  0.075), c("Paternal", "Maternal", "N's (rDNAs)", "Uncovered"), fill=c("blue", "green", "red", "grey"), title="Aligned Regions")

dev.off()

